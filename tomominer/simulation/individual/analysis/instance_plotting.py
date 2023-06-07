#!/usr/bin/env python

# given json file, load its subtomograms and save the center slices
# to combined slices into a single figure, use montage, for example     montage -trim -tile 6x50 -geometry +40+0 g/slice-g-* /tmp/combined.png

import os
import json
import pickle

import numpy as N

import matplotlib.pyplot as PLT
import matplotlib.cm as MCM

import pylab

import tomominer.io.file as IF


def save_slices(v, out_file, directions=[0,1,2]):
    assert len(directions) > 0

    s = [None] * 3

    i = v.shape[0]/2
    s[0] = N.squeeze(v[i,:,:])

    i = v.shape[1]/2
    s[1] = N.squeeze(v[:,i,:])

    i = v.shape[2]/2
    s[2] = N.squeeze(v[:,:,i])


    f, sp = PLT.subplots(nrows=1, ncols=len(directions))

    if len(directions) == 1:
        # in such case, sp is a single matplotlib.axes.AxesSubplot object
        assert  type(sp) != list
        sp = [sp]

    for i in range(len(directions)):
        sp[i].imshow( s[directions[i]], cmap = MCM.Greys_r )

    for i in range(len(directions)):
        sp[i].axis('off') # clear x- and y-axes

    #PLT.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)       # reduce margins
    PLT.tight_layout()      # reduce margins

    PLT.draw()

    PLT.savefig(out_file, bbox_inches='tight')
    
    PLT.close("all")





if __name__ == '__main__':
    
    with open('instance_plotting__config.json') as f:     op = json.load(f)

    with open(op['data file']) as f:    dj = json.load(f)

    if not os.path.isdir(op['out dir']):        os.makedirs(op['out dir'])


    for i, r in enumerate(dj):
        print i

        v = IF.get_mrc(r['subtomogram'])

        save_slices(v=v, out_file=os.path.join(op['out dir'], 'slice-%05d.png'%(r['id'] if ('id' in r) else i,)), directions=op['directions'])

