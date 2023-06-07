#!/usr/bin/env python

'''
given a data json file, randomly select instances, then plot xz planes of center slices

to combined slices into a single figure, use montage, for example:
montage -trim -tile 14x3 -geometry +10+10 /tmp/instances/slice-* ~/tmp2/combined.png
'''

import os
import sys
import json
import random 

import numpy as N

import matplotlib.pyplot as PLT
import matplotlib.cm as MCM

import tomominer.io.file as IF

def save_slices(v, out_file):
    i = v.shape[1]/2
    s1 = N.squeeze(v[:,i,:])

    f, sp = PLT.subplots(nrows=1, ncols=1)


    sp.imshow( s1, cmap = MCM.Greys_r )
    sp.axis('off') # clear x- and y-axes


    PLT.draw()

    PLT.savefig(out_file, bbox_inches='tight')
    
    PLT.close("all")



if __name__ == '__main__':
    
    with open('instance_plotting__random__config.json') as f:     op = json.load(f)

    with open(op['data_json_file']) as f:   dj = json.load(f)
    dj = random.sample(dj, op['sample_num'])



    for i, r in enumerate(dj):
        print '\r', i, '           ',           ;           sys.stdout.flush()

        v = IF.get_mrc(r['subtomogram'])
        if ('inverse_intensity' in op) and op['inverse_intensity']:     v = -v

        save_slices(v, os.path.join(op['out_dir'], 'slice-%05d.png'%(i,)))

