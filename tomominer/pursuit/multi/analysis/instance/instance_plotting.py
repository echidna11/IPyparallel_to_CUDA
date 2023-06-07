#!/usr/bin/env python

"""
given a cluster and its average, order these instances according to alignment scores, align all its instances to the average
then save the center slices, and isosurfaces of individual instances
to combined slices into a single figure, use montage, for example     montage -trim -tile 6x50 -geometry +40+0 g/slice-g-* /tmp/combined.png

~/ln/tomominer/tomominer/pursuit/multi/analysis/instance/instance_plotting.py
"""


import os, sys, json,  pickle

import numpy as N

import matplotlib.pyplot as PLT
import matplotlib.cm as MCM

import pylab

import tomominer.io.file as IF
import tomominer.geometry.rotate as GR
import tomominer.filter.gaussian as FG
import tomominer.geometry.ang_loc as GA


def save_slices(v, out_file):
    i = v.shape[0]/2
    s0 = N.squeeze(v[i,:,:])

    i = v.shape[1]/2
    s1 = N.squeeze(v[:,i,:])

    i = v.shape[2]/2
    s2 = N.squeeze(v[:,:,i])


    f, sp = PLT.subplots(nrows=1, ncols=3)


    sp[0].imshow( s0, cmap = MCM.Greys_r )
    sp[1].imshow( s1, cmap = MCM.Greys_r )
    sp[2].imshow( s2, cmap = MCM.Greys_r )

    for i in range(sp.shape[0]):
            sp[i].axis('off') # clear x- and y-axes

    PLT.draw()

    PLT.savefig(out_file, bbox_inches='tight')
    
    PLT.close("all")



def main():
    
    with open('instance_plotting__op.json') as f:     op = json.load(f)

    if ('template rotate file' in op) and os.path.isfile(op['template rotate file']):
        with open(op['template rotate file']) as f:     tr = json.load(f)
    else:
        tr = None

    if not os.path.isdir(op['out dir']):        os.makedirs(op['out dir'])

    v = IF.read_mrc_vol(op['template file'])
    if tr is not None:      v = GR.rotate_pad_mean(v, angle=tr['angle'], loc_r=tr['loc'])
    save_slices(v, os.path.join(op['out dir'], 'template.png'))

    with open(op['data file']) as f:       dj = json.load(f)
    dj = sorted(dj, key=lambda _ : _['score'], reverse=True)

    if ('top num' in op) and (len(dj) > op['top num']):     dj = dj[:op['top num']]
    
    for i, r in enumerate(dj):
        print '\r%05d'%(i,),            ;       sys.stdout.flush()

        v = IF.read_mrc_vol(r['subtomogram'])

        if tr is not None:
            (rm3, loc3) = GA.combined_transform(rm1=GA.rotation_matrix_zyz(r['angle']), loc_r1=N.array(r['loc']), rm2=GA.rotation_matrix_zyz(tr['angle']), loc_r2=N.array(tr['loc']) )
            vr = GR.rotate_pad_mean(v, rm=rm3, loc_r=loc3)
        else:
            vr = GR.rotate_pad_mean(v, angle=r['angle'], loc_r=r['loc'])

        save_slices(vr, os.path.join(op['out dir'], 'slice-%05d.png'%(i,)))

        if 'smooth sigma' in op:
            vrg = FG.smooth(vr, sigma=op['smooth sigma'])
            save_slices(vrg, os.path.join(op['out dir'], 'slice-g-%05d.png'%(i,)))


if __name__ == '__main__':
    main()


