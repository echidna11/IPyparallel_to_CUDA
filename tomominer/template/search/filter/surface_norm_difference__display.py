#!/usr/bin/env python



'''

show the vectors of surface normal and normal of template, using mayavi


~/ln/tomominer/tomominer/template/search/filter/surface_norm_difference__display.py

'''


import json

import numpy as N
from mayavi import mlab


def main():

    with open('surface_norm_difference__display__op.json') as f:        op = json.load(f)


    with open(op['data file']) as f:        dj = json.load(f)


    x = N.array([_['peak']['loc'] for _ in dj])
    n0 = N.array([_['surface_norm']['n0'] for _ in dj])
    n = N.array([_['surface_norm']['n'] for _ in dj])


    scale_factor = op['scale factor']

    mlab.figure(bgcolor=(0,0,0))
    if op['plot surface norm']:           mlab.quiver3d(x[:,0], x[:,1], x[:,2], n[:,0], n[:,1], n[:,2], scale_factor=scale_factor, color=(1,0,0))      # plot surface norm
    if op['plot template norm']:             mlab.quiver3d(x[:,0], x[:,1], x[:,2], n0[:,0], n0[:,1], n0[:,2], scale_factor=scale_factor, color=(0,0,1))      # plot template norm

    mlab.show()



if __name__ == '__main__':
    main()

