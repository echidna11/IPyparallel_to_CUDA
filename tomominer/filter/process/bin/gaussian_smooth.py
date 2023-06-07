#!/usr/bin/env python

'''
command line program for gaussian smoothing

~/ln/tomominer/tomominer/filter/process/bin/gaussian_smooth.py
'''

import os, json
import cPickle as pickle

import tomominer.io.file as IF
import tomominer.filter.gaussian as FG


def main():
    with open('gaussian_smooth__op.json') as f:     op = json.load(f)


    s = float(op['sigma'])
    if 'voxel spacing' in op:
        # convert sigma to the unit of voxels
        s /= float(op['voxel spacing'])

    assert not os.path.isfile(op['mrc out'])

    print 'reading', op['mrc in']
    v = IF.read_mrc_vol(op['mrc in'])
    print 'smoothing'
    vg = FG.smooth(v, sigma=s)

    if ('reverse intensity' in op) and op['reverse intensity']:
        # sometimes, the gaussian smoothing will result in a reverse of intensity, dont know why, in such case, we reverse it back
        print 'reverse intensity of smoothed map'
        vg = -vg

    assert not os.path.isfile(op['mrc out'])
    print 'writing', op['mrc out']
    if ('writing method' in op) and (op['writing method'] == 'chunk'):
        IF.write_mrc_by_chunk(vg, op['mrc out'])
    else:
        IF.put_mrc(vg, op['mrc out'], overwrite=True)


if __name__ == '__main__':
    main()


