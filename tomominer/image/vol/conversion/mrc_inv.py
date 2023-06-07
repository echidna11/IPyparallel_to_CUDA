#!/usr/bin/env python

'''
load one or several mrc files, invert its image intensities
'''

import os
import sys
import numpy as N

import tomominer.io.file as IF

if __name__ == '__main__':

    for i, mrc in enumerate(sys.argv):
        if i == 0:  continue

        v = IF.read_mrc_vol(mrc)

        mrc_dir = os.path.dirname(mrc)
        mrc_base = os.path.basename(mrc)
        mrc_root, mrc_ext = os.path.splitext(mrc)
        
        mrc_inv = mrc_root + '-inv' + mrc_ext
        if os.path.isfile(mrc_inv):
            print 'ignoring', mrc_inv
            continue

        print 'writing', mrc_inv
        IF.put_mrc(-v, mrc_inv, overwrite=True)

