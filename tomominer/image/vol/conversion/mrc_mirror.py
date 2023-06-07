#!/usr/bin/env python

'''
generate three mirrors of one mrc volume file

~/ln/tomominer/tomominer/image/vol/conversion/mrc_mirror.py
'''

import os, sys, json
import numpy as N

import tomominer.io.file as IF

if __name__ == '__main__':


    mrc = sys.argv[1]

    v = IF.read_mrc_vol(mrc)

    mrc_dir = os.path.dirname(mrc)
    mrc_base = os.path.basename(mrc)
    mrc_root, mrc_ext = os.path.splitext(mrc_base)

    for d in range(3):
        if d == 0:
            vt = v[::-1, :, :]
        elif d == 1:
            vt = v[:, ::-1, :]
        elif d == 2:
            vt = v[:, :, ::-1]
        else:
            assert False

        mrc_out = os.path.join(mrc_dir, mrc_root + '--mirror-%d'%(d,) + mrc_ext)
        assert  not os.path.isfile(mrc_out)

        print 'writing', mrc_out
        IF.put_mrc(vt, mrc_out, overwrite=True)

