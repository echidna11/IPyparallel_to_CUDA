#!/usr/bin/env python

'''
load a volume in npy format, then save it to mrc format
'''

import os
import sys
import numpy as N

import tomominer.io.file as IF

if __name__ == '__main__':
    npy_file = sys.argv[1]
    mrc_file = sys.argv[2]

    print 'loading', npy_file
    v = N.load(npy_file)

    print 'writing', mrc_file
    IF.put_mrc(v, mrc_file, overwrite=True)

