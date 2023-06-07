#!/usr/bin/env python


'''
simply fill up one volume's nan values with zeros
'''

import os
import sys
import numpy as N

import tomominer.io.file as IF

if __name__ == '__main__':
    file_in = sys.argv[1]
    file_out = sys.argv[2]

    print 'loading', file_in
    v = IF.read_mrc_vol(file_in)

    v[N.isnan(v)] = 0

    print 'writing', file_out
    IF.put_mrc(v, file_out, overwrite=True)

