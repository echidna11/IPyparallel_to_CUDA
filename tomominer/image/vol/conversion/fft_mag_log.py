#!/usr/bin/env python

'''
FFT magnitude of a image
'''

import os
import sys
import numpy as N
from numpy.fft import fftn, ifftn, fftshift, ifftshift


import tomominer.io.file as IF

if __name__ == '__main__':

    in_file = sys.argv[1]
    out_file = sys.argv[2]

    if os.path.isfile(out_file):
        print 'output exists, delete it and try again'
        sys.exit(0)

    v = IF.read_mrc_vol(in_file)
    v = fftn(v)
    v = N.log(N.abs(fftshift(v)))

    IF.put_mrc(v, out_file)


