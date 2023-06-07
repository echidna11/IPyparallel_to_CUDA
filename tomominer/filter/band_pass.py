

import numpy as N
from numpy.fft import fftn, ifftn, fftshift, ifftshift

import tomominer.image.vol.util as GV


'''
band pass filtering from given curve
'''


def filter_given_curve(v, curve):

    grid = GV.grid_displacement_to_center(v.shape, GV.fft_mid_co(v.shape))
    rad = GV.grid_distance_to_center(grid)
    rad = N.round(rad).astype(N.int)

    b = N.zeros(rad.shape)
    for i, a in enumerate(curve):        b[rad == i] = a

    vf = ifftn(ifftshift(      fftshift(fftn(v)) * b       ))
    vf = N.real(vf)

    return vf

