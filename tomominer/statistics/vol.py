
# statistics of volumes

import scipy.stats as ST

import numpy as N
import numpy.fft as NF

import tomominer.image.vol.util as GV


# signal to noise ratio given two volumes representing two realizations
def snr(v1, v2):

    pr = ST.pearsonr(v1.flatten(), v2.flatten())
    cor = pr[0]
    
    snr = cor / (1 - cor)
    
    return snr
    


# Fourier Shell correlation between two volumes
# translated from    Resolution.fsc()
def fsc(v1, v2, band_width_radius=1.0):
    
    siz = v1.shape
    assert(siz == v2.shape)

    origin_co = GV.fft_mid_co(siz)
    
    x = N.mgrid[0:siz[0], 0:siz[1], 0:siz[2]]
    for dim_i in range(3):      x[dim_i] -= origin_co[dim_i]

    rad = N.sqrt( N.square(x).sum(axis=0) )

    vol_rad = int( N.floor( N.min(siz) / 2.0 ) + 1)

    v1f = NF.fftshift( NF.fftn(v1) )
    v2f = NF.fftshift( NF.fftn(v2) )

    fsc_cors = N.zeros(vol_rad)

    # the interpolation can also be performed using scipy.ndimage.interpolation.map_coordinates()
    for r in range(vol_rad):

        ind = ( abs(rad - r) <= band_width_radius )

        c1 = v1f[ind]
        c2 = v2f[ind]

        fsc_cor_t = N.sum( c1 * N.conj(c2) ) / N.sqrt( N.sum( N.abs(c1)**2 ) * N.sum( N.abs(c2)**2) )
        fsc_cors[r] = N.real( fsc_cor_t )

    return fsc_cors


'''
# test commands



import sys
sys.path.append('/home/rcf-47/mxu/ln/tomominer/tomominer')
sys.path.append('/home/rcf-47/mxu/ln/tomominer/build/lib.linux-x86_64-2.7/tomominer')

import tomominer.model.util as MU
v1 = MU.generate_toy_model()

import numpy as N
v2 = v1 + N.random.random(v1.shape) * (1 * v1.std())

import tomominer.statistics.vol as SV
fsc_t = SV.fsc(v1, v2)

print fsc_t

resolution_t = SV.fft_resolution(min(v1.shape), 0.5)
print resolution_t
print resolution_t[0].shape



'''

def band_to_resolution(voxel_spacing, band_num, band):
    frequency = float(band)            # band: The current resolution in pixels (python convention, i.e., start from 0)
    return      (2.0  * voxel_spacing * band_num) /  (frequency + 1)


# this function calculates resolution from bands
# it can handle both case that band number and voxel spacing are same or different for different dimensions
# see bandToAngstrom() in pytom/basic/resolution.py
def fft_resolution_list(band_num, voxel_spacing):
    if type(band_num) == int:       band_num = N.array([band_num])

    if type(voxel_spacing) == float:        voxel_spacing = N.ones(len(band_num)) * voxel_spacing

    assert len(voxel_spacing) == len(band_num)

    re = {}
    for dim_i in range(len(band_num)):
        bands =         N.array( range(band_num[dim_i]) )
        re[dim_i] = [band_to_resolution(voxel_spacing=voxel_spacing[dim_i], band_num=band_num[dim_i], band=_) for _ in bands]

    return re



'''

import sys
sys.path.append('/home/rcf-47/mxu/ln/tomominer/tomominer')


import tomominer.statistics.vol as SV

print SV.fft_resolution_list(50, 1.2)



'''






def resolution(v1, v2, voxel_spacing, cutoff=0.5):

    n = v1.shape[0]
    assert  N.all( v1.shape == (n, n, n) )

    c = fsc(v1, v2)

    return fsc_resolution(c, voxel_spacing=voxel_spacing, cutoff=cutoff)


def fsc_resolution(fsc, voxel_spacing, cutoff):
    bi = band_interpolation(fsc, cutoff=cutoff)
    return  band_to_resolution(voxel_spacing=voxel_spacing, band_num=len(fsc), band=bi)

# from a list of FSCs of different bands, find the interpolated band number given a cutoff
def band_interpolation(c, cutoff=0.5):

    if c[0] < cutoff:       return 0
    if N.min(c) >= cutoff:      return len(c) - 1

    for i in range(len(c) -1):
        if (c[i] >= cutoff) and (c[i+1] < cutoff):            break

    assert      (c[i] >= cutoff) and (c[i+1] < cutoff)

    c0 = c[i] - cutoff
    c1 = c[i+1] - cutoff

    return      i + (c0 / (c0 - c1))



