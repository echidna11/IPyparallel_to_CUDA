
# functions for gaussian filtering


import gc as GC

import numpy as N
from numpy.fft import fftn, ifftn, fftshift, ifftshift

import tomominer.image.vol.util as CV
import tomominer.model.util as MU

import scipy.ndimage as SN
import scipy.ndimage.filters as SNF

# 3D gaussian filtering of a volume (v)

# smoothing using scipy.ndimage.gaussian_filter
def smooth(v, sigma):
    assert  sigma > 0
    return SN.gaussian_filter(input=v, sigma=sigma)

'''
# check whether scipy.ndimage.gaussian_filter is 3D

import numpy as N
v = N.zeros(    (30,30,30)      )

# add three central slices. If the smoothing is 3D but not 2D, we should see neighbor slices gets intensities.
v[5:25, 5:25, 15] = 1.0
v[5:25, 15, 5:25] = 1.0
v[15, 5:25, 5:25] = 1.0

import scipy.ndimage as SN
v = SN.gaussian_filter(input=v, sigma=2)

import tomominer.io.file as IF
IF.put_mrc(v, '/tmp/v.mrc')

'''


# smooth using fft, reflect padding
def smooth_fft(v, sigma):

    # generate a larger padded volume, in order to avoid the boundary effect in covolution
    pad_width = int(N.round(sigma*2))
    vp = N.pad(array=v, pad_width=pad_width, mode='reflect')

    g = MU.gauss_function(size=vp.shape, sigma=sigma)

    g_fft = fftn(ifftshift(g));     # use ifftshift(g) to move center of gaussian to origin

    del g

    v_conv = N.real( ifftn( fftn(vp) * N.conj( g_fft ) ) )

    v_conv = v_conv[pad_width:(pad_width+v.shape[0]), pad_width:(pad_width+v.shape[1]), pad_width:(pad_width+v.shape[2])]
    assert      v.shape == v_conv.shape

    return v_conv

# using convolution with mirroring
def smooth__convolve(v, sigma):
    raise   Exception('somehow this does not work well')

    g = MU.gauss_function(size=v.shape, sigma=sigma)

    mid_co = (N.array(g.shape) - 1) / 2          # IMPORTANT: following python convension, in index starts from 0 to size-1!!! So (siz-1)/2 is real symmetry center of the volume
    return SNF.convolve(input=v, weights=g, mode='reflect')



'''

# test scripts

# put a single 1 randomly in any position in the map, see after gaussian smoothing, the highest value is still at same position



import os,sys
sys.path.append(os.path.join( os.getenv('HOME'), 'ln/tomominer/tomominer' ))
sys.path.append(os.path.join( os.getenv('HOME'), 'ln/tomominer/build/lib.linux-x86_64-2.7/tomominer' ))




import numpy as N
siz = [30+N.random.randint(10), 40+N.random.randint(10), 50+N.random.randint(10)]
v = N.zeros(siz)


x = [  N.random.randint(siz[0]),  N.random.randint(siz[1]), N.random.randint(siz[2])  ]

v[x[0], x[1], x[2]] = 1.0


import tomominer.filter.gaussian as FG
if True:
    vg = FG.smooth(v, sigma = 2.0)
else:
    vg = FG.dog_smooth__large_map(v, s1=2.0)

print x, N.where(vg == vg.max())





# construct a toy unsymmetric structure, then see if the smoothed map has same orientation



v = N.zeros(siz)
s0=10;  s1=30

v[s0:s1, s0, s0] = 1.0
v[s0, s0:s1, s0] = 1.0
v[s0, s0, s0:s1] = 1.0


vg = FG.smooth(v, sigma = 3.0)


import matplotlib
matplotlib.use('Qt4Agg')

import tomominer.image.vol.util as CV

CV.dsp_cub(v)
CV.dsp_cub(vg)

import matplotlib.pyplot as plt
plt.show()





'''


def dog_smooth(v, s1, s2=None):
    if s2 is None:      s2 = s1 * 1.1       # the 1.1 is according to a DoG particle picking paper
    assert      s1 < s2

    return  smooth(v, s1) - smooth(v, s2)



# convolute with a dog function, delete unused data when necessary in order to save memory for large maps
def dog_smooth__large_map(v, s1, s2=None):

    if s2 is None:      s2 = s1 * 1.1       # the 1.1 is according to a DoG particle picking paper
    assert      s1 < s2

    size = v.shape


    pad_width = int(N.round(s2*2))
    vp = N.pad(array=v, pad_width=pad_width, mode='reflect')

    v_fft = fftn(vp).astype(N.complex64)
    del v;      GC.collect()


    g_small = MU.difference_of_gauss_function(size=N.array([int(N.round(s2 * 4))]*3), sigma1=s1, sigma2=s2)
    assert      N.all(N.array(g_small.shape) <= N.array(vp.shape))       # make sure we can use CV.paste_to_whole_map()

    g = N.zeros(vp.shape)
    CV.paste_to_whole_map(whole_map=g, vol=g_small, c=None)

    g_fft_conj = N.conj(   fftn(ifftshift(g)).astype(N.complex64)   )    # use ifftshift(g) to move center of gaussian to origin
    del g;      GC.collect()

    prod_t = (v_fft * g_fft_conj).astype(N.complex64)
    del v_fft;      GC.collect()
    del g_fft_conj;      GC.collect()

    prod_t_ifft = ifftn( prod_t ).astype(N.complex64)
    del prod_t;      GC.collect()

    v_conv = N.real( prod_t_ifft )
    del prod_t_ifft;      GC.collect()
    v_conv = v_conv.astype(N.float32)

    v_conv = v_conv[(pad_width+1):(pad_width+size[0]+1), (pad_width+1):(pad_width+size[1]+1), (pad_width+1):(pad_width+size[2]+1)]
    assert      size == v_conv.shape

    return v_conv







