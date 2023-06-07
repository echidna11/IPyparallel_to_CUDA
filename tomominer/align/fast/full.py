
import numpy as N
from numpy.fft import fftn, ifftn, fftshift, ifftshift

import tomominer.geometry.ang_loc as AL
import tomominer.image.vol.util as GV
import tomominer.align.fast.util as AFU

import tomominer.core as tomo

# fast alignment without missing wedge

def align_vols_no_wedge(v1, v2, L=36, peak_spacing=8, normalize=True, true_ang_loc=None):

    v1m = (v1 - v1.mean());       del v1
    v2m = (v2 - v2.mean());       del v2

    radii = 1.0 + N.array(range(int( N.floor(max(v1m.shape) / 2.0)) ), dtype=N.float)

    v1fa = N.abs(fftshift( fftn(v1m) ));     v1fa = N.array(v1fa, order='F')
    v2fa = N.abs(fftshift( fftn(v2m) ));     v2fa = N.array(v2fa, order='F')
    cor = tomo.adp_em_cor(v1fa, v2fa, L,  radii);        cor = N.real(cor);      cor = N.array(cor, order='F')

    angs = tomo.local_max_angles_py(cor, peak_spacing)

    best = {}
    best['cor'] = (-N.inf)

    for ang_i, (ang, score) in enumerate(angs):
        ang = N.array(ang)      # convert from a tuple to array
        v2mr = tomo.rotate_vol_pad_mean_py(v2m, ang, N.array([0.0, 0.0, 0.0]))
        loc, cor_max = AFU.translation_align(v1m, v2mr)

        if cor_max > best['cor']:
            if true_ang_loc is not None:
                print cor_max, AL.angle_difference_zyz(ang, true_ang_loc['angle']), AL.angle_zyz_translation_difference(ang1=ang, loc1_r=loc, ang2=true_ang_loc['angle'], loc2_r=true_ang_loc['loc'])

            best['cor'] = cor_max
            best['angle'] = ang
            best['loc'] = loc


    if normalize:
        best['cor'] = (best['cor'] / v1m.size) / (N.std(v1m) * N.std(v2m))        # normalize to Pearson correlation, this is used to facilate comparisons with different alignments

    return best


"""
# test code


import sys
sys.path.append('/home/rcf-47/mxu/ln/tomo/py')

import tomominer.model.util as mu
v2 = mu.generate_toy_model(dim_siz=64)




import numpy as N
loc = N.round( (N.random.random(3) - 0.5) * (N.array(v2.shape) * 0.5) )
if False:
    ang = N.round( 36 * (N.random.random(3) * N.pi * 2.0) ) / 36
else:
    ang = (N.random.random(3) * N.pi * 2.0) 


import tomominer.core as tomo
v1 = tomo.rotate_vol_pad_zero_py(N.array(v2, order='F'), ang, loc)

import tomominer.align.fast_full as aff
a_re = aff.align_vols_no_wedge(v1=v1, v2=v2, L=36, peak_spacing=4, true_ang_loc={'angle':ang, 'loc':loc})

print a_re



print ang, a_re['angle']
print loc, a_re['loc'], loc - a_re['loc']


import geometery.ang_loc as AL
print AL.angle_zyz_translation_difference(ang1=ang, loc1_r=loc, ang2=a_re['angle'], loc2_r=a_re['loc'])


import pdb; pdb.set_trace()

"""




