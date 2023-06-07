
'''
~/ln/tomominer/tomominer/align/fast/util.py
'''



from numpy.fft import fftn, ifftn, fftshift, ifftshift
import numpy as N 

import tomominer.geometry.rotate as GR
import tomominer.model.util as MU
import tomominer.image.vol.util as IVU

import tomominer.core as core


# fast alignment according to JSB 2012 paper
# given two subtomograms and their masks, perform populate all candidate rotational angles, 
# with missing wedge correction

# according to SphericalHarmonicsUtil.combined_search(), opt.search_method=8

def fast_roatation_align(v1, m1, v2, m2, max_l):
    raise Exception('core wrapper not implementated yet')

    radius = int( max(v1.shape) / 2 )

    radii = list(range(1,radius+1))     # radii must start from 1, not 0!

    # fftshift breaks order='F'
    v1fa = abs(fftshift(fftn(v1)))
    v2fa = abs(fftshift(fftn(v2)))

    v1fa = v1fa.copy(order='F')
    v2fa = v2fa.copy(order='F')

    m1sq = N.square(m1)
    m2sq = N.square(m2)

    a1t = v1fa * m1sq
    a2t = v2fa * m2sq

    cor12 = core.adp_em_cor(a1t, a2t, max_l, radii)
    angs = core.local_max_angles(N.real(cor12).copy(order='F'), 8)


    sqt_cor11 = N.sqrt( N.real( core.adp_em_cor( N.square(v1fa) * m1sq, m2sq, max_l, radii ) ) )
    sqt_cor22 = N.sqrt( N.real( core.adp_em_cor( m1sq, N.square(v2fa) * m2sq, max_l, radii ) ) )

    cors = cor12 / (sqt_cor11 * sqt_cor22)

    # N.real breaks order='F' by not making explicit copy.
    cors = N.real(cors)
    cors = cors.copy(order='F')

    angs = core.local_max_angles(cors, 8)

    return angs



'''
# test code

v1f = '/panfs/cmb-panasas2/mxu/proj/imaging/electron/frequent_structure/data/out/method/classification-system/pytom/thermosomes/particle/test_SNR1wedge30/test_SNR1wedge30_100.mrc'
m1f = '/panfs/cmb-panasas2/mxu/proj/imaging/electron/frequent_structure/data/out/method/classification-system/pytom/thermosomes/particle/test_SNR1wedge30/wedge-mask.mrc'
v2f = '/panfs/cmb-panasas2/mxu/proj/imaging/electron/frequent_structure/data/out/method/classification-system/pytom/thermosomes/particle/test_SNR1wedge30/test_SNR1wedge30_101.mrc'
m2f = '/panfs/cmb-panasas2/mxu/proj/imaging/electron/frequent_structure/data/out/method/classification-system/pytom/thermosomes/particle/test_SNR1wedge30/wedge-mask.mrc'

import tomominer.core as core

v1 = core.parse_mrc(v1f)
m1 = core.parse_mrc(m1f)

if False:
    v2 = core.parse_mrc(v2f)
    m2 = tomo.parse_mrc(m2f)
else:
    v2 = v1
    m2 = m1

fast_rotation_align({'v1':v1, 'm1':m1, 'v2':v2, 'm2':m2, 'max_l':36})

'''


'''
assume v2, m2 are already rotated, calculate translation alignment
'''
def translation_align_with_mask(v1, m1, v2, m2):
    m = m1 * m2

    v1f = fftn(v1)      ;       v1f[0,0,0] = 0.0        ;           v1f = fftshift(v1f)     ;       v1f *= m      ;        v1f /= N.sqrt(N.square(N.abs(v1f)).sum())
    v2f = fftn(v2)      ;       v2f[0,0,0] = 0.0        ;           v2f = fftshift(v2f)     ;       v2f *= m      ;        v2f /= N.sqrt(N.square(N.abs(v2f)).sum())

    return translation_align__given_unshifted_fft(v1f=ifftshift(v1f), v2f=ifftshift(v2f))



def translation_align(v1, v2):

    v1f = fftn(v1)
    v2f = fftn(v2)

    return translation_align__given_unshifted_fft(v1f, v2f)


def translation_align__given_unshifted_fft(v1f, v2f):
    cor = fftshift( N.real( ifftn( v1f * N.conj(v2f) ) ) )
    
    mid_co = IVU.fft_mid_co(cor.shape)
    loc = N.unravel_index( cor.argmax(), cor.shape )

    return {'loc': (loc - mid_co + 1), 'cor': cor[loc[0], loc[1], loc[2]]}



"""
# test commands


import sys
sys.path.append('/home/rcf-47/mxu/ln/tomo/py')

import tomominer.model.util as mu
v2 = mu.generate_toy_model(dim_siz=32)




import numpy as N
loc = N.round( (N.random.random(3) - 0.5) * (N.array(v2.shape) * 0.4) )

import tomo
v1 = tomo.rotate_vol_pad_zero_py(N.array(v2, order='F'), N.array([0.0, 0.0, 0.0]), loc)

import tomominer.align.fast_full as aff
loc_a, cor = aff.translation_align(v1, v2)

print loc
print loc_a
print loc - loc_a
print cor



"""



# for each angle, do a translation search
def translation_align_given_rotation_angles(v1, m1, v2, m2, angs):

    raise Exception('wrong code')       # for each angle, need to first rotate m2 to get m2r, then use joint mask m = m1 * m2r to mask out values, then perform translation align

    # first mask v1 and v2 using masks, then rotate v2 and find translation, also calculate Forster score
    v1f = fftn(v1)      ;       v1f[0,0,0] = 0.0        ;           v1f = fftshift(v1f)     ;       v1f *= m1      ;        v1f /= N.sqrt(N.square(N.abs(v1f)).sum())       ;       v1fi = N.real(ifftn(ifftshift(v1f)))
    v2f = fftn(v2)      ;       v2f[0,0,0] = 0.0        ;           v2f = fftshift(v2f)     ;       v2f *= m2      ;        v2f /= N.sqrt(N.square(N.abs(v2f)).sum())       ;       v2fi = N.real(ifftn(ifftshift(v2f)))

    a = []
    for ang in angs:
        v2fir = GR.rotate_pad_mean(v2fi, ang=ang)
        lc = translation_align(v1fi, v2fir)

        a.append({'ang':ang, 'loc':lc['loc'], 'score':lc['cor']})

    return a


def fast_align(v1, m1, v2, m2, max_l=36):
    raise Exception('not tested yet')

    angs = fast_roatation_align(v1, m1, v2, m2, max_l)
    a = translation_align_given_rotation_angles(v1, m1, v2, m2, angs)

    a = sorted(a, key=lambda _:(-_['score']))

    return a

'''

# test code, generate two density maps, perform backprojection reconstruction, and perfrom fast alignment, then see how the alignment result compared to ground truth


import tomominer.align.fast.util as AFU
tp = AFU.generate_toy_subtomogram_pair()
als = AFU.fast_align(v1=tp['v1b'], m1=tp['m1'], v2=tp['v2b'], m2=tp['m2'], max_l=36)

for a in als:    print a['score'], GA.angle_zyz_translation_difference(ang1=tp['ang'], loc1_r=tp['loc'], ang2=a['ang'], loc2_r=a['loc'])


'''







'''
Fast real space rotation alignment according to Xu ISMB2013
Parameters: v1:subtomogram1, m1:mask1, v2:subtomogram2, m2:mask2, max_l:band width for rotation alignment

Requirement: both v1 and v2 are translated so that their rotation center is at center of subtomogram, and they both have zero mean
'''
def fast_real_space_rotation_align(v1, m1, v2, m2, max_l=36, constant_correction_method=None):

    # warning: fftshift breaks order='F'
    v1f = fftshift(fftn(v1)) * m1
    v2f = fftshift(fftn(v2)) * m2

    v1fi = N.real(ifftn(ifftshift(v1f)))
    v2fi = N.real(ifftn(ifftshift(v2f)))

    m1sq = N.square(m1)
    m2sq = N.square(m2)


    radius = int( max(v1.shape) / 2 )

    radii = list(range(1,radius+1))     # radii must start from 1, not 0!
    radii = N.array(radii, dtype=N.float)

    cor12 = core.rot_search_cor(v1=v1fi, v2=v2fi, radii=radii, L=max_l)

    sqt_cor11 = N.sqrt( core.rot_search_cor( v1=N.square(N.abs(v1f)), v2=m2sq, L=max_l, radii=radii ) )
    sqt_cor22 = N.sqrt( core.rot_search_cor( v1=m1sq, v2=N.square(N.abs(v2f)), L=max_l, radii=radii ) )

    cor_array = cor12 / (sqt_cor11 * sqt_cor22)

    if constant_correction_method is not None:
        if constant_correction_method == 0:
            cor_array *= cor_array.size
        elif constant_correction_method == 1:
            cor_array *= MU.sphere_mask(cor12.shape).sum()
        else:
            raise Exception('constant_correction_method')

    (cors, angs) = core.local_max_angles(cor=cor_array)
    i = N.argsort(-cors)
    cors = cors[i]
    angs = angs[i,:]

    ca = []
    for i in range(len(cors)):      ca.append({'score':cors[i], 'ang':angs[i,:].flatten()})

    return ca






'''

# test code, generate two density maps, perform backprojection reconstruction, and perfrom real space fast rotational alignment, then see how the alignment result compared to ground truth

import tomominer.align.fast.util as AFU
tp = AFU.generate_toy_subtomogram_pair(loc=[0,0,0], missing_wedge_angle=30)

tp['v1b'] -= tp['v1b'].mean()
tp['v2b'] -= tp['v2b'].mean()

import pickle
with open('/tmp/tp.pickle', 'wb') as f:     pickle.dump(tp, f)



import pickle
with open('/tmp/tp.pickle', 'rb') as f:     tp = pickle.load(f)

import tomominer.align.fast.util as AFU
als = AFU.fast_real_space_rotation_align(v1=tp['v1b'], m1=tp['m1'], v2=tp['v2b'], m2=tp['m2'], constant_correction_method=0)

import tomominer.geometry.ang_loc as GA
for a in als:    print a['score'], GA.angle_zyz_translation_difference(ang1=tp['ang'], ang2=a['ang']), '\t',


'''



def generate_toy_subtomogram_pair(ang=None, loc=None, dim_siz=64, missing_wedge_angle=30, snr=float('inf')):

    bp_op = {}
    bp_op['matlab_code_path'] = '/home/rcf-47/mxu/ln/frequent_structure/code'
    bp_op['matlab_tom_2008_code_path'] = '/home/rcf-47/mxu/proj/imaging/electron/matlab_tom/2008/TOM_Release_2008'
    bp_op['matlab_tom_2005_code_path'] = '/home/rcf-47/mxu/proj/imaging/electron/util/tom_2005'

    bp_op['model'] = {'missing_wedge_angle': missing_wedge_angle, 'titlt_angle_step': 1, 'SNR': snr}
    bp_op['ctf'] = {'pix_size': 1.0, 'Dz': -15.0, 'voltage': 300, 'Cs': 2.2, 'sigma': 0.4}



    import pymatlab
    session = pymatlab.session_factory(options='-nodesktop -nodisplay')
    session.run( 'addpath(\'%s\')'%(bp_op['matlab_code_path']) )
    session.run( 'addpath(\'%s\')'%(bp_op['matlab_tom_2008_code_path']) )
    session.run( 'addpath(\'%s\')'%(bp_op['matlab_tom_2005_code_path']) )



    import tomominer.model.util as MU
    v2 = MU.generate_toy_model()

    import tomominer.geometry.ang_loc as GA
    if ang is None:     ang = GA.random_rotation_angle_zyz()
    ang = N.array(ang)

    if loc is None:     loc = GA.random_translation(v2.shape, 0.1)
    loc = N.array(loc)

    import tomominer.geometry.rotate as GR
    v1 = GR.rotate(v2, angle=ang, loc_r=loc, default_val=0.0)


    import tomominer.simulation.back_projection_reconstruction as SB
    v1b = SB.do_reconstruction(session, v1, bp_op, verbose=False).astype(N.float)
    v2b = SB.do_reconstruction(session, v2, bp_op, verbose=False).astype(N.float)


    import tomominer.image.vol.wedge.util as IVWU
    m1 = IVWU.wedge_mask(v1b.shape, bp_op['model']['missing_wedge_angle']) * MU.sphere_mask(v1b.shape)      ;       m1 = m1.astype(N.float)
    m2 = N.copy(m1)


    return {'v1b':v1b, 'm1':m1, 'v2b':v2b, 'm2':m2, 'ang':ang, 'loc':loc}

