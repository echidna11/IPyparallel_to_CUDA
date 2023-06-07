'''
functions for rotating volumes
~/ln/tomominer/tomominer/geometry/rotate.py
'''

import numpy as N
import scipy.ndimage.interpolation as SNI

import tomominer.geometry.ang_loc as AA
import tomominer.image.vol.util as IVU


'''

use scipy.ndimage.interpolation.affine_transform() to rotate a volume.
 
rotation at center mid_co1, then move to center mid_co2
see http://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.interpolation.affine_transform.html#scipy.ndimage.interpolation.affine_transform


This function is written according to VolumeUtil.rotate_vol(). The major difference between scipy.ndimage.interpolation.affine_transform() and tformarray() is that affine_transform()'s rotation matrix is at LEFT side, not right side

 
Following should be corrected accordingly : use rotation matrix rm.  Let g2 and g1 be functions representing vol_t and vol, then we % have, g2(x) = g1( rm_t^(-1) * (x-mid_co2) + mid_co1), ie x1 = rm_t^(-1) (x2 -mid_co2) + mid_co1, where x1 is for g1, x2 is for g2. ie, x2 = rm_t * x1 - rm_t*mid_co1 + mid_co2. It is important to know that here we use right product, i.e. after transpose of ALL vectors into row vectors, x2 = x1 * rm  - mid_co1*rm + mid_co2, in this case rm' = rm^{-1} = rm_t  ,  note, the values outside vol are treated as NaN, you need to decide how to deal with them after each time you run rotate_vol !!!!!    

WARNING:    If v contains NaN values, then vr all become NaN

'''

def rotate(v, angle=None, rm=None, c1=None, c2=None, loc_r=None, siz2=None, default_val=float('NaN')):
    

    if angle is not None:
        assert      rm is None
        angle = N.array(angle, dtype=N.float).flatten()
        rm = AA.rotation_matrix_zyz(angle)

    if rm is None:      rm = N.eye(v.ndim)

    siz1 = N.array( v.shape, dtype=N.float )
    if c1 is None:
        c1 = (siz1-1) / 2.0                  # IMPORTANT: following python convension, in index starts from 0 to size-1!!! So (siz-1)/2 is real symmetry center of the volume
    else:
        c1 = c1.flatten()
    assert  c1.shape == (3,)

    if siz2 is None:    siz2 = siz1
    siz2 = N.array(siz2, dtype=N.float)

    if c2 is None:      
        c2 = (siz2-1) / 2.0               # IMPORTANT: following python convension, in index starts from 0 to size-1!!! So (siz-1)/2 is real symmetry center of the volume
    else:
        c2 = c2.flatten()
    assert c2.shape == (3,)

    if loc_r is not None:
        loc_r = N.array(loc_r, dtype=N.float).flatten()
        assert  loc_r.shape == (3,)
        c2 += loc_r


    
    c = -rm.dot(c2) + c1

    #rm_ext = N.hstack( (rm, c) )

    vr = SNI.affine_transform(input=v, matrix=rm, offset=c, output_shape=siz2.astype(N.int), cval=default_val)          # note: output_shape need to be integers, otherwise N.zeros will raise a warning
        
    return vr




'''

# test code


# -----------------------------------
# initializations


import sys
sys.path.append('/home/rcf-47/mxu/ln/tomominer/tomominer')
sys.path.append('/home/rcf-47/mxu/ln/tomominer/build/lib.linux-x86_64-2.7/tomominer')



import pymatlab
session = pymatlab.session_factory(options='-nodisplay')

session.run( 'addpath(\'%s\')'%('/home/rcf-47/mxu/ln/frequent_structure/code') )





# -----------------------------------
# random tests

import numpy as N
ang = 2 * N.pi * N.random.random(3)
loc_r = (N.random.random(3) - 0.5) * 10


import tomominer.model.util as MU
v = MU.generate_toy_model(dim_siz=32)


c1 = ( (N.array(v.shape)-1) / 2.0 ) + (N.random.random(3)-0.5) * 10
c2 = ( (N.array(v.shape)-1) / 2.0 ) + (N.random.random(3)-0.5) * 10



import tomominer.geometry.rotate as GR
vr = GR.rotate(v=v, angle=ang, c1=c1, c2=c2+loc_r, default_val=0.0)






session.run('clear all;')

session.putvalue('v', v)

session.putvalue('ang', ang)
session.putvalue('loc_r', loc_r)

session.putvalue('c1', c1)
session.putvalue('c2', c2)


session.run('vr = VolumeUtil.rotate_vol_pad0(v,  AngLocUtil.rotation_matrix_zyz(ang), c1, c2+loc_r);')


vr1 = session.getvalue('vr')

print N.corrcoef(vr.flatten(), vr1.flatten())





import tomominer.image.vol.util as CV

import matplotlib
matplotlib.use('Qt4Agg')




CV.dsp_cub(vr)
CV.dsp_cub(vr1)




#----------------------------------
# test only rotation

vr = GR.rotate(v=v, angle=ang, default_val=0.0)


session.run('clear all;')

session.putvalue('v', v)

session.putvalue('ang', ang)

session.run('vr = VolumeUtil.rotate_vol_pad0(v,  AngLocUtil.rotation_matrix_zyz(ang));')

vr1 = session.getvalue('vr')

print N.corrcoef(vr.flatten(), vr1.flatten())



'''


def rotate_pad_mean(v, angle=None, rm=None, c1=None, c2=None, loc_r=None, siz2=None):
    vr = rotate(v, angle=angle, rm=rm, c1=c1, c2=c2, loc_r=loc_r, siz2=siz2, default_val=float('NaN'))

    if False:
        vr[N.logical_not(N.isfinite(vr))] = vr[N.isfinite(vr)].mean()           # this may lead to problem when loc_r is too large and pushes the whole v outside the volume region
    else:
        vr[N.logical_not(N.isfinite(vr))] = v.mean()

    return vr

    

def rotate_pad_zero(v, angle=None, rm=None, c1=None, c2=None, loc_r=None, siz2=None):
    vr = rotate(v, angle=angle, rm=rm, c1=c1, c2=c2, loc_r=loc_r, siz2=siz2, default_val=float('NaN'))

    vr[N.logical_not(N.isfinite(vr))] = 0.0

    return vr


def rotate_mask(v, angle=None, rm=None):
    c1 = IVU.fft_mid_co(v.shape)
    c2 = N.copy(c1)

    vr = rotate(v, angle=angle, rm=rm, c1=c1, c2=c2, default_val=float('NaN'))

    vr[N.logical_not(N.isfinite(vr))] = 0.0

    return vr


