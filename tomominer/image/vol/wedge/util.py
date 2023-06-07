
'''
~/ln/tomominer/tomominer/image/vol/wedge/util.py
'''


#--------------------------------------------------------------------------------
# functions for checking missing wedge


import warnings

import numpy as N
import tomominer.image.vol.util as IVU
import tomominer.model.util as MU


# this is used to find the minimum value in a (averaged) missing wedge mask region
def wedge_mask_min(mask, radius_ratio=0.5):

    grid_dist = IVU.grid_distance_to_center( grid_displacement_to_center(mask.shape) )
    radius = min(mask.shape) / 2
    mask_min = mask[ grid_dist< (radius_ratio * radius) ].min()
    return mask_min


# construct a missing wedge mask, see tom_wedge,
# angle represents the angle range of MISSING WEDGE region, the larger, the more missing wedge region!!!
# Currently the parameter 'direction' refer to the tilt axis
def wedge_mask(size, ang1, ang2=None, direction=1, sphere_mask=True, verbose=False):

    warnings.warn("The definition of wedge mask is still ambiguous")        # should define both tilt axis and electron beam (missing wedge) direction

    if ang2 is None:
        ang2 = float(N.abs(ang1))
        ang1 = -ang2

    else:
        assert      ang1 < 0
        assert      ang2 > 0


    if verbose:     print 'image.vol.wedge.util.wedge_mask()', 'ang1', ang1, 'ang2', ang2, 'direction', direction, 'sphere_mask', sphere_mask

    ang1 = (ang1 / 180.0) * N.pi
    ang2 = (ang2 / 180.0) * N.pi

    g = IVU.grid_displacement_to_center(size=size, mid_co=IVU.fft_mid_co(siz=size))

    if direction==0:
        # y-z plane
        x0 = g[1]           # y axis
        x1 = g[2]           # z axis

    elif direction==1:
        # x-z plane
        x0 = g[0]       # x axis
        x1 = g[2]       # z axis

    elif direction==2:
        # x-y plane
        x0 = g[0]       # x axis
        x1 = g[1]       # y axis



    m = N.zeros(size, dtype=float)

    m[ N.logical_and(x0 >= (N.tan(ang2)*x1), x0 >= (N.tan(ang1)*x1)) ] = 1.0
    m[ N.logical_and(x0 <= (N.tan(ang1)*x1), x0 <= (N.tan(ang2)*x1)) ] = 1.0

    if sphere_mask:    m *= MU.sphere_mask(m.shape)

    return m


'''
# code for validating wedge_mask(), by comparing to tom_wedge()




import os,sys
sys.path.append(os.path.join( os.getenv('HOME'), 'ln/tomominer/tomominer' ))


import numpy as N
siz = [64, 64, 64]
angle = 10.0

import tomominer.image.vol.wedge.util as IVWU
m = IVWU.wedge_mask(siz, angle)



import matplotlib
matplotlib.use('Qt4Agg')

import tomominer.image.vol.util as CV
CV.dsp_cub(m, view_dir=1)




import pymatlab
session = pymatlab.session_factory(options='-nodesktop -nodisplay')

tom_code_path = '/home/rcf-47/mxu/proj/imaging/electron/matlab_tom/2008/TOM_Release_2008'
session.run( 'addpath( genpath(\'%s\') )'%(tom_code_path) )

session.run('clear all;')

session.putvalue('siz', siz)
session.putvalue('ang', [angle])
session.run('m = tom_wedge(zeros(reshape(siz,1,[])), ang);')

mt = session.getvalue('m')

print (mt-m).max(), (mt-m).min()

CV.dsp_cub(mt-m)


import matplotlib.pyplot as plt
plt.show()














#--------------------------------------------------------------------

import os,sys
sys.path.append(os.path.join( os.getenv('HOME'), 'ln/tomominer/tomominer' ))


import numpy as N
siz = [64, 64, 64]

import tomominer.image.vol.wedge.util as IVWU
m = IVWU.wedge_mask(siz, ang1=-20, ang2=45, direction=2)



import matplotlib
matplotlib.use('Qt4Agg')

import tomominer.image.vol.util as CV
CV.dsp_cub(m)





'''



