
'''

utility functions that supports point hull processing

~/ln/tomominer/tomominer/geometry/point_cloud/util.py

'''


# given a binary 3D mask, find the convex hull of corresponding points, then return a mask that represents the convex hull
def mask_convex_hull(m):
    x = N.where(m > 0.5)
    x = N.array(x).T

    g = N.mgrid[0:m.shape[0], 0:m.shape[1], 0:m.shape[2]]
    g = N.array([g[0].flatten(), g[1].flatten(), g[2].flatten()]).T

    from scipy import spatial
    hull = spatial.Delaunay(x)
    hull_flag = in_hull(g, hull)
    return N.reshape(hull_flag, m.shape)



'''
# test code

import numpy as N
v = N.zeros([10,11,12])
v[2,2,2] = 1
v[6,2,2] = 1
v[2,8,2] = 1
v[2,2,10] = 1


from tomominer.geometry.point_cloud.util import mask_convex_hull

vc = mask_convex_hull(v)

from tomominer.image.vol.util import dsp_cub
dsp_cub(v)
dsp_cub(vc)


'''



def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimension
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
    coordinates of `M` points in `K`dimension for which a Delaunay triangulation
    will be computed
    """
    from scipy.spatial import Delaunay
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)>=0






'''
rigid transform of points in x (each point is a row), with rotation matrix rm and middle coordinates mid_co1 and mid_co2
the defination is according to  VolumeUtil.rotate_vol()
Warning: instead of using matlab's own isosurface, must use Surface.vol_isosurface() to convert a density map to surface representation!!
Warning: because python passes references instead of copies, need to make a copy first for any array value that is going to change!!
'''
import tomominer.geometry.ang_loc as GA
import numpy as N
def rotate_translate(x0, angle=None, rm=None, c0=None, c1=None, loc_r=None):

    if angle is not None:
        assert      rm is None
        rm = GA.rotation_matrix_zyz(angle)

    if rm is None:      rm = N.eye(x0.shape[1])

    assert c0 is not None

    if c1 is None:      c1 = N.copy(c0)

    if loc_r is not None:       
        c1t = c1 + loc_r
    else:
        c1t = c1


    ct = c0.dot(rm) - c1t
    ct = ct.reshape((1,-1))
    ct = N.tile(ct, (x0.shape[0], 1))

    x1 = x0.dot(rm) - ct

    return x1


'''

# testing code, test the consistancy to the matlab version


import numpy as N
x0 = N.random.random( (10, 3) )
ang = N.random.random(3) * 2 * N.pi
loc = N.random.random(3) * 10
c0 = N.random.random(3) * 10
c1 = N.random.random(3) * 10


import tomominer.geometry.point_cloud.util as GP
xr = GP.rotate_translate(x0, angle=ang, c0=c0, c1=c1, loc_r=loc)
xr1 = GP.rotate_translate(x0, angle=ang, c0=c0, c1=c1+loc)
print xr-xr1


import pymatlab
session = pymatlab.session_factory(options='-nodesktop -nodisplay')
session.run( 'addpath(\'%s\')'%('/home/rcf-47/mxu/ln/frequent_structure/code',))
session.putvalue('ang', ang)
session.putvalue('loc', loc)
session.putvalue('c0', c0)
session.putvalue('c1', c1)
session.putvalue('x0', x0)

session.run('c0 = reshape(c0, 1,[])')
session.run('c1 = reshape(c1, 1,[])')
session.run('loc = reshape(loc, 1,[])')

session.run('rm = AngLocUtil.rotation_matrix_zyz(ang)')
session.run('xr = PointCloud.rotate_translate(x0, rm, c0, c1+loc)')

xrt = session.getvalue('xr')

print   xr - xrt



'''




'''
given a set of points sampled on a surface (each point is one row), and these points can be weighted. Calculate weighted PCA to get surface norm
'''
import copy
def surface_norm(x, w=None):
    assert  x.shape[1] == 3
    x = copy.deepcopy(x)

    if w is not None:
        for i in range(x.shape[1]):     x[:,i] *= w

    eig_w, eig_v = N.linalg.eig(N.dot(x.T, x))
    i = N.argmin(N.abs(eig_w))

    return eig_v[:,i]



'''
# test code

import numpy as N
x = [[0,0,0], [1,0,0], [0,1,0]]
x = N.array(x)

import tomominer.geometry.point_cloud.util as GPU
print GPU.surface_norm(x)

'''






