
import numpy as N

import tomominer.image.vol.util as gv
import tomominer.geometry.rotate as GR

def gauss_function(size, sigma):

    grid = gv.grid_displacement_to_center(size)
    dist_sq = gv.grid_distance_sq_to_center(grid)

    del grid


    g = (1 / ( (2 * N.pi)**(3.0/2.0)  * (sigma**3)) ) * N.exp( - (dist_sq)  / (2.0 * (sigma**2)))               # gauss function

    return g

def difference_of_gauss_function(size, sigma1, sigma2):

    grid = gv.grid_displacement_to_center(size)
    dist_sq = gv.grid_distance_sq_to_center(grid)

    del grid

    dog = (1 / ( (2 * N.pi)**(3.0/2.0)  * (sigma1**3)) ) * N.exp( - (dist_sq)  / (2.0 * (sigma1**2)))               # gauss function
    dog -= (1 / ( (2 * N.pi)**(3.0/2.0)  * (sigma2**3)) ) * N.exp( - (dist_sq)  / (2.0 * (sigma2**2)))

    return dog



# translated from SphericalHarmonicsUtil.generate_toy_vol()
def generate_toy_model(dim_siz=64, model_id=0):
    
    siz = N.array([dim_siz, dim_siz, dim_siz])

    mid = siz / 2.0

    xg = N.mgrid[0:siz[0], 0:siz[1], 0:siz[2]]

    if model_id == 0:
        # four gauss functions

        short_dia = 0.4
        mid_dia = 0.8
        long_dia = 1.2

        e0 = generate_toy_model__gaussian(dim_siz=dim_siz, xg=xg, xm=(mid + N.array([siz[0]/4.0, 0.0, 0.0])), dias=[long_dia, short_dia, short_dia])
        e1 = generate_toy_model__gaussian(dim_siz=dim_siz, xg=xg, xm=(mid + N.array([0.0, siz[1]/4.0, 0.0])), dias=[short_dia, long_dia, short_dia])
        e2 = generate_toy_model__gaussian(dim_siz=dim_siz, xg=xg, xm=(mid + N.array([0.0, 0.0, siz[2]/4.0])), dias=[short_dia, short_dia, long_dia])

        e3 = GR.rotate_pad_zero(N.array(e0, order='F'), angle=N.array([N.pi/4.0, 0.0, 0.0]), loc_r=N.array([0.0, 0.0, 0.0]))

        e = e0 + e1 + e2 + e3

    elif model_id == 101:
        # three axis with markers, useful to inspect rotation

        e = N.zeros(siz)
        e[int(siz[0]/2):, int(siz[1]/2), int(siz[2]/2)] = 50;
        e[int(siz[0]/2), int(siz[1]/2):, int(siz[2]/2)] = 100;          e[int(siz[0]*0.45):int(siz[0]*0.55), int(siz[1]*0.8), int(siz[2]*0.45):int(siz[2]*0.55)] = 100;
        e[int(siz[0]/2), int(siz[1]/2), int(siz[2]/2):] = 150;          e[int(siz[0]*0.45):int(siz[0]*0.55), int(siz[1]*0.45):int(siz[1]*0.55), int(siz[2]*0.8)] = 150;        e[int(siz[0]*0.45):int(siz[0]*0.55), int(siz[1]*0.45):int(siz[1]*0.55), int(siz[2]*0.6)] = 150;

    else:
        raise

    cm = gv.center_mass(e)
    e = GR.rotate_pad_zero(N.array(e, order='F'), angle=N.array([0.0, 0.0, 0.0]), loc_r=(mid - cm))


    return e



def generate_toy_model__gaussian(dim_siz, xg, xm, dias):
    x = N.zeros(xg.shape)
    for dim_i in range(3):      x[dim_i] =  xg[dim_i] - xm[dim_i]
    xs = N.array([ x[0] / (dim_siz*dias[0]), x[1] / (dim_siz*dias[1]), x[2] / (dim_siz*dias[2]) ])
    e = N.exp(- N.sum( xs * x, axis=0) )

    return e



"""
#test commands

import sys
sys.path.append('/home/rcf-47/mxu/ln/tomo/py')

import tomominer.model.util as mu

import tomominer.image.vol.util as gv

import matplotlib
matplotlib.use('Qt4Agg')

gv.dspcub(mu.generate_toy_model(),  block_draw=True)


"""

    



def sphere_mask(shape, center=None, radius=None, smooth_sigma=None):

    shape = N.array(shape)

    v = N.zeros(shape)

    if center is None:          center = (shape-1) / 2.0              # IMPORTANT: following python convension, in index starts from 0 to size-1 !!! So (siz-1)/2 is real symmetry center of the volume

    center = N.array(center)

    if radius is None:          radius = N.min(shape/2.0)

    
    grid = gv.grid_displacement_to_center(shape, mid_co=center)
    dist = gv.grid_distance_to_center(grid)

    v[dist <= radius] = 1.0

    if smooth_sigma is not None:
        '''
        smooth a bit, see
        /home/rcf-47/mxu/ln/electron/matlab_tom/2008/TOM_Release_2008/Geom/tom_spheremask.m
        '''

        assert smooth_sigma > 0
        v_s = N.exp(-((dist -radius)/smooth_sigma) ** 2)
        v_s[v_s < N.exp(-3)] = 0.0      # use a cutoff of -3 looks nicer, although the tom toolbox uses -2
        v[dist > radius] = v_s[dist > radius]


    return v


"""
#test commands

import sys
sys.path.append('/home/rcf-47/mxu/ln/tomominer/tomominer')
sys.path.append('/home/rcf-47/mxu/ln/tomominer/build/lib.linux-x86_64-2.7/tomominer')

import tomominer.model.util as mu
import tomominer.image.vol.util as gv

import matplotlib
matplotlib.use('Qt4Agg')

gv.dsp_cub(mu.sphere_mask( (30,40,50) ))
gv.dsp_cub(mu.sphere_mask( (30,40,50), radius=5, smooth_sigma=2 ))



"""





def rectangle_mask(shape, p0, p1):
    shape = N.array(shape)

    v = N.zeros(shape)

    v[p0[0]:p1[0], p0[1]:p1[1], p0[2]:p1[2]] = 1.0

    return v


# generate a 3D volume where only the boundary regions are 1, other places are zeros
def boundary_mask(shape):
    shape = N.array(shape)
    m = N.zeros(shape)

    m[0,:,:] = 1
    m[-1,:,:] = 1
    m[:,0,:] = 1
    m[:,-1,:] = 1
    m[:,:,0] = 1
    m[:,:,-1] = 1

    return m





def cylinder(shape, radius=None, half_height=None, center=None, direction=2):
    v = N.zeros(shape)
    shape_f = N.array(shape, dtype=N.float)

    if center is None:
        center = shape_f / 2

    if half_height == None:
        half_height = min(shape_f[direction]-center[direction], center[direction])

    if radius is None:
        radius = N.inf
        for d in range(3):
            if d == direction:  continue
            radius = min(radius, center[d], shape_f[d]-center[d])
    radius_sq = radius * radius

    g = N.mgrid[0:shape[0], 0:shape[1], 0:shape[2]]

    dis_sum_sq = N.zeros(shape)
    for d in range(3):
        if d == direction:  continue
        dis_sum_sq += N.square(g[d] - center[d])

    v[dis_sum_sq <= radius_sq] = 1.0

    v[N.abs(g[direction] - center[direction]) > half_height] = 0.0

    return v


'''

import matplotlib
matplotlib.use('Qt4Agg')

import tomominer.model.util as MU
t = MU.cylinder(shape=(20,40,60), direction=2)

import tomominer.image.vol.util as IVU
IVU.dsp_cub(t)


'''







    
