# according to AngLocUtil.m

import math
import numpy as N


# the total difference combining both angle and location
def angle_zyz_translation_difference(ang1=N.zeros(3), loc1_r=N.zeros(3), ang2=N.zeros(3), loc2_r=N.zeros(3)):
    if loc1_r is None:            loc1_r = N.zeros(ang1.shape)
    if loc2_r is None:            loc2_r = N.zeros(ang2.shape)

    if ang1.ndim == 2:

        n = len(ang1);
        assert(n == len(loc1_r));  assert(n == len(ang2));  assert(n == len(loc2_r));
        
        dif_d = N.zeros(n);
        for i in range(n): 
            rm1 = rotation_matrix_zyz(ang1[i,:])
            rm2 = rotation_matrix_zyz(ang2[i,:])
            loc1_r_t = N.tile(loc1_r[i, :], [3, 1])
            loc2_r_t = N.tile(loc2_r[i, :], [3, 1])
            
            dif_m = (rm1 * (N.eye(3) - loc1_r_t)).transpose() - (rm2 * (N.eye(3) - loc2_r_t)).transpose()     
            dif_d[i] = math.sqrt(N.square(dif_m).sum())

    elif ang1.ndim == 1:
        # when ang1 is only one dimension, it should be a single vector of lenth 3
        assert(len(ang1) == 3);     assert(len(loc1_r) == 3);       assert(len(ang2) == 3);         assert(len(loc2_r) == 3);
        rm1 = rotation_matrix_zyz(ang1)
        rm2 = rotation_matrix_zyz(ang2)
        loc1_r_t = N.array([loc1_r, loc1_r, loc1_r])
        loc2_r_t = N.array([loc2_r, loc2_r, loc2_r])

        dif_m = (rm1.dot(N.eye(3) - loc1_r_t)).transpose() - (rm2.dot(N.eye(3) - loc2_r_t)).transpose()
        dif_d = math.sqrt(N.square(dif_m).sum())
    
    else:
        raise    


    return dif_d    

def angle_difference_zyz(ang1, ang2):
    return angle_zyz_translation_difference(ang1=ang1, loc1_r=None, ang2=ang2, loc2_r=None)

 
def rotation_matrix_zyz(ang):
    phi = ang[0];       theta = ang[1];     psi_t = ang[2];
    
    a1 = rotation_matrix_axis(2, psi_t)       # first rotate about z axis for angle psi_t
    a2 = rotation_matrix_axis(1, theta)
    a3 = rotation_matrix_axis(2, phi)
    
    rm = a3.dot(a2).dot(a1)      # for matrix left multiplication
    
    rm = rm.transpose()       # note: transform because tformarray use right matrix multiplication

    return rm

'''

# test code

import numpy as N
ang = N.random.random(3) * 2 * N.pi

import tomominer.geometry.ang_loc as GA
rm = GA.rotation_matrix_zyz(ang)

import pymatlab
session = pymatlab.session_factory(options='-nodesktop -nodisplay')
session.run( 'addpath(\'%s\')'%('/home/rcf-47/mxu/ln/frequent_structure/code',))
session.putvalue('ang', ang)

session.run('rm = AngLocUtil.rotation_matrix_zyz(ang)')
rmt = session.getvalue('rm')

print   rm - rmt


'''



def rotation_matrix_axis(dim, theta):
    # following are left handed system (clockwise rotation)
    # IMPORTANT: different to MATLAB version, this dim starts from 0, instead of 1
    if dim == 0:        # x-axis
        rm = N.array(  [[1.0, 0.0, 0.0], [0.0, math.cos(theta), -math.sin(theta)], [0.0, math.sin(theta), math.cos(theta)]]  )
    elif dim == 1:    # y-axis
        rm = N.array(  [[math.cos(theta), 0.0, math.sin(theta)], [0.0, 1.0, 0.0], [-math.sin(theta), 0.0, math.cos(theta)]]  )
    elif dim == 2:        # z-axis
        rm = N.array(  [[math.cos(theta), -math.sin(theta), 0.0], [math.sin(theta), math.cos(theta), 0.0], [0.0, 0.0, 1.0]]  )
    else:
        raise    
    
    return rm


def rotation_matrix_zyz_normalized_angle(rm):

    assert(all(N.isreal(rm.flatten())));     assert(rm.shape == (3,3));
    
    cos_theta = rm[2, 2]
    if N.abs(cos_theta) > 1.0:
        # warning(sprintf('cos_theta %g', cos_theta));
        cos_theta = N.sign(cos_theta);
    
    theta = N.arctan2(N.sqrt(1.0 - (cos_theta*cos_theta) ), cos_theta);
    
    if N.abs(cos_theta) < (1.0 - (1e-10)) :          # use a small epslon to increase numerical stability when abs(cos_theta) is very close to 1!!!!
        phi = N.arctan2(rm[2,1], rm[2,0]);
        psi_t = N.arctan2(rm[1,2], -rm[0,2]);
    else:
        theta = 0.0
        phi = 0.0
        psi_t = N.arctan2(rm[0,1], rm[1,1])
    
    ang = N.array([phi, theta, psi_t], dtype=N.float)

    return ang


# translated form AngLocUtil.reverse_transform()
def reverse_transform(rm, loc_r):
    
    rev_rm = rm.T
    rev_loc_r = ( - N.dot(loc_r, rev_rm) )

    return (rev_rm, rev_loc_r)


def reverse_transform_ang_loc(ang, loc_r):
    rm = rotation_matrix_zyz(ang)

    rev_rm, rev_loc_r = reverse_transform(rm, loc_r)

    return ( rotation_matrix_zyz_normalized_angle(rev_rm), rev_loc_r )

'''
# testing commands, compare the input and output to the matlab version

import sys
sys.path.append('/home/rcf-47/mxu/ln/tomo/py')

import numpy as np
ang = np.random.random(3) * np.pi * 2
loc = np.random.random(3) * 10

import tomominer.geometry.ang_loc as aal
rev_ang, rev_loc = aal.reverse_transform_ang_loc(ang, loc)

print ang, loc
print rev_ang, rev_loc

# compare the input and output to the matlab version, for example
# ang = [ 5.17819545  0.39606097  5.48179727];   loc = [ 6.16511819  6.72420034  3.07804606];    [rev_rm, rev_loc] = AngLocUtil.reverse_transform(AngLocUtil.rotation_matrix_zyz(ang), loc);   [AngLocUtil.rotation_matrix_zyz_normalized_angle(rev_rm); rev_loc]

'''


'''
suppose two rigid transforms, T1, T2, get combined transform T3 such
that T3 x = T2 T1 x
the rotation matrix is a right multiplication.
This means that suppose v3t is obtained from applying T3 to v1, v2 is obtained from applying T1 to v1, and v3 is obtained by applying T2 to v2. Then we have v3t == v3 
see 110221_fourier_pose_norm.tex subsection "Compositing two rigid transforms"

code translated from AngLocUtil.combined_transform()
'''
def combined_transform(rm1, loc_r1, rm2, loc_r2):
    assert  rm2.shape == (3,3)
    if loc_r2.shape[0] != 1:   loc_r2 = loc_r2.reshape((1,-1)) 
    assert  loc_r2.shape == (1,3)

    assert  rm1.shape == (3,3)
    if loc_r1.shape[0] != 1:   loc_r1 = loc_r1.reshape((1,-1))
    assert  loc_r1.shape == (1,3)

    
    rm2tr = rm2.T
    loc_r2tr = loc_r2.T;
    rm1tr = rm1.T
    loc_r1tr = loc_r1.T
    
    rm3tr = N.dot(rm2tr, rm1tr)
    loc_r3tr = N.dot(rm2tr, loc_r1tr) + loc_r2tr;
    
    rm3 = rm3tr.T
    loc_r3 = loc_r3tr.T
    return (rm3, loc_r3)


'''
# test code

import tomominer.model.util as MU
v1 = MU.generate_toy_model()

import tomominer.geometry.ang_loc as GA

rm1 = GA.random_rotation_matrix()
l1 = GA.random_translation(v1.shape, 0.2)

rm2 = GA.random_rotation_matrix()
l2 = GA.random_translation(v1.shape, 0.2)

import tomominer.geometry.rotate as GR
v2 = GR.rotate(v1, rm=rm1, loc_r=l1, default_val=0.0)
v3 = GR.rotate(v2, rm=rm2, loc_r=l2, default_val=0.0)

rm3t, l3t = GA.combined_transform(rm1=rm1, loc_r1=l1, rm2=rm2, loc_r2=l2)

v3t = GR.rotate(v1, rm=rm3t, loc_r=l3t, default_val=0.0)

import scipy.stats as SS
SS.pearsonr(v3.flatten(), v3t.flatten())


'''

# generate a random rotation matrix using SVD on a random matrix
def random_rotation_matrix():
    m = N.random.random( (3,3) )
    u,s,v = N.linalg.svd(m)

    return u

def random_rotation_angle_zyz():
    rm = random_rotation_matrix()
    return rotation_matrix_zyz_normalized_angle(rm)


# generate a random translation
def random_translation(size, proportion):
    size = N.array(size)
    return  (N.random.random(len(size))-0.5) * size * proportion






'''
use cosin law to calculate angle between two vectors
if take_abs=True, then the resurned angle will be no more than 90 degrees
'''

def angle_between_vectors(v0, v1, take_abs=False):
    v0 = v0 / N.linalg.norm(v0)
    v1 = v1 / N.linalg.norm(v1)

    d = N.dot(v0,v1)
    if take_abs:    d = N.abs(d)

    return N.arccos(d)


