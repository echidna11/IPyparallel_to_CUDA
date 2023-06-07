#!/usr/bin/env python



'''
RASL collaborative alignment of 3D images (use inexact Augmented Lagrangian Mutiplier method)

see     Peng10 RASL Robust Alignment by Sparse and Low-rank Decomposition for Linearly Correlated Images
http://perception.csl.illinois.edu/matrix-rank/rasl.html

'''

import sys
import copy

import numpy as N
import numpy.linalg as NL

import tomominer.geometry.rotate as GR



def main_loop(vs, xi, lambdac, maxIter=1000, stoppingDelta=1e-5, inner_maxIter=1000, inner_tol=1e-10, rotate_default_val_mode=None, loc_rad_max=N.inf):

    lambda_t = lambdac / N.sqrt(vs[0].size)

    jdxi = N.array([0.02, 0.02, 0.02, 1, 1, 1])

    xi = copy.deepcopy(xi)
    prev_obj = N.inf
    for iter_i in range(maxIter):
        ost = outer_step(vs, xi, jdxi, lambda_t, inner_maxIter=inner_maxIter, inner_tol=inner_tol, rotate_default_val_mode=rotate_default_val_mode)

        xi = update_xi(xi, ost['delta_xi'], loc_rad_max=loc_rad_max)

        if N.abs(prev_obj - ost['obj']) < stoppingDelta:    break
        prev_obj = ost['obj']

        #jdxi = dxi_restrict(jdxi=jdxi, dxi=ost['delta_xi'])

        print ost['obj']

    return xi


# update transformation, restrict the translation withing certain radius
def update_xi(xi, dxi, loc_rad_max):
    for i in xi:
        xi_org = xi[i]
        assert NL.norm(xi_org[3:]) < loc_rad_max

        xi_t = xi_org + dxi[i]
        if NL.norm(xi_t[3:]) > loc_rad_max:  
            print 'warning: translation part of xi + dxi exceeds range', xi_t[3:]
            xi_t[3:] = xi_org[3:]

        xi[i] = xi_t

    return xi


# reduce the magnitude of jdxi, according to given dxi
def dxi_restrict(jdxi, dxi, ratio=2):
    jdxi = N.copy(jdxi)

    for img_i in dxi:
        for par_i in range(len(jdxi)):
            dxi_t = N.abs(dxi[img_i][par_i] / ratio)
            if N.abs(jdxi[par_i]) > dxi_t:     jdxi[par_i] = dxi_t

    return jdxi



'''
see ~/ln/electron/util/alignment/RASL/RASL_Code/rasl_main.m

parameters: v: images      xi: initial transformation parameters    djxi: incremental for calculating jacobian       maxIter: max iteration number      stoppingDelta: stopping tolerence
'''
def outer_step(vs, xi, jdxi, lambda_t, inner_maxIter=None, inner_tol=None, rotate_default_val_mode=None):
    D = N.zeros(   (vs[0].size, len(vs))     )
    J = {}
    for i in range(len(vs)):
        jr = jacobian_simple(vs[i], xi[i], dxi=jdxi, normalize=True, rotate_default_val_mode=rotate_default_val_mode)       # IMPORTANT: normalize=True is important to make the process converge!!
        D[:,i] = jr['vr'].flatten()
        J[i] = jr['J']

    Q = {};         R = {}
    for i in J:        Q[i], R[i] = NL.qr(J[i])

    A, E, delta_xi, numIterInnerEach, Y = inner_ialm(D, Q, lambda_t, tol=inner_tol, maxIter=inner_maxIter)

    for i in delta_xi:        delta_xi[i] = NL.pinv(R[i]).dot(delta_xi[i])          # the original matlab version uses inv(), we use psudo inverse pinv() to prevent inverse of sigular matrix

    assert      N.all(N.isfinite(A))
    A_s = NL.svd(A, compute_uv=False)
    obj = NL.norm(A_s, ord=1) + lambda_t * NL.norm(E.flatten(), ord=1)      # objective function value

    return {'delta_xi':delta_xi, 'obj':obj}





'''
# calculate the jacobian of a 3D image, given image rigid transform xi, and delta xi
# transform parameters xi are encoded as three rotational angles followed by three translations
'''
def jacobian_simple(v, xi, dxi, normalize=False, rotate_default_val_mode=None):
    
    if rotate_default_val_mode == 0:
        default_val = 0.0
    elif rotate_default_val_mode == 1:
        default_val = v.mean()
    else:
        raise Exception('rotate_default_val_mode')

    vr = GR.rotate(v, angle=xi[:3], loc_r=xi[3:], default_val=default_val)

    if normalize:   vr /= NL.norm(vr.flatten(), ord=2)

    J = N.zeros( (v.size, len(xi)) )
    for par_i in range(len(xi)):
        xit = N.copy(xi)
        xit[par_i] += dxi[par_i]
        vrt = GR.rotate(v, angle=xit[:3], loc_r=xit[3:], default_val=default_val)
        J[:,par_i] = vrt.flatten() - vr.flatten()

    return {'J':J, 'vr':vr}
        






'''
# see http://en.wikipedia.org/wiki/Chain_rule#Higher_dimensions
# todo: see if you can use the chain rule for Jacobian martix, J_{\mathbf{a}}(f \circ g) = J_{g(\mathbf{a})}(f) J_{\mathbf{a}}(g).
'''




'''
use inexact Augmented Lagrangian Mutiplier method to solve
min ||A||_* + \lambda |E|_1 + <Y_k, D + J*deltaTau -A-E> + \mu/2 ||D + J*deltaTau - A - E||_F^2
see ~/ln/electron/util/alignment/RASL/RASL_Code/rasl_inner_ialm.m
'''
def inner_ialm(D, J, lambda_t, tol=1e-7, maxIter=1000):

    m, n = D.shape

    Y = N.copy(D)
    norm_two = NL.norm(Y, ord=2)
    norm_inf = NL.norm(Y.flatten(), ord=N.inf) / lambda_t
    dual_norm = N.max([norm_two, norm_inf])
    Y /= dual_norm

    obj_v = N.dot(D.flatten(), Y.flatten())

    A_dual = N.zeros(   (m, n)  )
    E_dual = N.zeros(   (m, n)  )
    dt_dual = {}
    dt_dual_matrix = N.zeros(   (m, n)  )

    mu = 1.25 / NL.norm(D, ord=2)       ;       assert  N.isfinite(mu)
    rho = 1.25

    d_norm = NL.norm(D, ord='fro')

    for iter_i in range(maxIter):
        temp_T = D + dt_dual_matrix - E_dual + (1/mu)*Y
        U, diagS, V = NL.svd(temp_T, full_matrices=False)
        A_dual = U.dot(N.diag(pos(diagS-1/mu))).dot(V.T)        # see eq 6 of Lin10

        temp_T = D + dt_dual_matrix - A_dual + (1/mu)*Y
        E_dual = N.sign(temp_T) * pos( N.abs(temp_T) - lambda_t/mu )

        temp_T = D - E_dual - A_dual + (1/mu)*Y
        for i in range(n):
            dt_dual[i] =  - N.dot(J[i].T, temp_T[:,i])      # J*deltaTau should be close to -tempT
            dt_dual_matrix[:, i] = N.dot(J[i], dt_dual[i])

        Z = D + dt_dual_matrix - A_dual - E_dual
        Y = Y + mu * Z

        obj_v = N.dot(D.flatten(), Y.flatten())

        mu = mu * rho
        stoppingCriterion = NL.norm(Z, 'fro') / d_norm

        #print '\r', obj_v, '\t', stoppingCriterion        ;        sys.stdout.flush()
        

        if stoppingCriterion <= tol:
            break


    return (A_dual, E_dual, dt_dual, iter_i, Y)



def pos(A):
    return A * (A > 0)







'''

# test commands

import tomominer.model.util as MU

v = MU.generate_toy_model(dim_siz=64)

import tomominer.geometry.rotate as GR
import tomominer.geometry.ang_loc as GA

image_num = 10
translation_proportion = 0.2

import numpy as N

vs = {}
if False:
    for i in range(image_num):      vs[i] = GR.rotate(v, rm=GA.random_rotation_matrix(), loc_r=GA.random_translation(v.shape, translation_proportion), default_val=0.0)
else:
    for i in range(image_num):      vs[i] = GR.rotate(v, angle=N.random.random(3)*0.2, loc_r=GA.random_translation(v.shape, translation_proportion), default_val=0.0)



xi = {}
for i in range(image_num):
    xi[i] = N.zeros(6)
    #xi[i][:3] = GA.random_rotation_angle_zyz()
    #xi[i][3:] = GA.random_translation(v.shape, translation_proportion)



import tomominer.align.collaborative.nuclear.rasl as ACNR

import time
cur_time = time.time()
xi = ACNR.main_loop(vs=vs, xi=xi, lambdac=1, stoppingDelta=1e-4, inner_tol=1e-5, rotate_default_val_mode=0, loc_rad_max=(translation_proportion*N.max(v.shape)*10))
print time.time() - cur_time, 'sec'




import matplotlib
matplotlib.use('Qt4Agg')
import tomominer.image.vol.util as CV


for i in range(image_num):    CV.dsp_cub(GR.rotate(vs[i], angle=xi[i][:3], loc_r=xi[i][3:], default_val=0.0))


'''





