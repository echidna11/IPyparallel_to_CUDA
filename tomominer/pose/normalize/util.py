# functions for pose normalization

import numpy as N


import tomominer.geometry.rotate as GR


def center_mass(v):
    assert  N.all(v >= 0)
    
    m = v.sum()
    assert  m > 0

    v = v/m

    s = v.shape
    g = N.mgrid[0:s[0], 0:s[1], 0:s[2]]

    c = [None] * v.ndim
    for dim_i in range(v.ndim):
        c[dim_i] = (g[dim_i] * v).sum()


    return N.array(c)



'''
perform pose normalization through PCA, assume that the values in V is non-negative
'''
def pca(v, c, do_flip=False):
    assert      N.all(v >= 0)

    re = {'c':c}

    s = v.shape

    g = N.mgrid[0:s[0], 0:s[1], 0:s[2]]
    g = N.array(g, dtype=N.float)

    for i in range(len(g)):     g[i] -= c[i]

    gv = []
    
    gv = [g[0].flatten(), g[1].flatten(), g[2].flatten()]
    
    vv = v.flatten()
    gvw = [_*vv for _ in gv]

    wsm = N.dot(  N.array(gv),    N.array(gvw).T   )      # weighted sum matrix
    re['wsm'] = wsm


    # perform eigen decomposition, and order eigen values and eigen vectors according to decreasing order of magnitude of eignenvalues
    (eig_w, eig_v) = N.linalg.eig(wsm)

    i = N.argsort(-N.abs(eig_w))        # the resulting projection vectors in the rotation matrix will be ordered by the magnitude of eigenvalues
    eig_w = eig_w[i]
    eig_v = eig_v[:, i]


    re['w'] = eig_w

    if do_flip:
        re['v'] = flip_sign(v=v, c=c, r=eig_v)
    else:
        re['v'] = eig_v         # this is the rotation matrix. rotate vol by rm gives pose normalized vol

    return re




'''

# tests

import matplotlib
matplotlib.use('Qt4Agg')




import tomominer.model.util as MU
vo = MU.generate_toy_model()

import tomominer.geometry.ang_loc as GA
import tomominer.geometry.rotate as GR
vr = GR.rotate(vo, rm=GA.random_rotation_matrix(), loc_r=GA.random_translation(vo.shape, 0.2), default_val=0)
vr[vr < 0] = 0

import tomominer.pose.normalize.util as PNU
c = PNU.center_mass(vr)
p = PNU.pca(vr, c)

vp = GR.rotate(vr, rm=p['v'], c1=c, default_val=0)          # IMPORTANT: in rotate(), we have x*rm, therefore columnes of rm are vectors that projects row vector x

import tomominer.image.vol.util as CV

CV.dsp_cub(vr)
CV.dsp_cub(vp)


'''


'''
flip the sign of rotation matrix so that the positive sign of two largest eigenvectors point to the v of larger weight.
parameters: v: density map,  c: center, r: rotation matrix(columes are projection vectors)
'''
def flip_sign(v, c, r):
    
    mid = N.ceil(N.array(v.shape) / 2)

    best = None
    for s0 in [-1,1]:
        for s1 in [-1, 1]:
            p0 = N.array([s0, 0, 0])
            p1 = N.array([0, s1, 0])
            p2 = N.cross(p0, p1)

            rt = N.dot(r, N.array([p0, p1, p2]))

            vr = GR.rotate(v, rm=rt, c1=c, default_val=0)
            vr_sum = vr[mid[0]:, mid[1]:, mid[2]:].sum()

            if (best is None) or (vr_sum > best['vr_sum']):                best = {'s0':s0, 's1':s1, 'rt':N.copy(rt), 'vr_sum':vr_sum}

    assert      best is not None
    return best['rt']




