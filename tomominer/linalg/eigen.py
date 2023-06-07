

# functions for calculating eigen values and eigen vectors

import numpy as N

# given a batch of 3 by 3 symmetric matrixes, calculate eigenvalues for each matrix.  See paper   Eberly06 Eigensystems for 3 by 3 Symmetric Matrices
# IMPORTANT: A is preconditioned by dividing by its maximum magnitude entry when that maximum is larger than 1
def eigen_value_3_symmetric_batch(A):

    for At in A:
        for Att in At:
            if Att is None: continue
            assert  Att.max() <= 1.0

    inv3 = 1.0 / 3.0
    root3 = float(  N.sqrt(3.0) )

    a00 = A[0][0].flatten()
    a01 = A[0][1].flatten()
    a02 = A[0][2].flatten()
    a11 = A[1][1].flatten()
    a12 = A[1][2].flatten()
    a22 = A[2][2].flatten()

    c0 = a00*a11*a22 + 2.0*a01*a02*a12 - a00*a12*a12 - a11*a02*a02 - a22*a01*a01
    c1 = a00*a11 - a01*a01 + a00*a22 - a02*a02 + a11*a22 - a12*a12
    c2 = a00 + a11 + a22

    c2Div3 = c2*inv3
    aDiv3 = (c1 - c2*c2Div3)*inv3

    aDiv3[aDiv3 > 0.0] = 0.0

    mbDiv2 = 0.5*(c0 + c2Div3*(2.0*c2Div3*c2Div3 - c1))

    q = mbDiv2*mbDiv2 + aDiv3*aDiv3*aDiv3

    q[q > 0.0] = 0.0

    magnitude = N.sqrt(-aDiv3)
    angle = N.arctan2(N.sqrt(-q), mbDiv2) *inv3

    cs = N.cos(angle)
    sn = N.sin(angle)

    root = N.zeros( (a00.size, 3) )
    root[:, 0] = c2Div3 + 2.0*magnitude*cs
    root[:, 1] = c2Div3 - magnitude*(cs + root3*sn)
    root[:, 2] = c2Div3 - magnitude*(cs - root3*sn)

    root_i = N.argsort(-N.abs(root), axis=1)     # Sort the roots here to obtain abs(root[0]) >= (root[1]) >= (root[2]).
    root = N.array(     [root[_,root_i[_]] for _ in range(root.shape[0])]      )

    re = [None] * 3
    for dim in range(3):
        re[dim] = N.reshape(root[:, dim], A[0][0].shape)

    return re


# this function is used to verify def eigen_value_3_symmetric_batch()
def eigen_value_3_symmetric_batch__simple(A):

    root = []
    for dim in range(3):        root.append(N.zeros(A[0][0].shape))

    for i0 in range(A[0][0].shape[0]):
        for i1 in range(A[0][0].shape[1]):
            for i2 in range(A[0][0].shape[2]):
                a00 = A[0][0][i0, i1, i2]
                a01 = A[0][1][i0, i1, i2]       ;       a10 = a01
                a02 = A[0][2][i0, i1, i2]       ;       a20 = a02
                a11 = A[1][1][i0, i1, i2]
                a12 = A[1][2][i0, i1, i2]       ;       a21 = a12
                a22 = A[2][2][i0, i1, i2]

                mat = N.array(      [[a00, a01, a02], [a10, a11, a12], [a20, a21, a22]] )
                e = N.linalg.eigvalsh(mat)      ;       assert      e.size == 3
                ei = N.argsort(-N.abs(e))
                es = e[ei]           # Sort the roots here to obtain abs(root[0]) >= (root[1]) >= (root[2])

                for dim in range(3):                root[dim][i0,i1,i2] = es[dim]

    return root




'''

# testing both eigen_value_3_symmetric_batch and filtering.differntial.hessian_3d

import numpy as N
v = N.random.random( (5, 5, 5) )

import tomominer.filter.differential as FD
h = FD.hessian_3d(v)

m = FD.hessian_3d__max_magnitude(h)
ht = FD.hessian_3d__normalize(h, m)


import tomominer.linalg.eigen as LE
e = LE.eigen_value_3_symmetric_batch(ht)

r = N.array( [e[_].flatten() for _ in range(3)] ).T
print N.abs(r)      # use this to check if the magnitues are in decreasing order


es = LE.eigen_value_3_symmetric_batch__simple(ht)
rs = N.array( [es[_].flatten() for _ in range(3)] ).T
print rs

print N.abs(r-rs)

'''












'''

given a eigen value, get corresponding eigen vector by solving linear systems of equation:


MATLAB symbolic math derivation:

clear all;
syms a00 a01 a02 a11 a12 a22
A = [a00 a01 a02; a01 a11 a12; a02 a12 a22];
linsolve(A, [0; 0; 0])
null(A)

[U S V] = svd(A);
V(:,end)



-----------------------------------------------------
hand derivation:

step:
a_{00}x_{0} + a_{01}x_{1} + a_{02}x_{2}=0\\
a_{01}x_{0} + a_{11}x_{1} + a_{12}x_{2}=0\\
a_{02}x_{0} + a_{12}x_{1} + a_{22}x_{2}=0\\

step:
a_{12}a_{00}x_{0} + a_{12}a_{01}x_{1} + a_{12}a_{02}x_{2}=0\\
a_{02}a_{01}x_{0} + a_{02}a_{11}x_{1} + a_{02}a_{12}x_{2}=0\\

-->
(a_{12}a_{00}-a_{02}a_{01})x_{0} + (a_{12}a_{01} - a_{02}a_{11})x_{1} = 0

-->
x_{1} = \frac{a_{02}a_{01}-a_{12}a_{00}}{a_{12}a_{01} - a_{02}a_{11}} x_{0}


step:
a_{12}a_{01}x_{0} + a_{12}a_{11}x_{1} + a_{12}a_{12}x_{2}=0\\
a_{11}a_{02}x_{0} + a_{11}a_{12}x_{1} + a_{11}a_{22}x_{2}=0\\

-->
(a_{12}a_{01}-a_{11}a_{02})x_{0} + (a_{12}a_{12}-a_{11}a_{22})x_{2}=0

-->
x_{2}=\frac{a_{11}a_{02}-a_{12}a_{01}}{a_{12}a_{12}-a_{11}a_{22}} x_{0}


so the eigen vector should be normalized from
[1, \frac{a_{02}a_{01}-a_{12}a_{00}}{a_{12}a_{01} - a_{02}a_{11}}, \frac{a_{11}a_{02}-a_{12}a_{01}}{a_{12}a_{12}-a_{11}a_{22}}]




-----------------------------------------------------
another hand derivation:

step:
a_{00} + a_{01}x_{1} + a_{02}x_{2}=0\\
a_{01} + a_{11}x_{1} + a_{12}x_{2}=0\\
a_{02} + a_{12}x_{1} + a_{22}x_{2}=0\\


step:
a_{12}a_{00} + a_{12}a_{01}x_{1} + a_{12}a_{02}x_{2}=0\\
a_{02}a_{01} + a_{02}a_{11}x_{1} + a_{02}a_{12}x_{2}=0\\

-->
a_{12}a_{00} + a_{12}a_{01}x_{1} = a_{02}a_{01} + a_{02}a_{11}x_{1}

-->
x_{1} = \frac{a_{02}a_{01}-a_{12}a_{00}}{a_{12}a_{01}-a_{02}a_{11}}
x_{2} = \frac{a_{11}a_{00}-a_{01}a_{01}}{a_{01}a_{12}-a_{11}a_{02}}


alternate step:
a_{22}a_{01} + a_{22}a_{11}x_{1} + a_{22}a_{12}x_{2}=0\\
a_{12}a_{02} + a_{12}a_{12}x_{1} + a_{12}a_{22}x_{2}=0\\

-->
x_{1} = \frac{a_{12}a_{02}-a_{22}a_{01}}{a_{22}a_{11}-a_{12}a_{12}}
x_{2} = \frac{a_{11}a_{02}-a_{12}a_{01}}{a_{12}a_{12}-a_{11}a_{22}}


'''

# parameters: A: symmetric real matrics,  lamb eigen values
# matrix a = A - lamb*I
def eigen_vector_3_symmetric_batch__given_eigen_value(A, lamb, epsilon=1e-5):

    
    a00 = A[0][0] - lamb
    a01 = A[0][1]
    a02 = A[0][2]
    a11 = A[1][1] - lamb
    a12 = A[1][2]
    a22 = A[2][2] - lamb


    d = a12*a01 - a02*a11

    c1 = a02*a01 - a12*a00
    c2 = a11*a00 - a01*a01

    ind = N.abs(d) > epsilon

    x0 = N.zeros(a00.shape) + float('nan')
    x0[ind] = 1.0

    x1 = N.zeros(a00.shape) + float('nan')
    x1[ind] = c1[ind] / d[ind]

    x2 = N.zeros(a00.shape) + float('nan')
    x2[ind] = c2[ind] / d[ind]


    # normalize vectors
    s = N.sqrt(x0*x0 + x1*x1 + x2*x2)
    x0 = x0 / s
    x1 = x1 / s
    x2 = x2 / s


    return [x0, x1, x2]


# this function is used to verify eigen_vector_3_symmetric_batch__given_eigen_value()
def eigen_vector_3_symmetric_batch__simple(A):

    root = []
    for dim in range(3):        root.append(N.zeros(A[0][0].shape))

    evs = []
    for dim in range(3):
        evs.append([])
        for i0 in range(A[0][0].shape[0]):
            evs[dim].append([])

            for i1 in range(A[0][0].shape[1]):
                evs[dim][i0].append([])

                for i2 in range(A[0][0].shape[2]):
                    evs[dim][i0][i1].append([])


    for i0 in range(A[0][0].shape[0]):
        for i1 in range(A[0][0].shape[1]):
            for i2 in range(A[0][0].shape[2]):
                a00 = A[0][0][i0, i1, i2]
                a01 = A[0][1][i0, i1, i2]       ;       a10 = a01
                a02 = A[0][2][i0, i1, i2]       ;       a20 = a02
                a11 = A[1][1][i0, i1, i2]
                a12 = A[1][2][i0, i1, i2]       ;       a21 = a12
                a22 = A[2][2][i0, i1, i2]

                mat = N.array(      [[a00, a01, a02], [a10, a11, a12], [a20, a21, a22]] )
                ew, ev = N.linalg.eig(mat)      ;       assert      ew.size == 3
                ei = N.argsort(-N.abs(ew))
                ew = ew[ei]           # Sort the roots here to obtain abs(root[0]) >= (root[1]) >= (root[2])
                ev = ev[:,ei]

                for dim in range(3):                root[dim][i0,i1,i2] = ew[dim]
                for dim in range(3):                evs[dim][i0][i1][i2] = ev[:,dim].flatten()

    return {'root':root, 'evs':evs}


'''

# testing eigen_vector_3_symmetric_batch__given_eigen_value

import numpy as N
v = N.random.random( (5, 5, 5) )

import tomominer.filter.differential as FD
h = FD.hessian_3d(v)

m = FD.hessian_3d__max_magnitude(h)
ht = FD.hessian_3d__normalize(h, m)


import tomominer.linalg.eigen as LE
e = LE.eigen_value_3_symmetric_batch(ht)
c = LE.eigen_vector_3_symmetric_batch__given_eigen_value(ht, e[0])

c_t = LE.eigen_vector_3_symmetric_batch__simple(ht)['evs']

for i0 in range(c[0].shape[0]):
    for i1 in range(c[0].shape[1]):
        for i2 in range(c[0].shape[2]):
            c_t_v = c_t[0][i0][i1][i2]
            if c_t_v[0] < 0:        c_t_v = -c_t_v
            print i0, i1, i2, N.linalg.norm(N.array([c[0][i0,i1,i2], c[1][i0,i1,i2], c[2][i0,i1,i2]]) - c_t_v)



'''


