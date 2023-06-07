




'''
functions for fitting one dimension gaussian functions, using Levenberg-Marquardt algorithm

see
https://docs.google.com/document/d/1TL8rcoY7zUlmFka9opuEc26fPPP-4E4eOUOh1T5kkZM/edit#

~/ln/tomominer/tomominer/fitting/gaussian/one_dim.py

'''



import numpy as N
import warnings
import exceptions


'''
fit gaussian function with following form

$ f(x, a, c) = a \exp{\left(- { \frac{x^2 }{ 2 c^2} } \right)} $
$ \frac{\partial f(x,a,c) }{\partial a} = \exp{\left(- { \frac{x^2 }{ 2 c^2} } \right)} $
$  \frac{\partial f(x,a,c) }{\partial c} = a \exp{\left(- { \frac{x^2 }{ 2 c^2} } \right)} \left (-\frac{x^2 }{ 2} \right ) (-2) c^{-3} $


as for parameters, a damping factor lambda_t >= 10 looks can get pretty stable results

'''

def fit__zero_mean(x, y, a0=None, c0=None, tolerance=1e-3, lambda_t=10.0, max_iter_num=1000, c_bound_max=100, verbose=False):

    assert  len(x) == len(y)


    i = N.argsort(x)
    x = x[i]
    y = y[i]

    x_abs_max = N.abs(x).max()

   
    # following initial guesses does not work
    if a0 is None:        a0 = y.max()           # just an adhoc start, such estimation is more robust than using      a0 = y[0]
    if c0 is None:        c0 = N.median(x)      ;           warnings.warn("such initial guess of c0 is not stable", exceptions.Warning)     # just an adhoc start. 

    a = a0
    c = c0
    if verbose:     print 'a0',a0, '   ', 'c0',c0

    iter_n = 0
    e = None
    e0 = None
    e_old = None
    while True:

        # calculate residue
        e_old = e

        exp_t = N.exp(- N.square(x) / (2 * c * c))
        yp = a * exp_t
        e = N.sqrt(  N.square(y - yp).sum()   )

        if e0 is None:            e0 = e

        # check if we can stop
        if e_old is not None:
            e_rate = (N.abs(e - e_old) / e0)
            if e_rate < tolerance:            break

        
        # update a and c using LM
        ga = exp_t
        gc = a * exp_t * N.square(x) / (c**3)

        J = N.zeros((len(x),2))
        J[:,0] = ga
        J[:,1] = gc

        Jq = N.dot(J.T, J)
        d = N.linalg.solve((Jq + lambda_t * N.diag(Jq)), N.dot(J.T, y - yp))
        if verbose and (e_old is not None):   print '\r', 'd',d, '    ', 'a',a, '    ', 'c',c, '    ', 'e',e, '    ', 'e_old',e_old, '    ', 'e0',e0, '    ', 'e_rate',e_rate, '           ',

        a += d[0]
        c += d[1]

        if N.abs(c) > (x_abs_max * c_bound_max):        return {'a':a, 'c':N.abs(c), 'e':e}             # if c is becomming too large, then this is a bad try, we do not improve any more

        iter_n += 1
        if iter_n > max_iter_num:       break           # just in case the process does not converge....


    return {'a':a, 'c':N.abs(c), 'e':e}



def fit__zero_mean__gaussian_function(x, a, c):
    return a * N.exp(- N.square(x) / (2 * c * c))


'''
it seems the fit__zero_mean() often stuck at local

test show that this gives pretty stable results

'''

def fit__zero_mean__multi_start(x, y, tolerance=1e-5, lambda_t=100.0, verbose=False):

    c0_s = N.linspace(x.min(), x.max(), num=10)

    best = None
    for i in range(1, len(c0_s)-1):
        try:
            p = fit__zero_mean(x, y, a0=y.max(), c0=c0_s[i], tolerance=tolerance, lambda_t=lambda_t, verbose=verbose)

            if (best is None) or (best['e'] > p['e']):            best = p
        except N.linalg.LinAlgError:
            pass

    return best





'''
# test code

import numpy as N
a = 3.0
c = 10.0
x = N.random.uniform(high=c*3.0, size=20)
y = a * N.exp(  - N.square(x) / (2 *c*c)  )
y += N.random.normal(scale=c/10.0)
y -= y.min()            # just make sure all values are non-negative


from tomominer.fitting.gaussian.one_dim import fit__zero_mean, fit__zero_mean__multi_start

p = fit__zero_mean(x=x, y=y, tolerance=1e-3, lambda_t=10.0, verbose=False)
print p['a'], p['c']

p = fit__zero_mean__multi_start(x=x, y=y, tolerance=1e-3, lambda_t=10.0, verbose=False)
print p['a'], p['c']


'''

 
