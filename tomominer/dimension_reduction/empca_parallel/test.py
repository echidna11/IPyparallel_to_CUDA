#!/usr/bin/env python


# adapted from sparse_pca_test.py
# perform profiling of empca, see how to make some of the components parallel to reduce computation cost


'''

for profiling

python -m cProfile -o /tmp/p ~/ln/tomominer/tomominer/dimension_reduction/empca_parallel/test.py


import pstats
p = pstats.Stats('/tmp/p')
p.strip_dirs().sort_stats('cumulative').print_stats(30)



observations: the most time consuming functions are
solve_eigenvectors()
R2vec()
dot()
outer()

'''



if __name__ == '__main__':



    if False:
        n_sample = 1000
        n_dim = 10
        n_rand_dim = 10000
    elif False:
        n_sample = 100
        n_dim = 10
        n_rand_dim = 10000
    elif True:
        n_sample = 20000
        n_dim = 10
        n_rand_dim = 10000
    elif False:
        n_sample = 100
        n_dim = 10
        n_rand_dim = 100



    import numpy as np
    if False:
        var = np.diag([n_dim-i for i in range(n_dim)])
        # simply generate a random rotation and to apply to the diagonal covarience matrix var
        svd_u, svd_s, svd_v = np.linalg.svd(np.random.random((n_dim, n_dim)))
        var = svd_u.dot(var)
    else:
        var = np.ones((n_dim,n_dim))


    x0 = np.random.multivariate_normal(np.zeros(n_dim), var, size=n_sample)
    x1 = np.random.multivariate_normal(np.zeros(n_dim), var, size=n_sample)
    xn = np.random.normal(scale=np.diag(var).min(), size=(n_sample, n_rand_dim))
    x = np.hstack((x0, x1, xn))
    print x.shape

    w = np.isfinite(x)
    x1 = np.copy(x)
    x1[np.logical_not(w)] = 0

    if False:
        import tomominer.dimension_reduction.empca_parallel.empca as DEE
    else:
        import tomominer.dimension_reduction.empca as DEE


    import time
    cur_time = time.time()
    m = DEE.empca(data=x1, weights=w, nvec=50, niter=2)
    print '%.2f sec'%(time.time() - cur_time)

    if True:
        dre = m.coeff       # m.coeff is approximately   np.dot(m.model, m.eigvec.T)
    else:
        dre = np.dot(m.model, m.eigvec.T)


    d = np.dot(x1, m.eigvec.T) - m.coeff

