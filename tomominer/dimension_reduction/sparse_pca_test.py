
# test the ability to select informative dimension from high dimension data, when there are very large amount of pure random variables

if __name__ == '__main__':


    if True:
        n_sample = 1000
        n_dim = 10
        n_rand_dim = 100000
    elif False:
        n_sample = 100
        n_dim = 10
        n_rand_dim = 10000
    elif False:
        n_sample = 30000
        n_dim = 10
        n_rand_dim = 500



    import numpy as N
    if False:
        var = N.diag([n_dim-i for i in range(n_dim)])
        # simply generate a random rotation and to apply to the diagonal covarience matrix var
        svd_u, svd_s, svd_v = N.linalg.svd(N.random.random((n_dim, n_dim)))
        var = svd_u.dot(var)
    else:
        var = N.ones((n_dim,n_dim))



    x0 = N.random.multivariate_normal(N.zeros(n_dim), var, size=n_sample)
    x1 = N.random.multivariate_normal(N.zeros(n_dim), var, size=n_sample)
    xn = N.random.normal(scale=N.diag(var).min(), size=(n_sample, n_rand_dim))
    x = N.hstack((x0, x1, xn))
    print x.shape


    if True:
        # randomly add some missing value to see if the code can handle missing values
        missing_val_num = 100
        row_i = N.int32( N.floor(N.random.random(missing_val_num) * n_sample) )
        col_i = N.int32( N.floor(N.random.random(missing_val_num) * n_dim) )
        for i in range(missing_val_num):    x[row_i, col_i] = float('nan')


    if True:

        import os
        import sys
        sys.path.append(os.path.join(os.getenv('HOME'), 'ln/tomominer/tomominer'))

        w = N.isfinite(x)
        x1 = N.copy(x)
        x1[N.logical_not(w)] = 0
        import tomominer.dimension_reduction.empca as drempca

        import time
        cur_time = time.time()
        m = drempca.empca(data=x1, weights=w, nvec=20)
        print '%.2f sec'%(time.time() - cur_time)

        if True:
            dre = m.coeff       # m.coeff is approximately   N.dot(m.model, m.eigvec.T)
        else:
            dre = N.dot(m.model, m.eigvec.T)



    if False:

        import multiprocessing
        from sklearn.decomposition import PCA

        sp = PCA(n_components=5)

        import time
        cur_time = time.time()
        drt=sp.fit(x)
        print '%.2f sec'%(time.time() - cur_time)


    if False:

        import multiprocessing
        from sklearn.decomposition import SparsePCA

        sp = SparsePCA(n_components=5, alpha=4, n_jobs=12, verbose=True)

        import time
        cur_time = time.time()
        drt=sp.fit(x)
        print '%.2f sec'%(time.time() - cur_time)




    if False:

        import multiprocessing
        from sklearn.decomposition import MiniBatchSparsePCA

        mbsp = MiniBatchSparsePCA(n_components=5, alpha=4, batch_size=1000, shuffle=True, n_jobs=8, verbose=True)
        import time
        cur_time = time.time()
        drt=mbsp.fit(x)
        print '%.2f sec'%(time.time() - cur_time)

    
    if False:

        import multiprocessing
        from sklearn.decomposition import RandomizedPCA

        rp = RandomizedPCA(n_components=10)
        import time
        cur_time = time.time()
        drt=rp.fit(x)
        print '%.2f sec'%(time.time() - cur_time)


    if False:

        import multiprocessing
        from sklearn.decomposition import DictionaryLearning

        dl = DictionaryLearning(n_components=10, alpha=4,  n_jobs=multiprocessing.cpu_count(), verbose=True)

        import time
        cur_time = time.time()
        drt=dl.fit(x)
        print '%.2f sec'%(time.time() - cur_time)





    print drt.components_.shape
    print drt.components_[:,:n_dim*3]
    print drt.components_[1][:]

    print N.dot(drt.components_[0], drt.components_[1])

    print sum(abs(drt.components_[0]) > 1) 


