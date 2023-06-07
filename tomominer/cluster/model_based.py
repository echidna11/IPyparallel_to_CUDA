
# functions for model based clustering

import numpy as np

# model based clustering using spherical co-variance matrix
# important, unlike MATLAB, the cluster labels in python starts from 0!!
def spherical_model(op):
    x = op['x']
    n = x.shape[0]      # data point number
    d = x.shape[1]      # dimension

    G = op['G']            # cluster number
    
    lbl_init = op['lbl_init']     # initialize initial cluster labels, which can be obtained from random assignment or kmeans clustering
    assert(lbl_init.size == n)

    z_init = np.zeros([n, G])
    for clus_i in range(G):
        z_init[(lbl_init == clus_i), clus_i] = 1.0

    
    tau_k = None
    mu_k = None
    lambda_k = None

    log_lik_model_initial = None
    log_lik_model_old = None


    # iterative EM
    for iter_i in range(op['max_iteration_num']):
        #--------------------------------------
        # E step
        if tau_k is not None:
            z = np.zeros([n, G]) + float('NaN')
            for k in range(G):
                z[:,k] = tau_k[k] * mvnpdf(x, mu_k[k,:], lambda_k[k]*np.eye(d))
            assert(all(np.isfinite(z.flatten())))
            
            z_sum = z.sum(axis=1)
            z /= np.tile(z_sum.reshape(-1,1), [1, G])
            z[np.logical_not(np.isfinite(z))] = 0.0
        else:
            z = z_init
        
        
 

        # -----------------------------------
        # M step
        n_k = z.sum(axis=0)
        if any(n_k < 1e-10):
            return None
        
        tau_k = n_k / n_k.sum()
        
        mu_k = np.zeros([G,d]);
        for k in range(G):
            mu_k[k,:] = (np.tile(z[:,k].reshape(-1,1), [1, d]) * x).sum(axis=0) / n_k[k]        # .reshape(-1, 1) is used to convert to column vector
    
        lambda_k = np.zeros([G,1])
        for k in range(G):
            trace_t = 0
            for d_i in range(d):
                trace_t = trace_t + ( ((x[:,d_i] - mu_k[k, d_i]) ** 2) * z[:, k]).sum();
                assert(np.isfinite(trace_t))
            
            lambda_k[k] = trace_t / (d * n_k[k]);
            
            assert(np.isfinite(lambda_k[k]))
        
        lambda_k_small_cutoff = 1e-10
        if any(lambda_k <= lambda_k_small_cutoff) and any(lambda_k > lambda_k_small_cutoff):
            # occationally there may be some lambda == 0
            # just use following to prevent numberical problem
            lambda_k[lambda_k <= lambda_k_small_cutoff] = min(lambda_k[lambda_k > lambda_k_small_cutoff]) / 10;
        
        
        log_lik_complete = spherical_model__log_complete_likelihood(x, z, tau_k, mu_k, lambda_k);          assert(np.isfinite(log_lik_complete));
        log_lik_model = spherical_model__log_likelihood(x, tau_k, mu_k, lambda_k);     assert(np.isfinite(log_lik_model));
        bic = spherical_model__BIC(n, d, G, log_lik_model)
        #print '%d  \t %g \t %g \t %g '%(iter_i, log_lik_complete,  log_lik_model, bic)
        
        
        if log_lik_model_initial is None:
            log_lik_model_initial = log_lik_model
        
        if log_lik_model_old is not None:
            if abs(log_lik_model - log_lik_model_old) < (abs(log_lik_model_initial) * op['stop_tolorence']):
                break;
        
        log_lik_model_old = log_lik_model

    
    labels = np.argmin(z, axis=1)

    return {'labels':labels, 'z':z, 'mu':mu_k, 'lambda':lambda_k, 'iter_i':iter_i, 'bic':bic}





def mvnpdf(x, mu, sigma):
    if np.linalg.det(sigma) == 0:
        return None 
    n = x.shape[0]      # data point number
    k = x.shape[1]      # dimension number
    dev = x - np.tile(mu, (n, 1))
    inv_sigma = np.linalg.inv(sigma)
    desc = ( dev.dot(inv_sigma) * dev ).sum(axis=1)
    c0 = (2*np.pi)**(-0.5*k)
    c1 = np.linalg.det(sigma) ** (-0.5*k)
    return c0*c1*np.exp(-0.5*desc)





def spherical_model__log_complete_likelihood(x, z, tau_k, mu_k, lambda_k):
    n = x.shape[0]      # data point number
    d = x.shape[1]      # data dimension number
    G = z.shape[1]      # cluster number
   
    assert(z.shape == (n, G))
    tau_k = tau_k.reshape(1, -1);   assert(tau_k.shape == (1, G))
    assert(mu_k.shape == (G, d))
    assert(lambda_k.shape == (G, 1))
    
    
    lik_m = z * np.tile(np.log(tau_k), [n, 1])
    lik_m = lik_m + z
    
    lik_gau = np.zeros([n, G]);
    for k in range(G):
        lik_gau[:, k] = -(1/2) * (  (x - np.tile(mu_k[k,:], [n, 1])) ** 2  ).sum(axis=1) / (lambda_k[k]**3) - (3/2)*(np.log(2 * np.pi * lambda_k[k]));
    
    lik_m = lik_m + z * lik_gau
    
    lik = lik_m.sum()
    return lik


def spherical_model__log_likelihood(x, tau_k, mu_k, lambda_k):
    
    assert(all(np.isfinite(x.flatten())))
    assert(all(np.isfinite(tau_k.flatten())))
    assert(all(np.isfinite(mu_k.flatten())))
    assert(all(lambda_k.flatten() > 0))
    
    n = x.shape[0]     # data point number
    d = x.shape[1]     # data dimension number
    G = len(tau_k)     # cluster number
    
    assert(tau_k.size == G)
    assert(mu_k.shape == (G, d))
    assert(lambda_k.shape == (G, 1))

    sum_v = np.zeros([n, 1]);
    for k in range(G):
        sum_v = sum_v + tau_k[k] * mvnpdf(x, mu_k[k,:], lambda_k[k] * np.eye(d)).reshape(-1,1)
        assert(all(sum_v.flatten() >= 0));        assert(all(np.isfinite(sum_v.flatten())));
    
    small_threshold = 1e-20
    sum_v[sum_v < small_threshold] = small_threshold
    
    lik = sum(np.log(sum_v))
    return lik




# calculate bayesian information criterion given a cluster
def spherical_model__BIC(n, d, G, log_lik):
    # see paper    Gaussian Parsimonious Clustering Models,   Table 1
    m_m = G + (G*d+G-1)
    bic = 2*log_lik - m_m * np.log(n)
    return float(bic)


import matplotlib.pyplot as mpl
def spherical_model__test():
    print 'spherical_model__test()'
    n2 = 100
    x1 = np.random.multivariate_normal([0,0], 1*np.eye(2), n2)
    x2 = np.random.multivariate_normal([10,10], 2*np.eye(2), n2)
    x3 = np.random.multivariate_normal([-10,10], 3*np.eye(2), n2)
    x = np.vstack((x1, x2, x3))
    
    if False:
        print x
        import matplotlib.pyplot as mpl
        mpl.plot(x[:,0], x[:,1], 'x')
        mpl.show()

        print 'show done'
    G = 3   
    lbl_init = np.floor(np.random.rand(x.shape[0]) * G) 
    re = spherical_model({'x':x, 'lbl_init':lbl_init, 'G':G, 'max_iteration_num':1000, 'stop_tolorence':1e-5})
    print re['mu']
    print re['lambda']
    print re['iter_i']



# testing
if __name__ == '__main__':

    spherical_model__test()


