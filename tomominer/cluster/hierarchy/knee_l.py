
'''
calculate knee using the L method, see Salvador04
'''

import numpy as N



'''
the L method. 
parameters:         nc: the numbers of clusters, in ascending order (does not have to be consecqutive!!!)              em: evaluation metric
'''

def l_knee(nc, em):
    bi = len(nc)
    assert len(em) == bi

    best = None
    for ci in range(2, (bi-2)+1):
        rmse_l = regression_rmse(x=nc[:ci], y=em[:ci])
        rmse_r = regression_rmse(x=nc[ci:], y=em[ci:])

        c = nc[ci-1]
        b = nc[-1]

        rmse_c = ((c-1) / float(b - 1)) * rmse_l         +       ((b-c) / float(b-1)) * rmse_r

        if (best is None) or (rmse_c < best['rmse_c']):            best = {'rmse_c':rmse_c, 'rmse_l':rmse_l, 'rmse_r':rmse_r, 'ci':ci, 'c':c, 'b':b}

    return best



from sklearn.linear_model import LinearRegression
def regression_rmse(x, y):
    lr = LinearRegression()

    x = N.reshape(N.array(x), (-1,1))
    lr.fit(X=x, y=y)
    yp = lr.predict(x)

    return N.sqrt(N.square(yp-y).mean())            # RMSE



'''
refinement according to Figure 5 of Salvador04
parameters:         nc: the numbers of clusters, in ascending order (does not have to be consecqutive!!!)              em: evaluation metric
'''
def refine(nc, em, cut_off_min=6):
    print 'knee_l.refine()'
    nc = N.array(nc)
    em = N.array(em)

    assert len(nc) == len(em)
    for i in range(1, len(nc)):             assert nc[i-1] < nc[i]          # correctness check of monotonicity

    cut_off = nc[-1]
    current_knee = nc[-1]

    while True:
        last_knee = current_knee
        current_best = l_knee(nc[nc <= cut_off], em[nc <= cut_off])
        print current_best
        current_knee = current_best['c']
        cut_off = current_knee * 2
        if N.sum(nc <= cut_off) < cut_off_min:       break           # we require a minimum number of cutoffs to evaluate, in our application, since sometimes the real number of classes can be very small, we set a small cut_off_min, say 6
        if current_knee >= last_knee:        break

    return current_best




'''
given output of cluster.hierarchy.hierarchy.clusters_vs_dist_cutoff(), record the different cluster numbers given different cutoffs, this is used for calculate knee in knee_l.py
parameter:  lr: output of clusters_vs_dist_cutoff()
'''
def cluster_number_vs_dist_cutoff(lr):
    nc = []         # the numbers of clusters, in ascending order (does not have to be consecqutive!!!)
    md = []         # max_dist cutoff, i.e. a kind of evaluation metric
    for r in lr:
        nc.append(r['cluster_num'])
        md.append(r['dist_cutoff'])

    nc = N.array(nc)
    md = N.array(md)

    nci = N.argsort(nc)
    nc = nc[nci]
    md = md[nci]

    return {'nc':nc, 'md':md}



