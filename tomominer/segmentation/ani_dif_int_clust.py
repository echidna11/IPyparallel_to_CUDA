
'''
functions mainly used for finding out structural region in order to segment out target complexes

first perform anistropic diffusion filtering (apply gaussian filtering if necessary). 
Then use either kmeans or gaussian mixture to seperate intensity into three / five clusters. 
The top cluster will define a structural region.
'''

import numpy as N
from sklearn.mixture.gmm import GMM
from sklearn.cluster import KMeans

import tomominer.filter.anistropic_diffusion.util as FADU

def do_segmentation(v, op):

    # perform anistropic diffusion filtering
    dif_op = op['anistropic diffusion filtering']
    vf = FADU.anisodiff3(v, niter=dif_op['niter'], kappa=dif_op['kappa'], gamma=dif_op['gamma'], gauss_smooth_sigma=(dif_op['gauss_smooth_sigma'] if 'gauss_smooth_sigma' in dif_op else None), option=dif_op['option'])
    
    # perform gaussian mixture model seperation
    vfv = N.reshape(vf, (-1,1))     # change into a N by 1 matrix
    vfv /= vfv.std()        # intensities must be rescaled. If the values are too small, GMM does not work

    if 'mixture model' in op:
        assert      'kmeans' not in op
        g = GMM(n_components=op['mixture model']['n_components'], n_init=op['mixture model']['n_init'])
        g.fit(vfv)

        vfp = g.predict(vfv)
    elif 'kmeans' in op:
        assert      'mixture model' not in op
        km = KMeans(n_clusters=op['kmeans']['n_clusters'])
        vfp = km.fit_predict(vfv)

    vfp = N.reshape(vfp, vf.shape)

    vf_mean = N.zeros(vfp.shape) + N.nan
    means = N.zeros(vfp.max()+1) + N.nan
    for l in range(vfp.max()+1):
        if (vfp == l).sum() == 0:       continue
        means[l] = vf[vfp == l].mean()
        vf_mean[vfp == l] = means[l]

    vfp_order = N.argsort(means)
    lbl = N.zeros(vf.shape)
    for i,o in enumerate(vfp_order):
        if not N.isfinite(means[o]):    continue
        lbl[vfp==o] = i

    return {'vf':vf, 'vf_mean':vf_mean, 'lbl':lbl}

