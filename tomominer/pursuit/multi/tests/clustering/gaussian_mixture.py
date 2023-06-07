
data_dir = '/home/rcf-47/mxu/ln/frequent_structure/data/out/method/imputation-classification/multi/classification/30/0.005/0170/classify/pass_006'

import os
import pickle

with open(os.path.join(data_dir, 'dimension_reduction.pickle'), 'rb') as f:     dr = pickle.load(f)


x = dr['cfp_re']['red']
print x.shape


from sklearn.mixture import GMM

covariance_type='spherical'
#covariance_type='diag'
#covariance_type='full'

bic_s = []
aic_s = []
for n_components in range(1, 60):
    g = GMM(n_components=n_components, covariance_type=covariance_type)
    f = g.fit(x)
    bic = g.bic(x)
    bic_s.append(bic)
    aic = g.aic(x)
    aic_s.append(aic)
    print 'n_components', n_components, 'bic',bic, 'aic', aic


