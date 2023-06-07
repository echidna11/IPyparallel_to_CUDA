

'''

inspect the scaling properties of hierarchical clustering

'''



data_dir = '/home/rcf-47/mxu/ln/frequent_structure/data/out/method/imputation-classification/multi/classification/30/0.005/0130/classify/pass_029'

import os, sys, pickle
import numpy as N


with open(os.path.join(data_dir, 'dimension_reduction.pickle'), 'rb') as f:     dr = pickle.load(f)

with open(os.path.join(data_dir, 'cluster-hierarchy.pickle'), 'rb') as f:     ch = pickle.load(f)


hi = ch['hc_re']['info']

size_s = []
dist_s = []
fsc_s = []

for i in hi:
    hit = hi[i]
    if 'ssnr' not in hit:   continue
    size_s.append(hit['size'])
    dist_s.append(hit['dist'])
    fsc_s.append(hit['ssnr']['fsc'].sum())

import scipy.io as SI
SI.savemat(     '/tmp/t.mat', {  'size_s':size_s, 'dist_s':dist_s, 'fsc_s':fsc_s, 'dim_n':int(dr['cfp_re']['red'].shape[1])    }     )


'''
% matlab code

clear all
load /tmp/t

dim_n = double(dim_n);
size_s = single(size_s);

volume_ball = ( (pi ^ (dim_n/2)) / gamma(dim_n/2 + 1) ) * ( (dist_s/2).^dim_n );        % the volume of a ball in high dimension, given the radius (dist_s/2)

scatter3(log10(size_s), log10(volume_ball), fsc_s)
xlabel('log size')
ylabel('log volume')
zlabel('fsc')


scatter(log(size_s) - log(volume_ball), fsc_s)
xlabel('density')
ylabel('fsc')

'''



