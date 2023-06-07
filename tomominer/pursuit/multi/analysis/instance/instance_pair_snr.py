#!/usr/bin/env python

'''
given a cluster and its average, order these instances according to alignment scores, select top number of instances, align all its instances to the average
then randomly choose a pair of subtomograms and calculate correlation, then convert this correlation to SNR according to Jachohim Frank's paper, and use matlab to display the histogram of the SNRs

~/ln/tomominer/tomominer/pursuit/multi/analysis/instance/instance_pair_snr.py
'''

import os
import json
import cPickle as pickle
import numpy as N
import scipy.stats as SS
import tomominer.io.file as IF
import tomominer.geometry.rotate as GR


if __name__ == '__main__':
    
    with open('instance_pair_snr__op.json') as f:     op = json.load(f)
    with open(op['cluster_average_select__file'], 'rb') as f:   cs = pickle.load(f)

    st = cs['selected_templates']
    ti = cs['tk_info']

    tic = ti[st[op['cluster_id']]['subtomogram']]

    dj = tic['data_json']
    dj = sorted(dj, key=lambda _ : (-_['score']))

    if 'top_num' in op:     dj = dj[:op['top_num']]

    vrs = []
    
    for i, r in enumerate(dj):
        v = IF.get_mrc(r['subtomogram'])
        vrs.append( GR.rotate(v, angle=N.array(r['angle']), loc_r=N.array(r['loc']), default_val=v.mean()) )


    print 'sampling', op['pair_num'], 'pairs from', len(vrs), 'aligned subtomograms'

    with open(op['out_file'], 'w') as f:
        for i in range(op['pair_num']):
            while True:
                i0 = N.random.randint(len(dj))
                i1 = N.random.randint(len(dj))
                if i0 != i1:        break

            c, _ = SS.pearsonr(vrs[i0].flatten(), vrs[i1].flatten())        # pearson correlation between two aligned subtomograms
            s = c / (1 - c)         # SNR according to Equation 3.61 of book     frank Three-dimensional Electron Microscopy of Macromolecular Assemblies
            print >>f, s

'''

for analysis using R
t = read.table('instance_pair_snr__out.txt')
t = as.vector(t[,1])
mean(t)
sd(t)
pdf('instance_pair_snr__out__hist.pdf')
hist(t, main='SNR Distribution', xlab='SNR')
dev.off()
'''

'''
For analysis using matlab:
~/ln/tomominer/tomominer/pursuit/multi/analysis/instance_pair_snr.py

In matlab:
clear all;
load -ASCII /tmp/snr.mat
size(snr)
hist(snr, 100)
mean(snr)
std(snr)
'''
