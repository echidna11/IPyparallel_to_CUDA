#!/usr/bin/env python

'''
given a collection of json files, and they are aligned

then randomly choose a pair of subtomograms and calculate correlation, then convert this correlation to SNR according to Jachohim Frank's paper, and use matlab to display the histogram of the SNRs

~/ln/tomominer/tomominer/average/genetic_algorithm/analysis/instance_pair_snr.py
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

    with open(op['pmpg_file']) as f:       pmpg = pickle.load(f)


    dj = pmpg['dj']
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
related code and analysis

~/ln/tomominer/tomominer/average/genetic_algorithm/analysis/instance_pair_snr.py
'''




