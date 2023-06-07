#!/usr/bin/env python



'''

across multiple multiple pursuit experiments,
given the corresponding file stat file,
collect all FSC sum score

this is used to analysis background distribution of fsc sum scores


~/ln/tomominer/tomominer/pursuit/multi/analysis/statistics/collect_fsc_across_multiple_experiments__mpp_all.py

'''



import json
import cPickle as pickle

import numpy as N


import tomominer.statistics.vol as SV




def collect_pursuit(fs, pass_i=None):
    ps = fs['passes']
    ps = {int(_):ps[_] for _ in ps}

    r = {}      # return info
    for pass_i in ps:        
        with open(ps[pass_i]['cluster_info_file']) as f:       ci = pickle.load(f)

        r[pass_i] = {}
        for c in ci:
            r[pass_i][c] = {'fsc_sum': ci[c]['fsc'].sum()}

    return r



def main():

    with open('collect_fsc_across_multiple_experiments__mpp_all__op.json') as f:         op = json.load(f)

    with open(op['file stats']) as f:         fs_s = json.load(f)

    co = {}
    for fs_fn in fs_s:
        print fs_fn
        with open(fs_fn['file']) as f:      fs = json.load(f)
        if 'best_pass_i' in fs:     continue            # in this case, fs is for GA averaging
        co[fs_fn['file']] = collect_pursuit(fs, pass_i=(fs_fn['pass_i'] if 'pass_i' in fs_fn else None))            # in this case, fs if for pursuit
    
    with open(op['collect out'], 'wb') as f:     pickle.dump(co, f, protocol=-1)


if __name__ == '__main__':
    main()



'''
related code

~/ln/tomominer/tomominer/pursuit/multi/analysis/statistics/collect_fsc_across_multiple_experiments.py

'''



'''
# R commands for plotting

z = rnorm(1000)
hist(z, prob=T, main="Histogram With Fitted Density Curve, bw=.5")
lines(density(z, bw=.5), col='yellow', lwd=2)
points(0,0, col='blue') 


'''




