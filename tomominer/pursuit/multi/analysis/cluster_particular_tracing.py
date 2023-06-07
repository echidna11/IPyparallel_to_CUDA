#!/usr/bin/env python

'''

use this code for diagnosis, given a particularlly interesting cluster, see what are the other clusters that have large overlap with this cluster, and print out the file name of their averages to see how the averages look like.

notice: in the input cluster_particular_tracing__op.json, the cluster is indexed by populated cluster. But in the output file cluster_particular_tracing__out.json, the clusters are indexed by their selection

'''


import cPickle as pickle
import json
from collections import defaultdict

import numpy as N

def main(op):

    with open(op['cluster info file'], 'rb') as f:       ci = pickle.load(f)
    with open(op['fsc stat file']) as f:            fs = json.load(f)
    fs = {int(_):fs[_] for _ in fs}

    s = set([_['subtomogram'] for _ in ci[op['pass_i']][op['cluster']]['data_json']])       # collect subtomograms of the cluster of interest

    stat = defaultdict(dict)
    for pass_i in range(N.max([_ for _ in fs]) + 1):
        if pass_i not in fs:        continue
        for c in fs[pass_i]:
            fst = fs[pass_i][c]
            cit = ci[fst['pass_i']][fst['cluster']]
            st = set([_['subtomogram'] for _ in cit['data_json']])
            
            stat[pass_i][c] = {'pass_i': fst['pass_i'], 'cluster': fst['cluster'], 'is_specific': cit['is_specific'], 'intersection':len(s & st), 'ratio':float(len(s & st))/len(s | st)}

    with open('cluster_particular_tracing__out.json', 'w') as f:        json.dump(stat, f, indent=2)

if __name__ == '__main__':


    with open('cluster_particular_tracing__op.json') as f:      op = json.load(f)
    main(op)

