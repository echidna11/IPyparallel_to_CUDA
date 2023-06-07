#!/usr/bin/env python


# given classification iterations, inspect the cluster average structural similarity between consective iterations

def pairwise_structure_similarity(d0, d1):
    import tomominer.io.file as IV
    import numpy as N

    from collections import defaultdict

    n0 = max( [_ for _ in d0] ) + 1
    n1 = max( [_ for _ in d1] ) + 1

    s = N.zeros( (n0, n1) ) + float('NaN')
    for i0 in d0:
        for i1 in d1:
            v0 = IV.get_mrc(d0[i0]['subtomogram'])
            v1 = IV.get_mrc(d1[i1]['subtomogram'])

            c = N.corrcoef(v0.flatten(), v1.flatten())[0,1]
            
            s[i0, i1] = c


    return s


if __name__ == '__main__':

    import os
    import pickle

    pass_i = 0
    while True:
        cp0f  = os.path.join(os.getcwd(), 'pass_%03d' % (pass_i), 'cluster_averaging.pickle')
        cp1f  = os.path.join(os.getcwd(), 'pass_%03d' % (pass_i + 1), 'cluster_averaging.pickle')

        if not os.path.exists(cp0f) :       break
        if not os.path.exists(cp1f) :       break

        with open(cp0f) as f:   cp0 = pickle.load(f)
        with open(cp1f) as f:   cp1 = pickle.load(f)

        #import pdb;     pdb.set_trace()

        print pass_i, pairwise_structure_similarity(cp0['template_keys'], cp1['template_keys']).max(axis=0).mean()

        pass_i += 1

