#!/usr/bin/env python






# calculate contegincy table for each cluster
# find its largest overlaping true cluster, 
# then see if the true cluster is fully covered by this cluster

# see if the predicted cluster also can correctly seperate other clusters

# perform pairwise alignment and calculate structual consistency, in terms of cross correlation, and fourier shell correlation, list together with contegincy table

# calculate and report overall cluster statistics using scikit-learn

import sys
import os
import json
import tomominer.io.file as IV
import tomominer.statistics.vol as SV
import numpy as N
import tomominer.core as core
import pickle


def fsc_stat(pa, voxel_spacing=1.0):

    from collections import defaultdict
    record = defaultdict(dict)
    for a in pa:
        v1k = a['v1_key']['subtomogram']
        v2k = a['v2_key']['subtomogram']

        v1 = IV.get_mrc(v1k)
        v2 = IV.get_mrc(v2k)

        v2r = core.rotate_vol_pad_mean(v2, a['angle'], a['loc'])

        IV.put_mrc(v2r, ('%d-%d.mrc'%(a['v1_key']['cluster_id'], a['v2_key']['cluster_id'])) )

        resolution, i = SV.resolution(v1, v2r, voxel_spacing=1.0, cutoff=0.5)

        id1 = a['v1_key']['cluster_id']
        id2 = a['v2_key']['cluster_id']

        assert id2 not in record[id1]
        record[id1][id2] = {}

        record[id1][id2]['v1k'] = v1k
        record[id1][id2]['v2k'] = v2k
        record[id1][id2]['resolution'] = resolution
        record[id1][id2]['fsc'] = SV.fsc(v1, v2r)
        record[id1][id2]['cor'] = N.corrcoef(v1.flatten(), v2r.flatten())[0,1]


    return record


if __name__ == '__main__':
    
    op_file = sys.argv[1]
    with open(op_file) as f:    op = json.load(f)

    with open(op['true_member']) as f:  true_member = json.load(f)
    with open(op['pred_member']) as f:  pred_member = json.load(f)

    count = N.zeros([max([int(_['cluster_label']) for _ in true_member])+1, max([int(_['cluster_label']) for _ in pred_member])+1], dtype=N.int)

    true_lbl = {_['subtomogram']:int(_['cluster_label']) for _ in true_member}
    pred_lbl = {_['subtomogram']:int(_['cluster_label']) for _ in pred_member}

    for s in pred_lbl:
        count[true_lbl[s], pred_lbl[s]] += 1



    with open(op['true_avg'], 'r') as f:        v1ks_v = json.load(f)

    v1ks = []
    for r in v1ks_v:  
        r['subtomogram'] = r['ground_truth']
        r['cluster_id'] = int(r['cluster_label'])
        v1ks.append(r)

    pid = {_['cluster_id']:str(_['pdb_id']) for _ in v1ks }



    with open(op['cluster_avg'], 'rb') as f:            v2ks_d = pickle.load(f)
    v2ks_d = v2ks_d['template_keys']

    v2ks = []
    for c in v2ks_d:
        r = v2ks_d[c]
        r['cluster_id'] = c
        v2ks.append(r)



    from tomominer.parallel.queue_master import QueueMaster
    runner = QueueMaster('localhost', 5011)

    align_op = {'with_missing_wedge':True, 'L':36}
    import classify.util as CU
    pa = CU.pairwise_align_keys__global(runner, v1ks, v2ks, op=align_op)

    fs = fsc_stat(pa)



    # print a table of statistics
    # IMPORTANT: because generated using MATLAB, true class label starts from 1

    with open('stat.out', 'w') as f:
        for ct in range(count.shape[0]):
            if ct not in pid:   continue

            f.write('%d\t%s\t%d\t\t'%( ct, pid[ct], count[ct, :].sum() ))
            for cp in range(count.shape[1]):
                f.write('\t')
                if cp in fs[ct]:    f.write( '%d, %.1f, %.2f'%(count[ct][cp], fs[ct][cp]['resolution'], fs[ct][cp]['cor']) )
            f.write('\n')





'''

example usage

~/ln/tomominer/tomominer/classify/analysis/cluster_consistency_with_ground_truth.py ./config.json


'''


