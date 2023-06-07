#!/usr/bin/env python

# for each cluster, random subsample subtomograms with a randomly selected proportion

if __name__ == '__main__':

    import json
    with open('subtomogram_cluster_subsample__op.json') as f:   op = json.load(f)

    with open(op['cluster_info_file']) as f:        ci = json.load(f)
    ci = {_[op['field']]:_ for _ in ci}

    with open(op['input_data']) as f:   d = json.load(f)


    # collect records according to cluster labels
    from collections import defaultdict
    clusters = defaultdict(list)
    for r in d:     clusters[ r[op['field']] ].append(r)


    high_freq_pids = set()
    if "high_freq" in op:       high_freq_pids = set(   op['high_freq']['pdb_ids']    )         # sometimes we want some complexes to have high frequency


    import random
    import math

    d_t = []
    stat = []

    for l in clusters:
        if ('selected_pdb_ids' in op) and (ci[l]['pdb_id'] not in set(op['selected_pdb_ids'])):     continue

        clus_size = len(clusters[l]) - 1
        if 'max_num' in op:     clus_size = min(clus_size, op['max_num'])

        min_num = 0
        if ci[l]['pdb_id'] in high_freq_pids:       min_num = int(op['high_freq']['minimum_freq'] * clus_size)

        cs = random.sample(clusters[l], random.randint(min_num, clus_size)  )

        d_t.extend( cs )

        ci[l]['size'] = len(cs)
        stat.append(ci[l])

    with open(op['output_data'], 'w') as f:   json.dump(d_t, f, indent=2)

    with open(op['output_stat'], 'w') as f:   json.dump(stat, f, indent=2)



'''   
example usage
cat /home/rcf-47/mxu/ln/frequent_structure/data/out/method/imputation-classification/ribosome/subtomograms/40-0.05/info.json | python ~/ln/tomominer/tomominer/common/subtomogram_cluster_subsample.py cluster_label 500 > data_config.json 



to check frequencies use:



'''
