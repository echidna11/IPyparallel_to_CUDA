#!/usr/bin/env python

# given ground truth label (of simulation data), and predicted label from tomominer.pick.pose_normalize_cluster_kmeans, calculate the membership consistancy between ground truth and prediction

if __name__ == '__main__':

    import json
    with open('pose_normalize_cluster_kmeans__member_consistancy__op.json') as f:       op = json.load(f)

    
    with open(op['classification op file']) as f:   cls_op = json.load(f)

    with open(cls_op['data_json_file']) as f:   dj_true = json.load(f)
    dj_true = {_['subtomogram']:_ for _ in dj_true}

    with open(op['classification out file']) as f:     dj_pred = json.load(f)

    cts = list(range(max(dj_true[_]['cluster_label'] for _ in dj_true)+1))
    cps = list(range(max(_['cluster_label'] for _ in dj_pred)+1))

    import numpy as N
    cnt = N.zeros([len(cts), len(cps)], dtype=int)

    for dj_pred_t in dj_pred:        cnt[dj_true[dj_pred_t['subtomogram']]['cluster_label'], dj_pred_t['cluster_label']] += 1

    from  pursuit.multi.analysis.cluster_avg_high_fsc_consistency_with_ground_truth import stat_table__arrange__ind

    i = stat_table__arrange__ind(cnt=cnt, cts=cts, cps=cps)


    with open('pose_normalize_cluster_kmeans__member_consistancy__out.txt', 'w') as f:
        for it in i['cts']:
            for ip in i['cps']:
                print >>f, cnt[it, ip], '\t',
            print >>f

