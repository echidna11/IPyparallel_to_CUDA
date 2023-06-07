#!/usr/bin/env python


'''
select clusters
'''

import os, json

if __name__ == '__main__':
    with open('cluster_selection__op.json') as f:     op = json.load(f)

    with open(op['original data file']) as f:     dj = json.load(f)
    with open(op['label data file']) as f:     lbl = json.load(f)
    lbl = {_['subtomogram']:_['cluster_label'] for _ in lbl}

    lbl_sel = set(op['selected clusters'])

    djf = []
    for d in dj:
        if d['subtomogram'] not in lbl:     continue
        d['cluster_label'] = lbl[d['subtomogram']]
        if d['cluster_label'] not in lbl_sel:       continue
        djf.append(d)


    with open(op['output data file'], 'w') as f:     json.dump(djf, f, indent=2)


