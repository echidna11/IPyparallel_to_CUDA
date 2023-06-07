#!/usr/bin/env python



'''
assume only one refernece is the true one, and others are controls, select subtomograms that only best matches the reference

~/ln/tomominer/tomominer/template/search/fast_align/reference/filter_/single_true.py

'''

import json, copy

import numpy as N

def main():

    with open('template_search_fast_align_reference_filter_single_true__op.json') as f:     op = json.load(f)

    with open(op['align in']) as f:     al = json.load(f)

    with open(op['data in']) as f:      dj = json.load(f)

    assert  len(al) == len(dj)

    dj_new = []
    for ai, a in enumerate(al):
        a = {int(_):a[_] for _ in a}
        
        score_true = a[op['true reference ind']]['score']

        scores_control = []
        for ri in a:
            if ri == op['true reference ind']:  continue
            scores_control.append(   a[ri]['score']   )
        scores_control = N.array(scores_control)

        if score_true <= max(scores_control):       continue

        if ('std threshold' in op) and (score_true <= (scores_control.mean() + float(op['std threshold']) * scores_control.std())):      continue

        dt = copy.deepcopy(dj[ai])
        dt['score'] = a[op['true reference ind']]['score']
        dt['loc'] = a[op['true reference ind']]['loc']
        dt['angle'] = a[op['true reference ind']]['angle']

        dj_new.append(dt)

    print len(dj_new), 'out of', len(dj), 'subtomograms selected'
    with open(op['data out'], 'w') as f:   json.dump(dj_new, f, indent=2)



if __name__ == '__main__':
    main()

