#!/usr/bin/env python


'''

given template search against a template of interest, and several control templates (can be fake templates), select only those subtomograms that achieve highest template matching score. 
in addition, a adhoc way, can also calculate mean and standard deviation of the control scores to set a score cutoff


~/ln/tomominer/tomominer/template/search/filter/max_against_controls.py

'''

import json
import numpy as N

def main():
    with open('template_search_filter_max_against_controls__op.json') as f:     op = json.load(f)

    with open(op['true match']) as f:       djt = json.load(f)
    djt = {_['subtomogram']:_ for _ in djt}

    
    djc = {}
    for c in op['control matches']:
        with open(c) as f:     djc_t = json.load(f)
        djc[c] = {_['subtomogram']:_ for _ in djc_t}
        del djc_t


    djtf = []
    for k in djt:
        score_true = djt[k]['score']
        score_control = N.array([djc[_][k]['score'] for _ in djc])

        if score_true <= score_control.max():       continue

        if 'std threshold' in op:
            if score_true <= (score_control.mean() + op['std threshold']*score_control.std()):  continue

        djtf.append(djt[k])


    print 'result match number', len(djtf)

    with open(op['out file'], 'w') as f:        json.dump(djtf, f, indent=2)


if __name__ == '__main__':
    main()


