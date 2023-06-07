#!/usr/bin/env python


'''

over iterations, calculate statistics of membership consistency between preticted patterns and ground truth

~/ln/tomominer/tomominer/pursuit/multi/analysis/iteration/pattern_ground_truth_member_consistency.py

'''

import sys, json
import cPickle as pickle

import numpy as N
import scipy.stats as SS

def collect_djs(fs):
    print 'collect_djs()'

    ps = fs['passes']
    ps = {int(_):ps[_] for _ in ps}

    djs = {}
    for pi in ps:
        print pi, '     \r',         ;       sys.stdout.flush()

        djs[pi] = {}

        with open(ps[pi]['cluster_average_select_file'], 'rb') as f:     cas = pickle.load(f)

        st = cas['selected_templates']
        for si in st:
            assert      st[si]['id'] == si
            k = st[si]['subtomogram']

            ti = cas['tk_info']
            assert      ti[k]['pass_i'] == st[si]['pass_i']
            assert      ti[k]['cluster'] ==st[si]['cluster']
            djs[pi][si] = ti[k]['data_json']
    

    return djs 


'''
calculate contingency table over iterations

'''
def contingency_count(lbl_true, djp):
    print   'contingency_count()'

    cont = {}          # contingency table c[pass_i][true_label][pred_label]
    for pi in djp:
        cont[pi] = {}
        for c in djp[pi]:
            for d in djp[pi][c]:
                k = d['subtomogram']
                l = lbl_true[k]
                if l not in cont[pi]:       cont[pi][l] = {}
                if c not in cont[pi][l]:    cont[pi][l][c] = 0
                cont[pi][l][c] += 1

    return cont




def contingency_tab(cont, remove_zero=True):

    tabs = {}

    for p in cont:
        lbl_pred = []
        for l in cont[p]:      lbl_pred.extend(cont[p][l].keys())

        lbl_true_max = N.max(cont[p].keys())
        lbl_pred_max = N.max(lbl_pred)

        cont_t = N.zeros((lbl_true_max+1, lbl_pred_max+1))
        for l in cont[p]:
            for c in cont[p][l]:
                cont_t[l][c] = cont[p][l][c]

        tabs[p] = cont_t


    if remove_zero:
        for p in tabs:
            t = tabs[p]
            t = t[t.sum(axis=1)>0, :]
            t = t[:, t.sum(axis=0)>0]
            tabs[p] = t

    return tabs





def main():
    with open('pattern_ground_truth_member_consistency__op.json') as f:     op = json.load(f)

    with open(op['data true']) as f:     djt = json.load(f)
    lbl_true = {_['subtomogram']:_['cluster_label'] for _ in djt}


    with open(op['file stat']) as f:        fs = json.load(f)
    djp = collect_djs(fs)

    cont = contingency_count(lbl_true=lbl_true, djp=djp)
    cont_tab = contingency_tab(cont)

    print 'chi2_contingency test'
    print 'pass\tstatistics\tp-val\tdegree_of_freedom'

    cs_chi2 = {_:SS.chi2_contingency(cont_tab[_]) for _ in cont_tab}
    for p in cs_chi2:        print p, '\t', cs_chi2[p][0], '\t', cs_chi2[p][1], '\t', cs_chi2[p][2]

    print 'power_divergence test'
    print 'pass\tstatistics\tp-val'
    cs_pd = {_:SS.power_divergence(cont_tab[_], axis=None) for _ in cont_tab}
    for p in cs_pd:        print p, '\t', cs_pd[p][0], '\t', cs_pd[p][1]




if __name__ == '__main__':
    main()




'''
related code


http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.chi2_contingency.html

'''

