#!/usr/bin/env python


'''

see
https://docs.google.com/document/d/1CcXs_gCoTPY8Cn0K3pY89cGMassjgVA8VZn_Yey0iAM/edit

also see 
~/ln/tomominer/tomominer/pursuit/multi/tests/clustering/refine/second_best_align_score_distribution.py

'''


import os, sys, json
import cPickle as pickle

from sklearn import svm
import numpy as N

def collect_items(template_id, al):
    # collect subtomograms / items that best align with this template
    alb = []
    for alt in al:
        if alt[template_id]['score'] < N.max([alt[_]['score'] for _ in alt]):            continue
        alb.append(alt)

    return alb

def get_best_two_scores(alt):
    s = [alt[_]['score'] for _ in alt]
    s = sorted(s, reverse=True)
    return (s[0], s[1])


'''
for each template T, collect subtomograms whose best alignment is with T
find the cutoff that best distinguish highest and second highest alignment scores
'''
def find_cutoff(al):

    cuts = {}
    for c in al[0].keys():
        alb = collect_items(c, al)
        if len(alb) == 0:       continue

        x = []
        y = []
        for albt in alb:
            s0, s1 = get_best_two_scores(albt)
            x.append(s0)            ;       y.append(1)
            x.append(s1)            ;       y.append(-1)

        x = N.array(x)
        y = N.array(y)

        ind = N.isfinite(x)
        x = x[ind]
        y = y[ind]

        if (y == 1).sum() == 0:         continue         # if for a particular class, we cannot decide the cut, then we keep all members of that class
        if (y == -1).sum() == 0:        continue         # if for a particular class, we cannot decide the cut, then we keep all members of that class

        x = x.reshape((-1,1))

        svc = svm.LinearSVC()
        svc.fit(x, y)

        cutoff = (- svc.intercept_) / svc.coef_
        cutoff = float(cutoff)

        assert      N.all(N.sign(x.flatten() - cutoff).astype(int) == svc.predict(x))     # verify if the calculated cutoff is correct

        cuts[c] = cutoff
   
    return cuts




'''
filter the data json according to given cutoffs
'''
def filter_dj_given_cuts(dj, cuts):
    djn = []
    for d in dj:
        if 'template' not in d:         continue
        if (d['template']['id'] in cuts) and (d['score'] < cuts[d['template']['id']]):      continue        # if for a particular class, we cannot decide the cut, then we keep all members of that class
        djn.append(d)

    return djn


def do_filter(al, dj):
    print 'pursuit.multi.recursive.filtering.second_largest_cut.do_filter()'

    cuts = find_cutoff(al)

    djf = filter_dj_given_cuts(dj, cuts)
    print 'filtered subtomogram number', len(djf)

    return djf


def load_data(fs):
    print 'pursuit.multi.recursive.filtering.second_largest_cut.do_filter()'
    fss = fs['passes'][fs['pass_i_current']]

    with open(fss['align_template_file'], 'rb') as f:       al = pickle.load(f)
    al = [_.result['align'] for _ in al]

    with open(fss['data_json_file']) as f:       dj = json.load(f)

    return {'al':al, 'dj':dj}
    


def main():
    with open('pursuit_multi_recursive_filter_second_largest_cut__op.json') as f:       op = json.load(f)

    with open(op['file_stat_file']) as f:   fs = json.load(f)
    fs['passes'] = {int(_):fs['passes'][_] for _ in fs['passes']}

    d = load_data(fs)
    djf = do_filter(d['al'], d['dj'])

    with open(op['data_json_file_out'], 'w') as f:       json.dump(djf, f, indent=2)




if __name__ == '__main__':
    main()

