#!/usr/bin/env python



# compare the pursuit result with template search result, see if the high score peaks of both cases tend to stay very close


if __name__ == '__main__':

    import json
    with open('template_colocalization__config.json') as f:     op = json.load(f)

    with open(op['template_peak_alignment_file']) as f:     tp = json.load(f)       # this file contains both peak location and alignment scores from tomominer.template search
    tp = sorted(tp, key=lambda _:(-_['score']))         # order according to scores
    if 'top_num' in op:     tp = tp[:op['top_num']]

    import cPickle as pickle
    with open(op['cluster_average_select_file'], 'rb') as f:    cas_re = pickle.load(f)
    st = cas_re['selected_templates']
    tk_info = {_:cas_re['tk_info'][st[_]['subtomogram']] for _ in st}
    del cas_re

    with open(op['pursuit_peak_file']) as f:     pp = json.load(f)     # this file contains peak locations
    pp = {_['subtomogram']:_ for _ in pp}


    import numpy as N
    tx = N.array([_['peak']['loc'] for _ in tp])      # location of peaks in template matching

    import os, sys
    import numpy as N
    import time

    import tomominer.geometry.point_cloud.colocalization.util as GPCU
    for avg_id in tk_info:
        print 'avg_id', avg_id,

        start_time = time.time()
        dj = tk_info[avg_id]['data_json']
        dj = sorted(dj, key=lambda _:(-_['score']))         # order according to scores
        if 'top_num' in op:     dj = dj[:op['top_num']]

        sx = N.array([pp[_['subtomogram']]['peak']['loc'] for _ in dj])      # location of peaks in pursuit

        print 'point number for matching', len(tx), len(sx),            ;   sys.stdout.flush()

        if op['exclusive']:
            m = GPCU.point_pair_match_exclusive(tx, sx)
        else:
            m = GPCU.point_pair_match_nonexclusive(tx, sx)

        tab = []
        for mi in range(len(m['i'])):
            d = m['d'][mi]
            i = m['i'][mi]

            tab.append([tp[i[0]]['score'], dj[i[1]]['score'], d])


        out_file = 'colocalization-%d.txt'%(avg_id,)
        N.savetxt(out_file, tab, fmt='%f', delimiter=' ')

        print 'time used', time.time() - start_time





'''

# plotting in matlab


clear all;
close all;
t = load('colocalization-1.txt');

figure;     plot3(t(:,1), t(:,2), t(:,3), 'k.');        xlabel('align score to template');       ylabel('align score to average');      zlabel('peak distance');

i = t(:,3)==0;
figure;     plot(t(i,1), t(i,2), 'k.');         xlabel('align score to template');       ylabel('align score to average');



'''

