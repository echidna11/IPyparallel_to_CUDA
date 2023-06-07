#!/usr/bin/env python



'''

across multiple experiments including multiple pursuit and GA averaging, 
given the corresponding file stat file,
collect FSC, FSC sum score, FSC resolution (with 0.5 cutoff). 
Can export a tab deliminated table so that we can order them using excel

maybe can also export averages for visual inspection


question: shall we deal with subtomogram membership overlap?????? In that case we can still use the greedy approach....



~/ln/tomominer/tomominer/pursuit/multi/analysis/statistics/collect_fsc_across_multiple_experiments.py

'''



import json, csv
import cPickle as pickle

import numpy as N


import tomominer.statistics.vol as SV


def collect_averaging(fs, pass_i=None):
    ps = fs['passes']
    ps = {int(_):ps[_] for _ in ps}

    if pass_i is not None:
        psb = ps[pass_i]
    else:
        psb = ps[fs['best_pass_i']]

    with open(psb['pmpg_file'], 'rb') as f:     p = pickle.load(f)

    pb = p['best']['e'][0]

    r = {}          # return info
    r[0] = {}
    r[0]['fsc'] = pb['fsc'].tolist()
    r[0]['dj'] = p['dj']

    with open(psb['cluster_averaging_file'], 'rb') as f:    ca = pickle.load(f)

    r[0]['template'] = ca['template_keys'][0]

    r[0]['mode'] = 'averaging'


    return r




def collect_pursuit(fs, pass_i=None):
    ps = fs['passes']
    ps = {int(_):ps[_] for _ in ps}

    ci = {}         # collect cluster info, because it contains FSC details
    for pass_i in ps:        
        with open(ps[pass_i]['cluster_info_file']) as f:       ci[pass_i] = pickle.load(f)

    if pass_i is not None:
        psc = ps[pass_i]
    else:
        psc = ps[fs['pass_i_current']]

    with open(psc['fsc_stat_file']) as f:       fscs = json.load(f)
    fscs = {int(_):fscs[_] for _ in fscs}
    fscs = fscs[fs['pass_i_current']]
    fscs = {int(_):fscs[_] for _ in fscs}


    r = {}      # return info

    for c in fscs:
        cit = ci[fscs[c]['pass_i']][fscs[c]['cluster']]
        assert cit['pass_i'] == fscs[c]['pass_i']
        assert cit['cluster'] == fscs[c]['cluster']

        r[c] = {'fsc': cit['fsc'].tolist(), 'dj':cit['data_json'], 'mode': 'pursuit'}


    with open(psc['cluster_average_select_file'], 'rb') as f:       st = pickle.load(f)['selected_templates_common_frame']
    
    for c in st:        r[c]['template'] = st[c]

    return r


def resolution(co, voxel_spacing, cutoff):
    
    for fs in co:
        for c in co[fs]:
            t = co[fs][c]
            t['resolution'] = SV.fsc_resolution(N.array(t['fsc']), voxel_spacing=voxel_spacing, cutoff=cutoff)


def export_csv(co, fn):
    import csv
    with open(fn, 'wb') as f:
        w = csv.writer(f)

        for fs_fn in co:
            for c in co[fs_fn]:
                t = co[fs_fn][c]
                w.writerow([N.sum(t['fsc']), t['resolution'], t['template']['subtomogram']])



def main():

    with open('collect_fsc_across_multiple_experiments__op.json') as f:         op = json.load(f)

    with open(op['file stats']) as f:         fs_s = json.load(f)

    co = {}
    for fs_fn in fs_s:
        print fs_fn
        with open(fs_fn['file']) as f:      fs = json.load(f)

        if 'best_pass_i' in fs:
            # in this case, fs is for GA averaging
            co[fs_fn['file']] = collect_averaging(fs, pass_i=(fs_fn['pass_i'] if 'pass_i' in fs_fn else None))
        else:
            # in this case, fs if for pursuit
            co[fs_fn['file']] = collect_pursuit(fs, pass_i=(fs_fn['pass_i'] if 'pass_i' in fs_fn else None))

    resolution(co, voxel_spacing=op['voxel spacing'], cutoff=op['resolution cutoff'])
    
    with open(op['collect out'], 'wb') as f:     pickle.dump(co, f, protocol=-1)

    export_csv(co, fn=op['summary csv out'])


if __name__ == '__main__':
    main()

