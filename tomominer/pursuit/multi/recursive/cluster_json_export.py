#!/usr/bin/env python

# for each selected cluster, create a subfolder and export the data_json field

import os
import json
import copy
import cPickle as pickle


if __name__ == '__main__':
    raise       Exception('deprecated, use config_prepare.py instead')

    with open('data_config.json') as f:     dj = json.load(f)

    with open('cluster_average_select.pickle', 'rb') as f:   cs = pickle.load(f)

    st = cs['selected_templates']
    ti = cs['tk_info']

    tif = {_:ti[st[_]['subtomogram']]   for _ in st}            # keep only selected cluster's information, indixed by cluster id

    del cs, st, ti


    for c in tif:
        
        out_dir = os.path.join(os.getcwd(), 'recursive', 'clus-%d'%(c,))
        if not os.path.isdir(out_dir):      os.makedirs(out_dir)
        
        dj_t = copy.deepcopy(tif[c]['data_json'])


        # remove previous alignment information so that the new pursuit won't go to local minimima
        for d in dj_t:          del d['angle'], d['loc'], d['score']

        with open(os.path.join(out_dir, 'data_config.json'), 'w') as f:     json.dump(djt, f, indent=2)

    
    if True:
        # generate a list of subtomograms uncovered by the clusters
        clus_cover_set = set()
        for c in tif:
            dj_t = tif[c]['data_json']
            clus_cover_set.update(  [_['subtomogram'] for _ in dj_t]  )


        dj_rest = []
        for d in dj:
            if d['subtomogram'] in clus_cover_set:      continue
            del d['angle'], d['loc'], d['score']

            dj_rest.append(d)

        out_dir = os.path.join(os.getcwd(), 'recursive', 'clus-rest')
        if not os.path.isdir(out_dir):      os.makedirs(out_dir)
        with open(os.path.join(out_dir, 'data_config.json'), 'w') as f:     json.dump(dj_rest, f, indent=2)

