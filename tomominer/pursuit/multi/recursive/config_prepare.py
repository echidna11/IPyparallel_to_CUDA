#!/usr/bin/env python


'''
prepare files for recursive multi-pursuit

this script starts at a folder that contains pursuit_op.json, data_config.json and classify

Given pursuit result, for each cluster, create a folder. 
Then create data_config.json file (also merge ground truth from original file); link to pursuit_op.json file; and create cluster_avg_high_fsc_consistency_with_ground_truth.json file if needed


~/ln/tomominer/tomominer/pursuit/multi/recursive/config_prepare.py

'''

import os
import sys
import json
import pickle
import copy

import tomominer.io.path_util as PMAP

def export_one_set(dj_clus, dj_dict, tem_dir, cah_json, pursuit_op__file=None, remove_alignment_info=True):


    print tem_dir, len(dj_clus)

    stat = {}


    dj_clus_mix = []
    for r in dj_clus:
        r = copy.deepcopy(r)
        if 'id' in dj_dict[r['subtomogram']]:   r['id'] = dj_dict[r['subtomogram']]['id']
        if 'cluster_label' in dj_dict[r['subtomogram']]:        r['cluster_label'] = dj_dict[r['subtomogram']]['cluster_label']        # replace predicted cluster label with ground truth for analysis purpose!!

        if remove_alignment_info:
            # remove alignment information to avoid the new pursuit stuck at local minima
            if 'angle' in r:    del r['angle']
            if 'loc' in r:  del r['loc']
            if 'score' in r: del r['score']
            if 'template' in r: del r['template']

        dj_clus_mix.append(r)

    if not os.path.isdir(tem_dir):   os.makedirs(tem_dir)
    stat['tem_dir'] = tem_dir


    tem__data_config_file = os.path.join(tem_dir, 'data.json')
    stat['tem__data_config_file'] = tem__data_config_file
    if not os.path.isfile(tem__data_config_file):
        with open(tem__data_config_file, 'w') as f:        json.dump(dj_clus_mix, f, indent=2)

    stat['size'] = len(dj_clus_mix)


    if cah_json is not None:
        cah_json_t = copy.deepcopy(cah_json)
        cah_json_t['true_avg'] = os.path.abspath(cah_json_t['true_avg'])        ;       assert os.path.isfile(cah_json_t['true_avg'])
        cah_json_t['true_member'] = tem__data_config_file           # this file actually may not contain the ground truth....
        cluster_avg_high_fsc_consistency_with_ground_truth__op__file = os.path.join(tem_dir, 'cluster_avg_high_fsc_consistency_with_ground_truth__op.json')
        stat['cluster_avg_high_fsc_consistency_with_ground_truth__op__file'] = cluster_avg_high_fsc_consistency_with_ground_truth__op__file
        with open(cluster_avg_high_fsc_consistency_with_ground_truth__op__file, 'w') as f:     json.dump(cah_json_t, f, indent=2)

    if pursuit_op__file is not None:
        tem__pursuit_op__file = os.path.join(tem_dir, 'pursuit_op.json')
        stat['tem__pursuit_op__file'] = tem__pursuit_op__file
        if not os.path.exists(tem__pursuit_op__file):           os.symlink(pursuit_op__file, tem__pursuit_op__file)

    return stat



def main():

    stat = {}
    stat['op_file'] = os.path.abspath('pursuit_multi_recursive_config_prepare__op.json')

    with open(stat['op_file']) as f:   op = json.load(f)


    if 'pursuit op file' in op:
        op['pursuit op file'] = os.path.abspath(op['pursuit op file'])
        assert  os.path.isfile(op['pursuit op file'])
    else:
        op['pursuit op file'] = None

    with open(op['file_stat file']) as f:   fs = json.load(f)
    fs['passes'] = {int(_):fs['passes'][_] for _ in fs['passes']}

    # to configure after particular pass or after last pass
    if 'pass_i' in op:
        pass_i = op['pass_i']
    else:
        pass_i = fs['pass_i_current']

    print 'use pass', pass_i

    p_sel = fs['passes'][pass_i]

    print 'loading', p_sel['data_json_file']

    with open(p_sel['data_json_file']) as f:     dj = json.load(f)
    dj_dict = {_['subtomogram']:_ for _ in dj}


    if ('cluster avg high fsc consistency with ground truth config file' in op) and os.path.isfile(op['cluster avg high fsc consistency with ground truth config file']):
        with open(op['cluster avg high fsc consistency with ground truth config file']) as f:       cah_json = json.load(f)
    else:
        cah_json = None

    print 'loading', p_sel['cluster_average_select_file']
    with open(p_sel['cluster_average_select_file']) as f:       cas = pickle.load(f)

    stat['clusters'] = {}

    # create folder for selected template key, then export files
    for c in cas['selected_templates']:
        tk = cas['selected_templates'][c]['subtomogram']
        dj_clus = cas['tk_info'][tk]['data_json']
        
        tem_dir = os.path.abspath(os.path.join(op['out dir'], 'clus-%03d'%(c,)))
        stat['clusters'][c] = export_one_set(dj_clus=dj_clus, dj_dict=dj_dict, tem_dir=tem_dir, cah_json=cah_json, pursuit_op__file=op['pursuit op file'], remove_alignment_info=op['remove_alignment_info'])


    # also create a folder for rest of subtomograms
    if True:
        selected_subtomograms = set()
        for c in cas['selected_templates']:
            tk = cas['selected_templates'][c]['subtomogram']
            dj_clus = cas['tk_info'][tk]['data_json']
            selected_subtomograms.update(   [_['subtomogram'] for _ in dj_clus]  )

        print 'data.json num', len(dj), '\t\t', 'selected_subtomograms num', len(selected_subtomograms)

        tem_dir = os.path.abspath(os.path.join(op['out dir'], 'clus-rest'))
        stat['rest'] = export_one_set(dj_clus=[_ for _ in dj if (_['subtomogram'] not in selected_subtomograms)], dj_dict=dj_dict, tem_dir=tem_dir, cah_json=cah_json, pursuit_op__file=op['pursuit op file'], remove_alignment_info=op['remove_alignment_info'])

    with open(op['stat file out'], 'w') as f:       json.dump(stat, f, indent=2)

    print 'WARNING: BE AWARE that the exported clusters are selected clusters, which could be from any past iterations!!  the exported alignment information is the alignment of last iteration!!'


if __name__ == '__main__':
    main()



