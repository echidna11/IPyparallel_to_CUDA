
'''
first perform pursuit based on kmeans mode, then cut the hierarchical clusters so that the result grouping is most consistent with hierarchical clustering.

see
https://docs.google.com/document/d/1xU1AWj--Uol-2fx_xhgToneeRjmOEtgAOW04_PF_uCU/edit

'''



import os, sys, json
import cPickle as pickle
import numpy as N
import copy



pass_dir = '/home/rcf-47/mxu/ln/frequent_structure/data/out/method/imputation-classification/multi/classification/30/0.005/0150/classify/pass_016'
pass_dir_last = '/home/rcf-47/mxu/ln/frequent_structure/data/out/method/imputation-classification/multi/classification/30/0.005/0150/classify/pass_015'


#--------------------------------------------------
print 'get selected non-specific clusters of last pass'
with open(os.path.join(pass_dir_last, 'cluster_average_select.pickle'), 'rb') as f:     cas = pickle.load(f)

with open(os.path.join(pass_dir_last, 'data_config.json')) as f:     dj = json.load(f)

with open(os.path.join(pass_dir_last, 'cluster_info_list.json')) as f:     cil = json.load(f)

ci = {}
for pass_i in cil:        
    with open(cil[pass_i], 'rb') as f:     ci[int(pass_i)] = pickle.load(f)

with open(os.path.join(pass_dir_last, 'cluster_info_stat.pickle'), 'rb') as f:     cis = pickle.load(f)


dj_sel = {}
for i in cas['selected_templates']:
    k = cas['selected_templates'][i]
    assert      ci[k['pass_i']][k['cluster']]['template_key']['subtomogram'] == k['subtomogram']
    if cis[k['pass_i']][k['cluster']]['is_specific'] is not None:       continue
    dj_sel[i] = cas['tk_info'][k['subtomogram']]['data_json']

print sorted([len(dj_sel[_]) for _ in dj_sel], reverse=True)      # size distribution



dj_sel_dict = {}            # the dictionary of subtomogram's set id, notice that it may have a little ambiguity when we allow a little bit of overlap between sets in the selection process
for i in dj_sel:
    for _ in dj_sel[i]:
        dj_sel_dict[_['subtomogram']] = i

label_sel = [N.nan] * len(dj)
for i, _ in enumerate(dj):
    if _['subtomogram'] not in dj_sel_dict:     continue
    label_sel[i] = dj_sel_dict[_['subtomogram']]



#--------------------------------------------------
print 'load dimension reduction'
with open(os.path.join(pass_dir, 'dimension_reduction.pickle'), 'rb') as f:     dr = pickle.load(f)
assert  len(dj) == len(dr['cfp_re']['red'])

print 'hierarchical clustering'
import tomominer.pursuit.multi.hierarchy_ssnr_calculate_select as SCS
hc_re = SCS.hierarchical_clustering(x=dr['cfp_re']['red'])


'''
import pickle
with open('/tmp/ws.p', 'wb') as f:       pickle.dump({'dr':dr, 'hc_re':hc_re}, f, protocol=-1)
'''


#--------------------------------------------------
print 'select cutoff for most consistent clustering'

import tomominer.cluster.hierarchy.hierarchy as CHH
best_cut = CHH.optimal_consistent_cut(z=hc_re['clusters'], l0=label_sel, find_min_dist=True, consistency_matric='nmi')


# inspect cluster size distribution
chh_cc = CHH.cluster_cut(r=hc_re['root_id'], hi=hc_re['info'], dist_cutoff=best_cut['t'])
print sorted([hc_re['info'][_]['size'] for _ in chh_cc], reverse=True)




'''
# when file_stat.json is avaliable, can use it to load data instead

file_stat_file = '/home/rcf-47/mxu/ln/frequent_structure/data/out/method/imputation-classification/multi/classification/30/0.005/0180/out/file_stat.json'
with open(file_stat_file) as f:     fs = json.load(f)
fs['passes'] = {int(_):fs['passes'][_] for _ in fs['passes']}

# load the and get class labels of selected (and non-specific) sets of last iteration
with open(fs['passes'][fs['pass_i_current']-1]['cluster_average_select_file'], 'rb') as f:      cas = pickle.load(f)
with open(fs['passes'][fs['pass_i_current']-1]['cluster_info_file'], 'rb') as f:      ci = pickle.load(f)
with open(fs['passes'][fs['pass_i_current']-1]['cluster_info_stat_file'], 'rb') as f:      cist = pickle.load(f)
with open(fs['passes'][fs['pass_i_current']-1]['data_json_file']) as f:      dj = json.load(f)
'''


