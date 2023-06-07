
'''

after performing hierarchical clustering, first calculate knee to (hopfully) obtain a optimal curoff, 
then under such cutoff find the non-overlaping clusters with highest FSC scores

test such idea by examing pursuit process of early aligned and late aligned passes

'''

import os, sys
import cPickle as pickle
import numpy as N
import copy

data_dir = '/home/rcf-47/mxu/ln/frequent_structure/data/out/method/imputation-classification/multi/classification/30/0.005/0130/classify/pass_029'         # we can expect a large number of clusters in such case. Result: the size distribution of clusters is more even. There is no single dominated cluster. However, the SCS.hierarchy_ssnr_cluster_select() selected clusters with almost similiar size as the ones cutted by knee from the hierarchy. So this means SSNR based selection does not further discriminate mixed clusters. Wuwuwuw. Also, if we do not use knee, say do_hi_filter = False, then SCS.hierarchy_ssnr_cluster_select() selected a huge cluster that contains all subtomograms
#data_dir = '/home/rcf-47/mxu/ln/frequent_structure/data/out/method/imputation-classification/multi/classification/30/0.005/0130/classify/pass_000'         # we can expect a smaller number of clusters in such case, only seperated by size of particles. Result: only one domainate cluster
#data_dir = '/home/rcf-47/mxu/ln/frequent_structure/data/out/analysis/jensen/tocheva/130401-complex/classification/classification/0080/classification/0055-rep/out/001/pass_000'     # we can expect one domainated cluster in such case
#data_dir = '/home/rcf-47/mxu/ln/frequent_structure/data/out/analysis/jensen/tocheva/130401-complex/classification/classification/0080/classification/0055-rep/out/001/pass_006'     # we can expect one domainated number of clusters in such case


print 'load dimension reduction'
with open(os.path.join(data_dir, 'dimension_reduction.pickle'), 'rb') as f:     dr = pickle.load(f)

print 'hierarchical clustering'
import tomominer.pursuit.multi.hierarchy_ssnr_calculate_select as SCS
hc_re = SCS.hierarchical_clustering(x=dr['cfp_re']['red'])

import pickle
with open('/tmp/ws.p', 'wb') as f:       pickle.dump({'dr':dr, 'hc_re':hc_re}, f, protocol=-1)

'''
import pickle
with open('/tmp/ws.p', 'rb') as f:       tmp = pickle.load(f)
dr = tmp['dr']
hc_re = tmp['hc_re']
'''

print 'calculate knee'
knee = {}
import tomominer.cluster.hierarchy.hierarchy as CHH
knee['chh_cnvdc'] = CHH.cluster_number_vs_dist_cutoff(r=hc_re['root_id'], hi=hc_re['info'])

import tomominer.cluster.hierarchy.knee_l as CHK
knee['chk_cnvdc'] = CHK.cluster_number_vs_dist_cutoff(knee['chh_cnvdc'])
knee['chk_r'] = CHK.refine(nc=knee['chk_cnvdc']['nc'], em=knee['chk_cnvdc']['md'])


# inspect cluster size distribution
knee['chh_cc'] = CHH.cluster_cut(r=hc_re['root_id'], hi=hc_re['info'], dist_cutoff=knee['chk_cnvdc']['md'][knee['chk_r']['ci']])
print [hc_re['info'][_]['size'] for _ in knee['chh_cc']]

knee['dist_cutoff'] = knee['chk_cnvdc']['md'][knee['chk_r']['ci']]



# calculate clusters of highest FSC scores

import json
import warnings
from multiprocessing.pool import Pool
from tomominer.parallel.queue_master import QueueMaster

import tomominer.common.obj as CO
from tomominer.io.cache import Cache

self = CO.Object()
self.pool = Pool()

tmp_dir = '/home/rcf-47/mxu/tmp-panasas2'
self.cache = Cache(tmp_dir=tmp_dir)

self.runner = QueueMaster('localhost', 5011)

with open(os.path.join(data_dir, '../../classify_config.json')) as f:       op = json.load(f)
with open(os.path.join(data_dir, 'data_config.json')) as f:       data_json = json.load(f)

clusters_populated = {}

op['cluster']['ssnr']['segmentation'] = True
op['segmentation'] = op['segmentation bak']
op['template']['guided_segmentation'] = op['template']['guided_segmentation bak']
if op['cluster']['ssnr']['segmentation']:
    segmentation_op = copy.deepcopy(op['segmentation'])
    segmentation_op['density_positive'] = op['density_positive']
    segmentation_tg_op = segmentation_op
    if ('guided_segmentation' in op['template']) and ('gaussian_smooth_sigma' in op['template']['guided_segmentation']):            segmentation_tg_op['gaussian_smooth_sigma'] =  op['template']['guided_segmentation']['gaussian_smooth_sigma']

else:
    segmentation_tg_op = None

SCS.hierarchy_ssnr_calculate(  self=self, hc=hc_re, dj=data_json, cluster_size_min=op['cluster']['size_min'], segmentation_op=segmentation_tg_op  )


do_hi_filter = False
if do_hi_filter:
    hi_filtered = {_:hc_re['info'][_] for _ in hc_re['info'] if (hc_re['info'][_]['dist'] <= knee['dist_cutoff'])}
else:
    hi_filtered = hc_re['info']

scs_csl_re = SCS.cluster_select_label(hi=hc_re['info'], ids=SCS.hierarchy_ssnr_cluster_select(hi=hi_filtered), dj_len=len(data_json), starting_label=(N.max([_ for _ in clusters_populated])+1 if (len(clusters_populated) > 0) else 0))
 

l = N.array(scs_csl_re['labels'])
print [(l == _).sum() for _ in scs_csl_re['label_to_id']]




