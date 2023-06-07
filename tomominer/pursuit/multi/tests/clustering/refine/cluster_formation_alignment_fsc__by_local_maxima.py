'''
test
pursuit.multi.util.cluster_formation_alignment_fsc__by_local_maxima()
'''


import os, sys, json
import cPickle as pickle

import numpy as N


from multiprocessing.pool import Pool
from tomominer.parallel.queue_master import QueueMaster

import tomominer.common.obj as CO
from tomominer.io.cache import Cache



data_dir = '/home/rcf-47/mxu/ln/frequent_structure/data/out/method/imputation-classification/multi/classification/30/0.005/0170/classify/pass_005'
with open(os.path.join(data_dir, 'data_config.json')) as f:     dj = json.load(f)

with open(os.path.join(data_dir, 'cluster_info_list.json')) as f:                    cluster_info_list = json.load(f)
cluster_info_list = {int(_):cluster_info_list[_] for _ in cluster_info_list}

cluster_info = {}
for pass_i_t, cluster_info_file_t in cluster_info_list.iteritems():
    with open(cluster_info_file_t, 'rb') as f:                    cluster_info[pass_i_t] = pickle.load(f)


# just get a generic object that can contain a cache
self = CO.Object()

self.pool = Pool()
self.cache = Cache(tmp_dir='/home/rcf-47/mxu/tmp-panasas2')

self.runner = QueueMaster('localhost', 5011)


'''
parameters according to
/home/rcf-47/mxu/ln/frequent_structure/data/out/analysis/jensen/tocheva/130401-complex/classification/classification/0080/classification/0055/classify_config.json
since it gets good results
'''

from tomominer.pursuit.multi.util import tomominer.cluster_formation_alignment_fsc__by_local_maxima
b = cluster_formation_alignment_fsc__by_local_maxima(self=self, dj=dj, ci=None, op={'min_expansion_size':50, 'max_expansion_size':2000, 'n_chunk':50, 'gaussian_smooth_sigma':10, 'debug':True, 'averaging':{'use_fft':True, 'with_missing_wedge':True, 'mask_count_threshold':50, 'mask_binarize':False, 'centerize_loc':True}, 'ssnr_sequential':{'segmentation':{'smooth_weight':1.0, 'image_weight':1.0, 'phi_propotion_cutoff':-2.0, 'density_positive':False}}})


import pdb          ;           pdb.set_trace()


