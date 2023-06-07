#!/usr/bin/env python


'''
simple average of a set of already aligned subtomograms. 

Modified from  
~/ln/tomominer/tomominer/average/genetic_algorithm/main.py



~/ln/tomominer/tomominer/average/aligned/simple_average.py

'''


import os, json, copy, random

import numpy as N


from tomominer.parallel.queue_master import QueueMaster
import tomominer.common.obj as CO
from tomominer.io.cache import Cache


import tomominer.pursuit.multi.util as PMU


def do_average(self, op, data_json):

    out_dir = os.path.abspath(op['out_dir'])
    if not os.path.isdir(out_dir):      os.makedirs(out_dir)


    n_worker = max(self.runner.work_queue.get_worker_number(), 1) * 2
    n_chunk = max(op['options']['min_chunk_size'], int( N.ceil(len(data_json) / n_worker) ))

    #-----------------------------------------------------------
    # calculate SSNR based FSC

    if 'smooth' in op['averaging']:
        cluster_ssnr_fsc__op = {}
        cluster_ssnr_fsc__op['ssnr'] = copy.deepcopy(op['ssnr'])
        cluster_ssnr_fsc = PMU.cluster_ssnr_fsc(self=self, clusters={0:data_json}, n_chunk=n_chunk, op=cluster_ssnr_fsc__op)
        print 'fsc_sum', N.sum(cluster_ssnr_fsc['fsc'][0])
 
    #--------------------------------------------------------------------
    # calculate global average
    avg_op = copy.deepcopy(op['averaging'])
    avg_op['mask_count_threshold'] = op['min_sample_num']
    avg_op['n_chunk'] = n_chunk
    avg_op['out_dir'] = op['out_dir']
    if 'smooth' in avg_op:                avg_op['smooth']['fsc'] = {0:cluster_ssnr_fsc['fsc'][0]}

    avg_re = PMU.cluster_averaging(self=self, clusters={0:data_json}, op=avg_op)



def main():

    with open('average__op.json') as f:     op = json.load(f)

    # just get a generic object that can contain a cache
    self = CO.Object()

    self.cache = Cache(tmp_dir=op['options']['tmp_dir'])
    self.runner = QueueMaster(op['options']['network']['qhost'], op['options']['network']['qport'])
    self.pool = None

    with open(op['data_file']) as f:    dj = json.load(f)

    if 'random_sample_num' in op:
        # mainly for testing
        dj = random.sample(dj, op['random_sample_num'])

    print 'averaging', len(dj), 'subtomograms'
    do_average(self=self, op=op, data_json=dj)



if __name__ == '__main__':
    main()

