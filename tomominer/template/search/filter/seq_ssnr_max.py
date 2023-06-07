#!/usr/bin/env python





'''

~/ln/tomominer/tomominer/template/search/filter/seq_ssnr_max.py

'''





import os, json, copy

import numpy as N

import tomominer.common.obj as CO
from tomominer.parallel.queue_master import QueueMaster
from tomominer.io.cache import Cache

import tomominer.statistics.ssnr as SS


def main():
    with open('template_search_filter_seq_ssnr_max__op.json') as f:       op = json.load(f)


    # just get a generic object that can contain a cache
    self = CO.Object()

    self.pool = None
  
    self.cache = Cache(tmp_dir=op['options']['tmp_dir'])

    self.runner = QueueMaster(op['options']['network']['qhost'], op['options']['network']['qport'])



    with open(op['data json in']) as f:  dj = json.load(f)


        
    dj = sorted(dj, key=lambda _:_['score'], reverse=True)

    # sequential SSNR
    n_worker = max(self.runner.work_queue.get_worker_number(), 1)
    n_chunk = max(op['options']['min_chunk_size'], int( N.ceil(len(dj) / n_worker) ))
   
    ssnr_sequential_op = {}
    ssnr_sequential_op['ssnr'] = copy.deepcopy(op['ssnr'])
    ssnr_sequential_op['n_chunk'] = n_chunk
    ssnr_s = SS.ssnr_sequential_parallel(self=self, data_json_dict={0:dj}, op=ssnr_sequential_op)
    fsc_sum = N.array([N.sum(_) for _ in ssnr_s[0]['fsc']])

    for i, d in enumerate(dj):        d['fsc'] = {'sum':float(fsc_sum[i])}

    fsc_sum_i = N.argmax(fsc_sum)

    print 'fsc_sum_i', fsc_sum_i

    with open(op['data json out'], 'w') as f:       json.dump(dj[:(fsc_sum_i+1)], f, indent=2)





if __name__ == '__main__':

    main()

