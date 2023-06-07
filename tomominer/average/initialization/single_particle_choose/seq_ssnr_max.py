#!/usr/bin/env python



'''
give a collection of subtomograms, select one, then align all subtomograms against this one. Then order them according to alignment score, then calculate sequential SSNR. Then calculate max FSC score.



~/ln/tomominer/tomominer/average/initialization/single_particle_choose/seq_ssnr_max.py

'''


import os, json, copy

import numpy as N

import tomominer.common.obj as CO
from tomominer.parallel.queue_master import QueueMaster
from tomominer.io.cache import Cache

import tomominer.statistics.ssnr as SS


def main():
    with open('average_initialization_single_particle_choose_seq_ssnr_max__op.json') as f:       op = json.load(f)


    # just get a generic object that can contain a cache
    self = CO.Object()

    self.pool = None
  
    self.cache = Cache(tmp_dir=op['options']['tmp_dir'])

    self.runner = QueueMaster(op['options']['network']['qhost'], op['options']['network']['qport'])



    with open(op['data json in']) as f:  dj = json.load(f)

    if os.path.isfile(op['stat file out']):
        print 'loading', op['stat file out']
        with open(op['stat file out']) as f:   stat = json.load(f)
        best = stat['best']
    else:
        stat = None
        best = None

    for d_i, d in enumerate(dj):
        if (stat is not None) and (d_i <= stat['current_i']):   continue
        
        dj_t = [_ for _ in dj if (_['subtomogram'] != d['subtomogram'])]            ;           assert  len(dj) == (len(dj_t) + 1)
        
        # perform alignment
        task_priority = 2000 + N.random.randint(100)        # add a random number so that this batch of tasks stay together
        tasks = []
        for d_t in dj_t:   tasks.append(self.runner.task( priority=task_priority, module='tomominer.pursuit.multi.util', method='align_to_templates', kwargs={'rec':d_t, 'tem_keys':{0:d}, 'align_op':op['align']} ))
        at = [_ for _ in self.runner.run__except(tasks)]
        at = [_.result for _ in at]

        dj_ta = [ {'subtomogram':_['vol_key']['subtomogram'], 'mask':_['vol_key']['mask'], 'score':_['align'][0]['score'], 'loc':_['align'][0]['loc'].tolist(), 'angle':_['align'][0]['angle'].tolist()}  for _ in at  ]
        dj_ta = sorted(dj_ta, key=lambda _:_['score'], reverse=True)

        # sequential SSNR
        n_worker = max(self.runner.work_queue.get_worker_number(), 1)
        n_chunk = max(op['options']['min_chunk_size'], int( N.ceil(len(dj) / n_worker) ))
       
        ssnr_sequential_op = {}
        ssnr_sequential_op['ssnr'] = copy.deepcopy(op['ssnr'])
        ssnr_sequential_op['n_chunk'] = n_chunk
        ssnr_s = SS.ssnr_sequential_parallel(self=self, data_json_dict={0:dj_ta}, op=ssnr_sequential_op)
        fsc_sum = N.array([N.sum(_) for _ in ssnr_s[0]['fsc']])
        fsc_sum_i = N.argmax(fsc_sum)

        print '\r', '%05d %0.3f '%(d_i, float(d_i)/len(dj)),

        if (best == None) or (best['sequential_ssnr']['max_fsc_sum'] < fsc_sum[fsc_sum_i]):
            best = d
            best['i'] = d_i
            best['sequential_ssnr'] = {'max_fsc_sum':fsc_sum[fsc_sum_i], 'i':fsc_sum_i}
            print '\r', '%05d %0.3f '%(d_i, float(d_i)/len(dj)), '       ', 'fsc_sum_i %5d'%(fsc_sum_i,), '    ', 'max_fsc_sum %5.3f'%(fsc_sum[fsc_sum_i],)

            d_t = copy.deepcopy(d)
            d_t['loc'] = [0.0, 0.0, 0.0]
            d_t['angle'] = [0.0, 0.0, 0.0]
            d_t['score'] = 1.0

            dj_ta_t = [d_t]
            dj_ta_t.extend(dj_ta)       ;       assert  len(set(_['subtomogram'] for _ in dj_ta_t)) == len(dj)

            with open(op['data json out'], 'w') as f:       json.dump(dj_ta_t, f, indent=2)

        with open(op['stat file out'], 'w') as f:   json.dump({'current_i':d_i, 'best': best}, f, indent=2)




if __name__ == '__main__':

    main()

