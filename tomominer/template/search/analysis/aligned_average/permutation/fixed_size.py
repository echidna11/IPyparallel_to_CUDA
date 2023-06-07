#!/usr/bin/env python




'''

given a list of subtomograms aligned to a template, randomly sample one set with a fixed size, 
calculate average, and the FSC score of this average to the template


~/ln/tomominer/tomominer/template/search/analysis/aligned_average/permutation/fixed_size.py

'''


import os, shutil, copy, json, random, uuid
import numpy as N

import tomominer.io.file as IF

from tomominer.parallel.queue_master import QueueMaster
import tomominer.common.obj as CO
from tomominer.io.cache import Cache

import tomominer.pursuit.multi.util as PMU
import tomominer.statistics.vol as SV


def main():

    with open('fixed_size__op.json') as f:      op = json.load(f)

    op['options']['tmp_dir'] = os.path.join(os.path.abspath(op['options']['tmp_dir']), str(uuid.uuid4()))
    print 'tmp dir', op['options']['tmp_dir'] 

    if os.path.isdir(op['options']['tmp_dir']):
        shutil.rmtree(op['options']['tmp_dir'])
        os.makedirs(op['options']['tmp_dir'])


    self = CO.Object()
    self.cache = Cache(tmp_dir=op['options']['tmp_dir'])
    self.runner = QueueMaster(op['options']['network']['qhost'], op['options']['network']['qport'])

    with open(op['data']) as f:     dj = json.load(f)

    with open(op['reference_set']) as f:    s_ref = json.load(f)
    s_ref = set([_['subtomogram'] for _ in s_ref])
    print 'reference set size', len(s_ref)


    avg_op = copy.deepcopy(op['averaging'])
    avg_op['mask_count_threshold'] = op['min_sample_num']
    avg_op['out_dir'] = os.path.join(op['options']['tmp_dir'], 'reference')
    avg_op['n_chunk'] = op['options']['min_chunk_size']

    t = IF.read_mrc_vol(op['template'])

    fsc_sum_ref = average_fsc_sum(self=self, djs={0:[_ for _ in dj if (_['subtomogram'] in s_ref)]}, t=t, avg_op=avg_op)
    
    djt = {}
    for i in range(op['repeat_num']):
        djt[i] = random.sample(dj, len(s_ref))



    n_worker = max(self.runner.work_queue.get_worker_number(), 1)
    n_chunk = max(op['options']['min_chunk_size'], int( N.ceil(len(s_ref) * op['repeat_num'] / n_worker) ))
    avg_op['out_dir'] = os.path.join(op['options']['tmp_dir'], 'sample')
    avg_op['n_chunk'] = n_chunk

    fsc_sum_permutation = average_fsc_sum(self=self, djs=djt, t=t, avg_op=avg_op)

    with open(op['stat_file'], 'w') as f:       json.dump({'fsc_sum':{'ref':fsc_sum_ref, 'permutation':fsc_sum_permutation}}, f, indent=2)

    #shutil.rmtree(op['options']['tmp_dir'])


def average_fsc_sum(self, djs, t, avg_op):

    # averaging
    avg_re = PMU.cluster_averaging(self=self, clusters=djs, op=avg_op)

    fsc_sum = {}
    for i in djs:
        v = IF.read_mrc_vol(avg_re['template_keys'][i]['subtomogram'])
        fsc_sum[i] = SV.fsc(t, v).sum()

    return fsc_sum

if __name__ == '__main__':
    main()




'''

# code for analysis

import json

with open('fixed_size__stat_out.json') as f:    s = json.load(f)

s_ref = s['fsc_sum']['ref']['0']
s_perm = s['fsc_sum']['permutation']

import numpy as N
s_perm = N.array([s_perm[_] for _ in s_perm])

print N.sum(s_perm > s_ref) / float(len(s_perm))


# assume  normal distribution
import scipy.stats as SS
1 - SS.norm.cdf(s_ref, loc=s_perm.mean(), scale=s_perm.std())





# export for R, for normality check
import csv
with open('/tmp/t.csv', 'w') as f:
    w = csv.writer(f)
    w.writerows(N.reshape(s_perm, (-1,1)).tolist())





# R commands for analysis

v = read.csv('/tmp/t.csv')
v = as.matrix(v)[,1]

hist(v, 100)
qqnorm(v)


# regarding to test of normality, see              http://stackoverflow.com/questions/7781798/seeing-if-data-is-normally-distributed-in-r


'''


