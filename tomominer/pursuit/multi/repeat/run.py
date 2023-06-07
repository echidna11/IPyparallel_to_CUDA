#!/usr/bin/env python

'''
load classify_config and data_config file, 
then for each repeat, create out dir, and run the repeat of pursuit process
collect and record final results for all repeats
'''

import os, sys, json, copy
import warnings
import cPickle as pickle
from collections import defaultdict
import numpy as N
from multiprocessing.pool import Pool
from tomominer.parallel.queue_master import QueueMaster

import tomominer.common.obj as CO
from tomominer.io.cache import Cache
from tomominer.pursuit.multi.main import tomominer.pursuit


def main():
    warnings.filterwarnings('error')

    with open('pursuit_multi_repeat__op.json') as f:    op = json.load(f)
    op['out_dir'] = os.path.abspath(op['out_dir'])
    
    with open(op['pursuit_op_file']) as f:       pursuit_op = json.load(f)

    # just get a generic object that can contain a cache
    self = CO.Object()

    self.pool = Pool()
  
    self.cache = Cache(tmp_dir=pursuit_op['options']['tmp_dir'])

    self.runner = QueueMaster(pursuit_op['options']['network']['qhost'], pursuit_op['options']['network']['qport'])

   
    with open(op['data_json_file']) as f:      data_json = json.load(f)

    stat = {}
    for repeat_i in range(op['repeat_num']):
        print '-----------------------------------------------------------------------'
        print 'repeat', repeat_i
        print '-----------------------------------------------------------------------'
        out_dir_t = os.path.join(op['out_dir'], '%03d'%(repeat_i,))
        if not os.path.isdir(out_dir_t):        os.makedirs(out_dir_t)

        pursuit_op_t = copy.deepcopy(pursuit_op)
        pursuit_op_t['out_dir'] = out_dir_t

        file_stat = pursuit(self=self, op=pursuit_op_t, data_json=data_json)
        file_stat_t = file_stat['passes'][file_stat['pass_i_current']]

        #------------------------------------------------------------
        # collect and record statistics of interest
        with open(file_stat_t['cluster_info_stat_file'], 'rb') as f:     cluster_info_stat = pickle.load(f)
        with open(file_stat_t['fsc_stat_file']) as f:     fsc_stat = json.load(f)

        fsc_stat_t = defaultdict(dict)
        for pass_i in fsc_stat:
            for c in fsc_stat[pass_i]:
                fsc_stat_t[int(pass_i)][int(c)] = fsc_stat[pass_i][c]

        fsc_stat = fsc_stat_t
        del fsc_stat_t
        for pass_i in fsc_stat:
            for c in fsc_stat[pass_i]:
                fsc_stat_ = fsc_stat[pass_i][c]
                fsc_stat_['cluster_size'] = cluster_info_stat[fsc_stat_['pass_i']][fsc_stat_['cluster']]['cluster_size']
                fsc_stat_['is_specific'] = (fsc_stat_['is_specific'] is None)

        pass_i_last = N.max([_ for _ in fsc_stat])

        stat[repeat_i] = {'pass_dir':file_stat_t['pass_dir'], 'fsc_stat':fsc_stat[pass_i_last], 'file_stat_file':file_stat['file_stat_file']}

        with open(os.path.join(op['out_dir'], op['stat_out_file']), 'w') as f:       json.dump(stat, f, indent=2)



if __name__ == '__main__':
    main()

