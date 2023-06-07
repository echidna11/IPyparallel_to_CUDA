

# given a set of subtomograms, test the robustness of averaging. 
import os

import json
import copy

import numpy as np

import classify.config as config

from work_queue.queue_master import QueueMaster

import tomominer.average.average as avg

# test whether can get same structure given different random initial orientations

def random_orientation_test(op, data, out_dir, runner):
    
    for repeat_i in range(op['test']['robust']['orientation']['repeat_num']):

        print '=================================================================================='
        print 'repeat %d'%(repeat_i)

        data_t = [None] * len(data)
        for i, (v,m,a,l) in enumerate(data):
            data_t[i] = ( v, m, np.random.random(3) * (np.pi * 2), np.array([0.0, 0.0, 0.0]) )

        repeat_dir = os.path.join(out_dir, 'rep_%04d'%(repeat_i))
        if not os.path.exists(repeat_dir) : 
            os.mkdir(repeat_dir);   
            #os.chmod(repeat_dir, 0775)    # mxu: enable writing access from other members in the group

        avg.average(op=op, data=data_t, out_dir=repeat_dir, runner=runner, dump_initial_data=True)




# test whether can get same structre given a random subsampling of subtomograms with random initial orientations




# test the change of averaging parameters, see if can still converge


# optional: if you increase or decrease subtomogram size, see if can still converge



if __name__ == '__main__':
    import sys
    qhost = sys.argv[1]
    qport = 5011

    param_file = sys.argv[2]
    data_file  = sys.argv[3]
    out_dir = sys.argv[4]

    with open(param_file) as f:
        op = json.load(f)

    with open(data_file) as f:
        data_json = json.load(f)
    

    data = config.parse_data(data_json)
    
    runner = QueueMaster(qhost, qport)

    random_orientation_test(op, data, out_dir, runner)

