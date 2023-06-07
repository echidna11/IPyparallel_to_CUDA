#!/usr/bin/env python



'''
use multiprocessing to run multiple MPP instances

~/ln/tomominer/tomominer/pursuit/multi/batch/pose_norm/run.py
'''




import os, sys, json, copy, warnings
from multiprocessing.pool import Pool as Pool

from tomominer.parallel.queue_master import QueueMaster

import tomominer.common.obj as CO
from tomominer.io.cache import Cache


def run_single(op, uuid):

    warnings.filterwarnings('error')
    
    # just get a generic object that can contain a cache
    self = CO.Object()

    self.pool = None
  
    self.cache = Cache(tmp_dir=op['options']['tmp_dir'])

    self.runner = QueueMaster(op['options']['network']['qhost'], op['options']['network']['qport'])

    with open(op['data_file']) as f:        dj = json.load(f)

    # convert relative path to absolute path, if needed
    for d in dj:
        if not os.path.isabs(d['subtomogram']):     d['subtomogram'] = os.path.abspath(os.path.join(os.path.dirname(op['data_file']), d['subtomogram']))
        if not os.path.isabs(d['mask']):     d['mask'] = os.path.abspath(os.path.join(os.path.dirname(op['data_file']), d['mask']))


    from tomominer.pursuit.multi.main import pursuit
    fs = pursuit(self=self, op=op, data_json=dj)['file_stat_file']

    return {'uuid':uuid, 'file_stat':fs}


def main():
    with open('run__op.json') as f:          op = json.load(f)
    with open(op['config prepare stat']) as f:      cps = json.load(f)

    with open(op['pursuit op']) as f:           pop = json.load(f)
    assert  'data_file' not in pop
    assert  'out_dir' not in pop

    #pool = Pool()
    results = []

    for grp_i in cps:
        for rep_i in cps[grp_i]['repeat']:
            pop_t = copy.deepcopy(pop)
            pop_t['data_file'] = cps[grp_i]['data']
            pop_t['out_dir'] = cps[grp_i]['repeat'][rep_i]['out_dir']
            if not os.path.isdir(pop_t['out_dir']):         os.makedirs(pop_t['out_dir'])
            results.append(run_single(op=pop_t, uuid=cps[grp_i]['repeat'][rep_i]['uuid']))
            #pre.append(pool.apply_async(func=run_single, kwds={'op':pop_t, 'uuid':cps[grp_i]['repeat'][rep_i]['uuid']}))

    #st = []
    #for i, r in enumerate(pre):
        #st.append(r.get(99999999))

        #print '\r ', i, float(i) / len(pre),
        #sys.stdout.flush()

    with open(op['stat out'], 'w') as f:        json.dump(results, f, indent=2)
    

if __name__ == '__main__':
    main()



'''
related code

/home/rcf-47/mxu/ln/tomominer/tomominer/average/genetic_algorithm/config/from_multi_pursuit__batch_run.py

~/ln/tomominer/tomominer/pursuit/multi/run.py

'''


