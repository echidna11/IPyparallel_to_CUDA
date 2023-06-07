#!/usr/bin/env python



# parallel calculation fo pairwise alignment between all subtomograms, in parallel


import os
import sys
import json

import classify.util as CU
from tomominer.parallel.queue_master import QueueMaster


def align(dj):
    runner = QueueMaster('localhost', 5011)


    align_op = {'with_missing_wedge':True, 'L':36}
    tasks = []
    for i1,v1k in enumerate(dj):
        for i2,v2k in enumerate(dj):
            if i1 == i2:    continue
            tasks.append(runner.task(module='tomominer.classify.util', method='align_keys', kwargs={'v1k':v1k, 'v2k':v2k, 'op':align_op}))


    for res in runner.run__except(tasks):
        re = res.result
        print '%d %d %f'%(re['v1_key']['id'], re['v2_key']['id'], re['score'])



if __name__ == '__main__':

    data_json_file = sys.argv[1]
    with open(data_json_file) as f:     data_json = json.load(f)

    for i, r in enumerate(data_json):
        assert 'id' not in r
        r['id'] = i

    id_file = sys.argv[2]
    with open(id_file, 'w') as f:       json.dump(data_json, f, indent=2)

    align(dj=data_json)


