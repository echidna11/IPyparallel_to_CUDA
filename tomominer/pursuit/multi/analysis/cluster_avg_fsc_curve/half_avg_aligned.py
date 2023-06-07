#!/usr/bin/env python

# split a cluster into two halves then calculate averages of aligned subtomograms


import os
import json
import cPickle as pickle

import tomominer.io.file as IF
import tomominer.statistics.vol as SV
import tomominer.common.obj as CO
from tomominer.parallel.queue_master import QueueMaster

import tomominer.pursuit.multi.util as PMU



def half_average(self, djs, out_dir, avg_op):

    labels = []
    dj = []
    for l in djs:
        dj.extend(djs[l])
        labels.extend([int(l)] * len(djs[l]))


    a = PMU.cluster_averaging(self=self, data_json=dj, labels=labels, n_chunk=10, out_dir=out_dir, op=avg_op, centerize_loc=True)

    return a['template_keys']


def calculate_fsc(ha):
    v0 = IF.read_mrc(ha[0]['subtomogram'])['value'];       os.remove(ha[0]['subtomogram']);          os.remove(ha[0]['mask'])
    v1 = IF.read_mrc(ha[1]['subtomogram'])['value'];       os.remove(ha[1]['subtomogram']);          os.remove(ha[1]['mask'])

    return SV.fsc(v0, v1)




if __name__ == '__main__':

    with open('half_avg_aligned__config.json') as f:    op = json.load(f)

    # just get a generic object that can contain a cache
    self = CO.Object()

    self.runner = QueueMaster(host='localhost', port=5011)

    with open('half_split__out.json') as f:   djs = json.load(f)

    fsc = {}
    for c in djs:
        print c
       
        # calculate standard FSC through spliting to two halves, keep subtomograms aligned
        ha = half_average(self=self, djs=djs[c], out_dir=op['tmp_dir'], avg_op=op['averaging'])

        fsc[c] = [float(_) for _ in calculate_fsc(ha)]

    with open('half_avg_aligned__out.json', 'w') as f:      json.dump({'fsc':fsc}, f, indent=2)

