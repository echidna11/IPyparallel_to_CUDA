#!/usr/bin/env python



'''
given output from config_prepare.py, batch calculate gold standard FSC

~/ln/tomominer/tomominer/pursuit/multi/analysis/statistics/gold_standard_fsc/batch_run.py

'''


import os, json, pickle, random, copy


from tomominer.pursuit.multi.analysis.cluster_avg_fsc_curve.half_avg_independent import half_average, calculate_fsc
from tomominer.parallel.queue_master import QueueMaster
import tomominer.common.obj as CO



def process(self, avg_op, conf_stat):

    fsc = {}
    for i, cs in enumerate(conf_stat):
        with open(cs['split'], 'rb') as f:      djs = pickle.load(f)
        
        print 'avg', cs['template']['subtomogram'], 'half sizes', {_:len(djs[_]) for _ in djs}

        # calculate standard FSC through spliting to two halves, keep subtomograms aligned
        ha = half_average(self=self, djs=djs, out_dir=cs['out_dir'], avg_op=avg_op)
        if len(ha) < 2:     continue

        fsc[i] = [float(_) for _ in calculate_fsc(ha, align_op=avg_op['align'], out_dir=cs['out_dir'])]
        print 'fsc', fsc[i]


    return fsc


def main():
    with open('batch_run__op.json') as f:       op = json.load(f)

    with open(op['config_stat']) as f:      conf_stat = json.load(f)

    with open(op['averaging_op']) as f:      avg_op = json.load(f)

    # just get a generic object that can contain a cache
    self = CO.Object()

    self.runner = QueueMaster(host='localhost', port=5011)

    fsc = process(self, avg_op=avg_op, conf_stat=conf_stat)

    with open(op['stat_out'], 'w') as f:      json.dump({'fsc':fsc}, f, indent=2)


if __name__ == '__main__':
    main()


'''
related code
~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_fsc_curve/half_avg_independent.py
'''


