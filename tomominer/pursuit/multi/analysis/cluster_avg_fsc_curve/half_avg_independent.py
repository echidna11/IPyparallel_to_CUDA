#!/usr/bin/env python


'''
split a cluster into two halves then calculate averages of aligned subtomograms

~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_fsc_curve/half_avg_independent.py
'''



import os
import json
import cPickle as pickle
import copy
import numpy as N

import tomominer.io.file as IF
import tomominer.statistics.vol as SV
import tomominer.common.obj as CO
from tomominer.parallel.queue_master import QueueMaster

import tomominer.average.main_routine as AA
import tomominer.pursuit.multi.util as PMU
import tomominer.geometry.rotate as GR

def half_average(self, djs, out_dir, avg_op):

    ha = {}
    for l in djs:
        out_dir_t = os.path.join(out_dir, str(l))
        if not os.path.isdir(out_dir_t):        os.makedirs(out_dir_t)

        # make a copy of data json and REMOVE existing alignment information
        dj = []
        for d in djs[l]:            dj.append(     {'subtomogram':d['subtomogram'], 'mask':d['mask']}     )

        tk = AA.do_average(self=self, op=avg_op, data_json=dj, out_dir=out_dir_t)
        if tk is None:      continue

        ha[int(l)] = tk

    return ha
    


def calculate_fsc(ha, align_op, out_dir):
    v0 = IF.read_mrc(ha[0]['subtomogram'])['value'].astype(N.float)
    m0 = IF.read_mrc(ha[0]['mask'])['value'].astype(N.float)

    v1 = IF.read_mrc(ha[1]['subtomogram'])['value'].astype(N.float)
    m1 = IF.read_mrc(ha[1]['mask'])['value'].astype(N.float)

    if True:
        a = PMU.align_vols_with_wedge(v0, m0, v1, m1, op={'L':align_op['L']})
    else:
        a = PMU.align_vols_no_wedge(v0, v1, op={'L':align_op['L']})     # Observation: alignment without wedge gives better FSC. But currently this function is not working because corresponding interfaces in cython is not added

    if a['err'] is not None:    raise Exception(a['err'])



    v1r = GR.rotate_pad_mean(v1, angle=a['angle'], loc_r=a['loc'])
    m1r = GR.rotate_mask(m1, angle=a['angle'])

    common_frame_dir = os.path.join(out_dir, 'common_frame')
    if not os.path.isdir(common_frame_dir):     os.makedirs(common_frame_dir)

    IF.put_mrc(v0, os.path.join(common_frame_dir, 'avg-0.mrc'), overwrite=True)
    IF.put_mrc(m0, os.path.join(common_frame_dir, 'mask-0.mrc'), overwrite=True)
    IF.put_mrc(v1r, os.path.join(common_frame_dir, 'avg-1.mrc'), overwrite=True)
    IF.put_mrc(m1r, os.path.join(common_frame_dir, 'mask-1.mrc'), overwrite=True)

    return SV.fsc(v0, v1r)


def process(self, op):

    djs = op['djs']

    fsc = {}
    for c in djs:

        if ('selected_clusters' in op) and (int(c) not in set(op['selected_clusters'])):      continue

        print 'cluster', c, 'half sizes', {_:len(djs[c][_]) for _ in djs[c]}

        cluster_dir = os.path.join(op['out_dir'], str(c))
       
        # calculate standard FSC through spliting to two halves, keep subtomograms aligned
        ha = half_average(self=self, djs=djs[c], out_dir=cluster_dir, avg_op=op)
        if len(ha) < 2:     continue

        fsc[c] = [float(_) for _ in calculate_fsc(ha, align_op=op['align'], out_dir=cluster_dir)]
        print 'fsc', fsc[c]


    return fsc


def main():
    with open('half_avg_independent__op.json') as f:    op = json.load(f)
    op['out_dir'] = os.path.abspath(op['out_dir'])


    # just get a generic object that can contain a cache
    self = CO.Object()

    self.runner = QueueMaster(host='localhost', port=5011)

    with open(op['data']) as f:   op['djs'] = json.load(f)

    fsc = process(self, op)

    with open(os.path.join(op['out_dir'], 'half_avg_independent__out.json'), 'w') as f:      json.dump({'fsc':fsc}, f, indent=2)


if __name__ == '__main__':
    main()

