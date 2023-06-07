#!/usr/bin/env python


'''
for the repeated tests, calculate Gold Standard FSC for all resulting clusters at last step, using code inside following programs

~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_fsc_curve/half_split.py
~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_fsc_curve/half_avg_independent.py

# todo: for resolutions, add both cutoff 0.5 and 0.143

'''


import os, json, copy
import cPickle as pickle
import numpy as N

from tomominer.parallel.queue_master import QueueMaster
import tomominer.common.obj as CO

import tomominer.pursuit.multi.analysis.cluster_avg_fsc_curve.half_split as PMACHS
import tomominer.pursuit.multi.analysis.cluster_avg_fsc_curve.half_avg_independent as PMACHAI
import tomominer.pursuit.multi.analysis.cluster_avg_fsc_curve.plot as PMACP


def process(self, op):
    op['out_dir'] = os.path.abspath(op['out_dir'])

    with open(op['stat_file']) as f:    rep_stat = json.load(f)
    rep_stat = {int(_):rep_stat[_] for _ in rep_stat}

    with open(op['plot']['tomogram_info_file'], 'rb') as f:     tom_info = json.load(f)
    tom_info = {  _['id']:_ for _ in tom_info   }
    voxel_spacing = tom_info[op['plot']['tomogram_id']]['voxel_spacing']

    stat = {}
    for rep_i in rep_stat:
        print '---------------------------------'
        print 'repeat', rep_i

        stat[rep_i] = {}

        with open(rep_stat[rep_i]['file_stat_file']) as f:    fs = json.load(f)
        fs['passes'] = {int(_):fs['passes'][_] for _ in fs['passes']}
        fsp = fs['passes'][fs['pass_i_current']]
        stat[rep_i]['average'] = {}
        stat[rep_i]['average']['file_stat_file'] = rep_stat[rep_i]['file_stat_file']
        stat[rep_i]['average']['pass_i_used'] = fsp['pass_i']

        stat[rep_i]['out_dir'] = os.path.join(op['out_dir'], '%04d'%(rep_i,))
        if not os.path.isdir(stat[rep_i]['out_dir']):            os.makedirs(stat[rep_i]['out_dir'])

        hai_op = copy.deepcopy(op['half_avg_independent'])
        
        #-------------------------------------------
        # split data
        stat[rep_i]['half_split__out__file'] = os.path.join(stat[rep_i]['out_dir'], 'half_split__out.json')
        print 'calculating', stat[rep_i]['half_split__out__file']

        if os.path.isfile(stat[rep_i]['half_split__out__file']):
            with open(stat[rep_i]['half_split__out__file']) as f:       hai_op['djs'] = json.load(f)

        else:
            with open(fsp['cluster_average_select_file'], 'rb') as f:   cas = pickle.load(f)
            st = cas['selected_templates']
            ti = cas['tk_info']
            del cas

            ti = {_:ti[st[_]['subtomogram']]   for _ in st}
            hai_op['djs'] = PMACHS.split(ti)
            with open(stat[rep_i]['half_split__out__file'], 'w') as f:  json.dump(hai_op['djs'], f, indent=2)

        del fsp

        #------------------------------------------
        # calculate FSC
        hai_op['out_dir'] = os.path.join(stat[rep_i]['out_dir'], 'avg')
        if not os.path.isdir(hai_op['out_dir']):            os.makedirs(hai_op['out_dir'])

        hai_op['out_file'] = os.path.join(hai_op['out_dir'], 'half_avg_independent__out.json')
        print 'calculating', hai_op['out_file']

        if os.path.isfile(hai_op['out_file']):
            with open(hai_op['out_file']) as f:        fsc = json.load(f)

        else:
            fsc = PMACHAI.process(self, hai_op)
            with open(hai_op['out_file'], 'w') as f:        json.dump(fsc, f, indent=2)

        stat[rep_i]['fsc_file'] = hai_op['out_file']

        del hai_op

        #--------------------------------------------
        # plot and calculate resolution
        stat[rep_i]['plot'] = {}
        stat[rep_i]['plot']['dir'] = os.path.join(stat[rep_i]['out_dir'], 'plot')
        if not os.path.isdir(stat[rep_i]['plot']['dir']):            os.makedirs(stat[rep_i]['plot']['dir'])

        stat[rep_i]['plot']['stat_file'] = os.path.join(stat[rep_i]['plot']['dir'], 'stat.json')
        print 'calculating', stat[rep_i]['plot']['stat_file']

        if os.path.isfile(stat[rep_i]['plot']['stat_file']):
            with open(stat[rep_i]['plot']['stat_file']) as f:   plot_stat = json.load(f)
        else:
            plot_stat = {}
            for c in fsc:
                if not N.all(N.isfinite(fsc[c])):   continue

                plot_stat[c] = {}
                plot_stat[c]['plot_file'] = os.path.join(stat[rep_i]['plot']['dir'], 'plot-%d.eps'%(int(c),))
                plot_stat[c]['resolution'] = PMACP.plot_one(fsc=fsc[c], voxel_spacing=voxel_spacing, x_tick_step=op['plot']['x_tick_step'], out_file=plot_stat[c]['plot_file'])
            with open(stat[rep_i]['plot']['stat_file'], 'w') as f:   json.dump(plot_stat, f, indent=2)

        stat[rep_i]['plot']['stat'] = plot_stat


        del plot_stat

    return stat


if __name__ == '__main__':

    self = CO.Object()
    self.runner = QueueMaster(host='localhost', port=5011)

    with open('half_avg_independent__op.json') as f:    op = json.load(f)
    stat = process(self, op)
    with open(op['stat_file_out'], 'w') as f:    json.dump(stat, f, indent=2)


