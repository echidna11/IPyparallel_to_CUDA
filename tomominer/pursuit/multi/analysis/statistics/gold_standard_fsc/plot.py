#!/usr/bin/env python




'''


~/ln/tomominer/tomominer/pursuit/multi/analysis/statistics/gold_standard_fsc/plot.py
'''


import os, json
import pickle

import tomominer.pursuit.multi.analysis.cluster_avg_fsc_curve.plot as PMACP


def main():
    with open('plot__op.json') as f:       op = json.load(f)

    op['out_dir'] = os.path.abspath(op['out_dir'])
    if not os.path.isdir(op['out_dir']):        os.makedirs(op['out_dir'])

    with open(op['config_stat_file']) as f:         conf_stat = json.load(f)
    with open(op['fsc_stat_file']) as f:         fsc_stat = json.load(f)['fsc']           # normally  batch_run__out_stat.json
    fsc_stat = {int(_):fsc_stat[_] for _ in fsc_stat}

    # collect file names of averages, find common prefix
    ks = []
    for cs in conf_stat:
        ks.append(cs['template']['subtomogram'])

    common_prefix = os.path.commonprefix(ks)


    stat = {}       # just record the file name of average, and corresponding figure file
    for cs in conf_stat:
        sfn_key = cs['template']['subtomogram']
        sfn = sfn_key[len(common_prefix):]
        sfn = sfn.replace('/', '--')
        sfn = sfn.replace('.', '--')
        sfn = '%04d----%s'%(cs['count'], sfn)

        out_file = os.path.join(op['out_dir'], sfn+'.eps')
        pr = PMACP.plot_one(fsc=fsc_stat[cs['count']], voxel_spacing=op['voxel_spacing'], fsc_cutoff=(op['fsc_cutoff'] if 'fsc_cutoff' in op else None), x_tick_rotation=op['x_tick_rotation'], out_file=out_file)

        stat[cs['count']] = {'file': out_file, 'plot_res':pr}

        print sfn


    with open(op['stat_out'], 'w') as f:        json.dump(stat, f, indent=2)



if __name__ == '__main__':
    main()






'''
replated code

~/ln/tomominer/tomominer/pursuit/multi/analysis/statistics/collect_fsc_across_multiple_experiments__plot.py

'''

