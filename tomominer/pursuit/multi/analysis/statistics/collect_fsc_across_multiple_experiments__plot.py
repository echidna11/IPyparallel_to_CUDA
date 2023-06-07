#!/usr/bin/env python



'''
plot fsc curves

~/ln/tomominer/tomominer/pursuit/multi/analysis/statistics/collect_fsc_across_multiple_experiments__plot.py
'''


import os, json
import pickle

import tomominer.pursuit.multi.analysis.cluster_avg_fsc_curve.plot as PMACP


def main():
    with open('collect_fsc_across_multiple_experiments__plot__op.json') as f:       op = json.load(f)

    op['out_dir'] = os.path.abspath(op['out_dir'])
    if not os.path.isdir(op['out_dir']):        os.makedirs(op['out_dir'])

    with open(op['fsc_collect_file'], 'rb') as f:         co = pickle.load(f)


    # collect file names of averages, find common prefix
    ks = []
    for fn in co:
        for c in co[fn]:
            ks.append(co[fn][c]['template']['subtomogram'])

    common_prefix = os.path.commonprefix(ks)


    stat = {}       # just record the file name of average, and corresponding figure file

    for fn in co:
        for c in co[fn]:
            sfn_key = co[fn][c]['template']['subtomogram']
            sfn = sfn_key[len(common_prefix):]
            sfn = sfn.replace('/', '--')
            sfn = sfn.replace('.', '--')
            out_file = os.path.join(op['out_dir'], sfn+'.eps')
            pr = PMACP.plot_one(fsc=co[fn][c]['fsc'], voxel_spacing=op['voxel_spacing'], fsc_cutoff=(op['fsc_cutoff'] if 'fsc_cutoff' in op else None), x_tick_rotation=op['x_tick_rotation'], out_file=out_file)

            stat[sfn_key] = {'file': out_file}

            print sfn


    with open(op['stat_out'], 'w') as f:        json.dump(stat, f, indent=2)



if __name__ == '__main__':
    main()




'''

related code

~/ln/tomominer/tomominer/pursuit/multi/analysis/statistics/collect_fsc_across_multiple_experiments.py

~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_ssnr_fsc__plot.py

'''

