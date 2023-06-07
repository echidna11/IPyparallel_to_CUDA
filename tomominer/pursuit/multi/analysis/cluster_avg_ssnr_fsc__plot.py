#!/usr/bin/env python

'''
plot SSNR based FSC curves
~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_ssnr_fsc__plot.py
'''

import os, json
import cPickle as pickle
import numpy as N

def main():
    with open('cluster_avg_ssnr_fsc__plot__op.json') as f:     op = json.load(f)

    op['out_dir'] = os.path.abspath(op['out_dir'])
    if not os.path.isdir(op['out_dir']):        os.makedirs(op['out_dir'])

    with open(op['file_stat_file']) as f:         fs = json.load(f)
    fsp = fs['passes']
    fsp = {int(_):fsp[_] for _ in fsp}

    with open(fsp[fs['pass_i_current']]['fsc_stat_file']) as f:       fsc_stat = json.load(f)
    from tomominer.pursuit.multi.main import fsc_stat_json_convert
    fsc_stat = fsc_stat_json_convert(fsc_stat)

    ci = {}
    for pass_i in fsp:
        with open(fsp[pass_i]['cluster_info_file'], 'rb') as f:     ci[pass_i] = pickle.load(f)


    stat = {}

    import tomominer.pursuit.multi.analysis.cluster_avg_fsc_curve.plot as PMACP
    pass_i = op['pass_i'] if 'pass_i' in op else fs['pass_i_current']
    for c in fsc_stat[pass_i]:
        cit = ci[fsc_stat[pass_i][c]['pass_i']][fsc_stat[pass_i][c]['cluster']]

        assert      fsc_stat[pass_i][c]['pass_i'] == cit['pass_i']
        assert      fsc_stat[pass_i][c]['cluster'] == cit['cluster']

        out_file = os.path.join(op['out_dir'], '%03d.eps'%(c,))
        pr = PMACP.plot_one(fsc=cit['fsc'], voxel_spacing=op['voxel_spacing'], fsc_cutoff=(op['fsc_cutoff'] if 'fsc_cutoff' in op else None), x_tick_rotation=op['x_tick_rotation'], out_file=out_file)
        print 'selected', 'pass', pass_i, 'cluster', c, 'original', 'pass', cit['pass_i'], 'cluster', cit['cluster'], 'fsc_sum', N.sum(cit['fsc']), 'resolution', (pr['resolution'] if 'resolution' in pr else None)

        stat[c] = {'subtomogram':cit['template_key']['subtomogram'], 'file':out_file, 'fsc_sum':N.sum(cit['fsc']), 'template_key':cit['template_key']}
        if 'resolution' in pr:      stat[c]['resolution'] = pr['resolution']

    with open(op['stat_out'], 'w') as f:       json.dump(stat, f, indent=2)



if __name__ == '__main__':
    main()

'''
related code
~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_ssnr_fsc_resolution.py
~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_fsc_curve/plot.py
'''
