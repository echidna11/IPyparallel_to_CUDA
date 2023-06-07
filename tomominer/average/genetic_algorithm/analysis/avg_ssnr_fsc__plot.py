#!/usr/bin/env python






'''
plot fsc curve of a particular average

~/ln/tomominer/tomominer/average/genetic_algorithm/analysis/avg_ssnr_fsc__plot.py
'''



import os, json
import cPickle as pickle



def main():
    with open('avg_ssnr_fsc__plot__op.json') as f:     op = json.load(f)


    with open(op['file_stat']) as f:     fs = json.load(f)
    pass_i = op['pass_i'] if 'pass_i' in op else fs['best_pass_i']
    print 'choose pass_i', pass_i


    ps = fs['passes']
    ps = {int(_):ps[_] for _ in ps}
    ps = ps[pass_i]

    with open(ps['pmpg_file'], 'rb') as f:       pmpg = pickle.load(f)
    with open(ps['cluster_averaging_file'], 'rb') as f:       cas = pickle.load(f)


    import tomominer.pursuit.multi.analysis.cluster_avg_fsc_curve.plot as PMACP

    out_file = os.path.abspath(op['plot_file_out'])
    pr = PMACP.plot_one(fsc=pmpg['best']['e'][0]['fsc'], voxel_spacing=op['voxel_spacing'], fsc_cutoff=op['fsc_cutoff'], x_tick_rotation=op['x_tick_rotation'], out_file=out_file)
    print 'fsc', pmpg['best']['e'][0]['fsc']

    tk = cas['template_keys'][0]

    stat = {}
    stat[0] = {'subtomogram':tk['subtomogram'], 'file':out_file, 'resolution':pr['resolution'], 'template_key':tk}

    with open(op['stat_out'], 'w') as f:       json.dump(stat, f, indent=2)



if __name__ == '__main__':
    main()


'''

code modified from

~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_ssnr_fsc__plot.py

'''

