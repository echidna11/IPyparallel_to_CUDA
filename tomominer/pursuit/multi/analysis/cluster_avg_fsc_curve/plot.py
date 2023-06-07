#!/usr/bin/env python


'''
plot the FSC curve of a given set of cluster averages, save these curves as EPS images

~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_fsc_curve/plot.py
'''


import os
import json
import cPickle as pickle
import shutil
import matplotlib.pyplot as PLT
import matplotlib.cm as MCM
import pylab
import tomominer.io.file as IF
import tomominer.statistics.vol as SV
import tomominer.common.obj as CO
from tomominer.parallel.queue_master import QueueMaster
import tomominer.pursuit.multi.util as PMU

def plot_all(ssnr_fsc, half_fsc, ind_half_fsc, voxel_spacing, out_file=None):

    assert      len(ssnr_fsc) == len(half_fsc)
    assert      len(ssnr_fsc) == len(ind_half_fsc)

    # get resolution list
    r = SV.fft_resolution_list(band_num=len(ssnr_fsc), voxel_spacing=voxel_spacing)[0].flatten()

    r_s = []
    for r_t in r:   r_s.append('1/%.1f'%(r_t,))

    # plot curve and save
    fig = PLT.figure()
    #fig.suptitle('bold figure suptitle', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)
    #ax.set_title('Fourier Shell Correlation')

    ax.plot(ssnr_fsc, '-', label='SSNR FSC')
    ax.plot(half_fsc, '--', label='HALF FSC')
    ax.set_xlabel('Spatial Frequency (1/nm)')
    ax.set_ylabel('Fourier Shell Correlation')
    ax.set_ylim(ymin=0, ymax=1.1)

    PLT.legend( ('SSNR FSC', 'HALF FSC') )

    PLT.xticks(range(len(ssnr_fsc)), r_s)

    for tick in ax.xaxis.get_major_ticks():
        # see http://stackoverflow.com/questions/6390393/matplotlib-make-tick-labels-font-size-smaller
        tick.label.set_fontsize('x-small')
        tick.label.set_rotation('vertical')


    PLT.draw()

    PLT.savefig(out_file, bbox_inches='tight')
    
    PLT.close("all")


def plot_one(fsc, voxel_spacing, x_tick_step=1, x_tick_rotation=None, fsc_cutoff=[0.5], out_file=None):
    
    print 'voxel_spacing', voxel_spacing
    if fsc_cutoff is not None: 
        # print out resolution
        resolution_s = []
        for fsc_cutoff_t in fsc_cutoff:
            band_i = SV.band_interpolation(fsc, cutoff=fsc_cutoff_t)
            resolution = SV.band_to_resolution(band_num=len(fsc), voxel_spacing=voxel_spacing, band=band_i)
            print 'cutoff', fsc_cutoff_t, 'resolution', resolution
            resolution_s.append({'fsc_cutoff':fsc_cutoff_t, 'band_i':band_i, 'resolution':resolution})
    else:
        resolution_s = None


    # get resolution list
    r = SV.fft_resolution_list(band_num=len(fsc), voxel_spacing=voxel_spacing)[0]
    r_s = []
    for r_t in r:   r_s.append('1/%.1f'%(r_t,))

    # plot curve and save
    fig = PLT.figure()
    #fig.suptitle('bold figure suptitle', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)
    #ax.set_title('Fourier Shell Correlation')

    ax.plot(fsc, '-', label='FSC')
    if fsc_cutoff is not None:
        for rt in resolution_s:
            ax.plot([0, rt['band_i']], [rt['fsc_cutoff'], rt['fsc_cutoff']], 'g--')
            ax.plot([rt['band_i'], rt['band_i']], [rt['fsc_cutoff'], 0], 'g--')

    ax.set_xlabel('Spatial Frequency (1/nm)')
    ax.set_ylabel('Fourier Shell Correlation')
    ax.set_ylim(ymin=0, ymax=1.1)

    #ax.set_xticks(range(len(fsc)), r_s)
    #ax.set_xticklabels(r_s, size='x-small')
    if x_tick_rotation is None:     x_tick_rotation = 0
    x_i = range(0, len(fsc), x_tick_step)
    PLT.xticks(x_i, [r_s[_] for _ in x_i], size='medium', rotation=x_tick_rotation)

    PLT.draw()
    PLT.savefig(out_file, bbox_inches='tight')
    PLT.close('all')

    re = {}
    if resolution_s is not None:        re['resolutions'] = resolution_s
    return re

if __name__ == '__main__':

    with open('plot__op.json') as f:       op = json.load(f)

    with open(op['tomogram_info_file'], 'rb') as f:     tom_info = json.load(f)
    tom_info = {  _['id']:_ for _ in tom_info   }

    with open(op['half_avg_independent__out']) as f:   hai_fsc = json.load(f)['fsc']
    
    out_dir = os.path.abspath(op['out_dir'])
    if os.path.isdir(out_dir):     shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    r = []
    for c in hai_fsc:
        print c
        resolutions = plot_one(fsc=hai_fsc[c], voxel_spacing=tom_info[op['tomogram_id']]['voxel_spacing'], x_tick_step=op['x_tick_step'], x_tick_rotation=(op['x_tick_rotation'] if 'x_tick_rotation' in op else None), fsc_cutoff=op['fsc_cutoff'], out_file=os.path.join(out_dir, 'plot-%d.eps'%(int(c),)))
        r.append(   {'cluster':c, 'resolutions':resolutions}  )

    with open(op['stat_out'], 'w') as f:      json.dump(r, f, indent=2)


