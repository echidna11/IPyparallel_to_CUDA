#!/usr/bin/env python


'''
plot isosurfaces of ground truth, and averages, given a fixed threshold.
notice that for each ground truth, a view angle need to be manually determined, using 
~/ln/tomominer/tomominer/geometry/surface/isosurface/plot_batch.py


it uses output from 
~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_high_fsc_consistency_with_ground_truth/cluster_avg_high_fsc_consistency_with_ground_truth.py



the plot of isosurface is through mayavi

~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_high_fsc_consistency_with_ground_truth/ground_truth_average_isosurface.py

'''

import os, json, pickle

import scipy.stats as SS
import numpy as N

import tomominer.io.file as IF
import tomominer.geometry.rotate as GR
import tomominer.plot.color as PC



def main():
    with open('ground_truth_average_isosurface__op.json') as f:     op = json.load(f)

    with open(op['ground truth']) as f:       true_avg = json.load(f)
    with open(op['ground truth views']) as f:       tv = json.load(f)         # read the view angles of ground truth, for adjusting view angles for plotting

    with open(op['stat file'], 'rb') as f:      st = pickle.load(f)      # this is the stat file dumped by cluster_avg_high_fsc_consistency_with_ground_truth.py

    import tomominer.pursuit.multi.analysis.cluster_avg_high_fsc_consistency_with_ground_truth.cluster_avg_high_fsc_consistency_with_ground_truth as PMACC

    tab = PMACC.stat_table(st)
    PMACC.stat_table__arrange(tab)

    op['isovalue'] = tv[0]['view']['contour']
    st_out = generate_plots(true_avg=true_avg, tab=tab, st=st, view_angs={_['subtomogram']:_ for _ in tv}, op=op)

    with open(op['stat file out'], 'w') as f:       json.dump(st_out, f, indent=2)


import matplotlib.cm as CM
def generate_plots(true_avg, tab, st, view_angs, op):
    op['out dir'] = os.path.abspath(op['out dir'])

    if not os.path.isdir(op['out dir']):      os.makedirs(op['out dir'])

    tab_stat = tab['tab_stat']

    fs = st['fs']

    color = {}
    if 'color in' not in op:
        # if colors are not pre-defined, automatically assign
        cm = PC.get_colormap(len(true_avg) + 1)
        for true_avg_t in true_avg:
            color[true_avg_t['pdb_id']] = cm(true_avg_t['cluster_label'] + 1)[:3]      #   we only get rgb values which are the first three values, not the gamma which is the fourth value
        with open(op['color out'], 'w') as f:       json.dump(color, f, indent=2)
    else:
        with open(op['color in']) as f:     color = json.load(f)
    

    stat = {}

    for true_avg_t in true_avg:
        print true_avg_t['subtomogram']
        assert true_avg_t['subtomogram'] == view_angs[true_avg_t['subtomogram']]['subtomogram'] 
        assert true_avg_t['cluster_label'] == view_angs[true_avg_t['subtomogram']]['cluster_label']

        vt = IF.read_mrc_vol(true_avg_t['subtomogram'])
        if not op['intensity positive']:      vt = -vt

        true_lbl = true_avg_t['cluster_label']

        stat[true_lbl] = {}
        stat[true_lbl]['plot'] = os.path.join(op['out dir'], '%03d--%s.png'%(true_lbl, true_avg_t['pdb_id']))

        if not os.path.isfile(stat[true_lbl]['plot']):
            isosurface(v=vt, isovalue=op['isovalue'], bgcolor=op['bgcolor'], color=color[true_avg_t['pdb_id']], view=view_angs[true_avg_t['subtomogram']]['view'], out_file=stat[true_lbl]['plot'])
        else:
            print 'ignoring', stat[true_lbl]['plot']

        stat[true_lbl]['pred'] = {}
        for pred_lbl in st['selected_templates']:
            if tab_stat[true_lbl][pred_lbl]['count_cluster'] < op['count cluster threshold']:    continue           # if there are very little membership overlap, ignore

            print st['selected_templates'][pred_lbl]['subtomogram']

            v = IF.read_mrc_vol(st['selected_templates'][pred_lbl]['subtomogram'])
            if not op['intensity positive']:      v = -v

            align_t = fs[true_lbl][pred_lbl]['align']
            assert align_t['v1_key']['subtomogram'] == true_avg_t['subtomogram']
            assert align_t['v2_key']['subtomogram'] == st['selected_templates'][pred_lbl]['subtomogram']
            vr = GR.rotate_pad_mean(v, angle=align_t['angle'], loc_r=align_t['loc'])

            isovalue_guess = guess_isovalue(v=vr, m=vt>op['isovalue'])
            if 'dust size cutoff' in op:
                # remove dust if any
                md = dust_identify(vr>=isovalue_guess, size=op['dust size cutoff'])
                vr[md] = isovalue_guess - N.abs(vr).max()*0.01

            stat[true_lbl]['pred'][pred_lbl] = {}

            stat[true_lbl]['pred'][pred_lbl]['plot'] = out_file=os.path.join(op['out dir'], '%03d--%s--%03d.png'%(true_lbl, true_avg_t['pdb_id'], pred_lbl))
            if not os.path.isfile(stat[true_lbl]['pred'][pred_lbl]['plot']):
                isosurface(v=vr, isovalue=isovalue_guess, color=color[true_avg_t['pdb_id']], bgcolor=op['bgcolor'], view=view_angs[true_avg_t['subtomogram']]['view'], out_file=stat[true_lbl]['pred'][pred_lbl]['plot'])
            else:
                print 'ignoring', stat[true_lbl]['pred'][pred_lbl]['plot']

            # in some cases, the pattern is a mixture of multiple true complexes, we may use a grey version of the plot instead
            stat[true_lbl]['pred'][pred_lbl]['plot-grey'] = out_file=os.path.join(op['out dir'], '%03d--%s--%03d--grey.png'%(true_lbl, true_avg_t['pdb_id'], pred_lbl))
            if not os.path.isfile(stat[true_lbl]['pred'][pred_lbl]['plot-grey']):
                isosurface(v=vr, isovalue=isovalue_guess, color=op['grey color'], bgcolor=op['bgcolor'], view=view_angs[true_avg_t['subtomogram']]['view'], out_file=stat[true_lbl]['pred'][pred_lbl]['plot-grey'])
            else:
                print 'ignoring', stat[true_lbl]['pred'][pred_lbl]['plot-grey']
           

            print true_lbl, pred_lbl


    return stat


'''
find optimal isovalue so that its cutoff is most consistant with the mask
'''
def guess_isovalue(v, m, bin_num=100):
    m = m.astype(N.float).flatten()

    best = None
    ts = N.linspace(start=v.min(), stop=v.max(), num=bin_num)
    
    for t in ts[1:-1]:
        vm = (v>t).astype(N.float).flatten()

        cor = SS.pearsonr(vm, m)[0]

        if (best is None) or (best['cor'] < cor):       best = {'cor':cor, 'isovalue':t}

    return best['isovalue']


from mayavi import mlab
def isosurface(v, isovalue, bgcolor, color, view, out_file, show=False):
    fig = mlab.figure(bgcolor=tuple(bgcolor), size=(400, 400))
    src = mlab.pipeline.scalar_field(v)
    mlab.pipeline.iso_surface(src, contours=[isovalue,], color=tuple(color))

    if False:
        mlab.view(azimuth=view['azimuth'], elevation=view['elevation'], distance=view['distance'], focalpoint=N.array(view['focalpoint']), roll=view['roll'])
    else:
        # somehow distance must set after change angle, otherswise some structures cannot be correctly displayed, don't know why
        mlab.view(azimuth=view['azimuth'], elevation=view['elevation'], roll=view['roll'])
        mlab.view(distance=view['distance'], focalpoint=N.array(view['focalpoint']))

    if show:        mlab.show()
    mlab.savefig(filename=out_file)
    mlab.close()



import tomominer.segmentation.util as SU
'''
find dust defined by smaller than certain size
'''
def dust_identify(m, size):

    sc = SU.connected_components(m)

    s = N.zeros(sc['label'].shape, dtype=bool)
    for l in range(1, sc['num']+1):
        if (sc['label'] == l).sum() > size:     continue
        s[sc['label'] == l] = True

    return s

if __name__ == '__main__':
    main()

