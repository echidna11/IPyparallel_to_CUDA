#!/usr/bin/env python



# given option file and a data json file, for each subtomogram in the json file, identify and exclude its neighbor complexes, then save the excluded subtomogram
import os, errno
import sys
import json
from multiprocessing.pool import Pool
import numpy as N
import scipy.ndimage.morphology as SNM

import tomominer.io.file as IF
import tomominer.filter.gaussian as FG
import tomominer.segmentation.active_contour.chan_vese.segment as SA
import tomominer.segmentation.watershed as SW
import tomominer.segmentation.connected_regions as SC
import tomominer.image.vol.util as CV
import tomominer.io.path as IP


# for a Gausaian smoothed subtomogram, use chan-vese to perform segmentation. Then identify those segments that touches the max ball of the subtomogram. Use watersheld segmentation on phi to exclude such regions
def neighbor_exclude(record, op):
    v_org = IF.get_mrc(record['subtomogram'])
    v = v_org - v_org.mean()

    # first perform segmentation
    vg = FG.dog_smooth(v, s1=op['smooth']['sigma1'], s2=op['smooth']['sigma1']*op['smooth']['sigma_k'])

    phi = N.sign(vg - N.abs(vg).mean())
    phi = SA.segment(vg, phi, op['ac_segmentation']['smooth_weight'], op['ac_segmentation']['image_weight'] / vg.var())

    if vg[phi > 0].mean() < vg[phi < 0].mean():    phi = -phi      # assume bright regions in subtomograms corresponding to high electron densities


    # identify those connected regions that touches the max ball of the subtomogram
    seg_msk = N.zeros(phi.shape, dtype=N.int);        seg_msk[phi > 0] = 1
    conn_rgn_lbl = SC.connected_regions(seg_msk);       conn_rgn_lbl = conn_rgn_lbl['lbl']
    if conn_rgn_lbl.max() == 0:     return None     # no connected region found

    
    # use connected regions as seed, perform watershed segmentation on phi
    sws = SW.segment(phi, conn_rgn_lbl)
    sws_seg = sws['vol_seg_lbl'];   sws_seg[sws_seg < 0] = 0

    '''
    # delete the phi-segments of neighbor structures
    for lbl in range(1, sws_seg.max()+1):
        if lbl not in selected_lbls:
            sws_seg[sws_seg == lbl] = 0

    # apply a little bit of dilation (closing) to fill up gaps between selected segments, if any
    sws_seg = SNM.binary_closing(sws_seg)
    if not N.any(sws_seg > 0):      return None
    '''

    mid_co = N.round( (N.array(sws_seg.shape)-1) / 2.0).astype(N.int)                     # IMPORTANT: following python convension, in index starts from 0 to size-1 !!!  So (siz-1)/2 is real symmetry center of the volume
    middle_label = sws_seg[mid_co[0], mid_co[1], mid_co[2]]     # get the label of the phi-segment that coveres the middle location of the subtomogram
    if middle_label <= 0:       return None

    '''
    ## exclude such connected regions as candidates of target complex
    dist = CV.grid_distance_to_center( CV.grid_displacement_to_center(size=v.shape) )
    sphere_msk =        dist >= (min(v.shape) / 2.0)

    selected_lbls = set()
    for lbl in range(1, conn_rgn_lbl.max()+1):
        if ((conn_rgn_lbl == lbl) & sphere_msk).sum() == 0:
            selected_lbls.add(lbl)

    if len(selected_lbls) == 0:     return None         # this means that no target complex region is found, all structural regions are excluded

    if middle_label not in selected_lbls:       return None     # if the corresponding connected region touches boundary of max sphere, quit
    '''

    sws_seg[sws_seg != middle_label] = 0            # we only keep the phi-segment of middle_label

    v_target = N.array(v_org, order='F')
    v_target[sws_seg == 0] = v_org[sws_seg > 0].mean()


    # record resulting subtomogram and segment etc...
    out_path = os.path.join(op['out_dir'], record['subtomogram'][len(op['common_path']):])

    out_path_dir = os.path.dirname(out_path)
    if not os.path.isdir(out_path_dir):
        try:
            os.makedirs(out_path_dir)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(out_path_dir):
                pass
            else: raise
    
    IF.put_mrc(v_target, out_path)

    if op['debug']:
        IF.put_mrc(vg, out_path+'-gau.mrc')
        IF.put_mrc(phi, out_path+'-ac-phi.mrc')
        IF.put_mrc(N.array(phi > 0, dtype=N.float, order='F'), out_path+'-ac-seg.mrc')
        IF.put_mrc(N.array(conn_rgn_lbl, dtype=N.float, order='F'), out_path+'-conn.mrc')
        IF.put_mrc(N.array(sws_seg, dtype=N.float, order='F'), out_path+'-ws-seg.mrc')

    record['subtomogram'] = out_path

    return record




if __name__ == '__main__':

    op_file = sys.argv[1]


    with open(op_file) as f:    op = json.load(f)

    # we use the same smooth as the one used in particle picking
    with open(op['smooth']['dog_config_file']) as f:        dog_op = json.load(f)
    op['smooth'] = dog_op[op['smooth']['param_field']]



    with open(op['data_config_in']) as f:           data = json.load(f)


    op['common_path'] = os.path.commonprefix([ _['subtomogram'] for _ in data ])

    # parallel processing
    pool = Pool()
    pool_apply = []
    for i,r in enumerate(data):
        r['sequential_id'] = i
        #neighbor_exclude(record=r, op=op)
        pool_apply.append(  pool.apply_async( func=neighbor_exclude, kwds={'record':r, 'op':op} )  )

    data_new = []
    for pa in pool_apply:
        r = pa.get()
        if r is not None:
            data_new.append(r)
            print '\r', r['sequential_id'], '      ',

    with open(op['data_config_out'], 'w') as f:           json.dump(data_new, f, indent=2)

