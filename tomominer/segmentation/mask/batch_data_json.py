#!/usr/bin/env python


'''
assume that each subtomogram only contain a single particle.
given option file and a data json file, for each subtomogram in the json file, perform gauss smoothing, then use chan-vese to segment to identify the region occupied by the target complex. Then mask out the outside region

calculate the max size of the regions occupied by segments, 
save the masked subtomograms according to the diagonal of such max size so that the subtomograms are tight and the particle is located at center of the subtomogram and the rotation of the subtomogram does not move the particle to outside
'''


import os
import sys
import errno
import json
from multiprocessing.pool import Pool
import numpy as N
import scipy.ndimage.morphology as SNM

import tomominer.io.file as IF
import tomominer.filter.gaussian as FG
import tomominer.segmentation.active_contour.chan_vese.segment as SA
import tomominer.segmentation.watershed as SW
import tomominer.image.vol.util as CV
import tomominer.model.util as MU

import tomominer.core as core

def mask(record, op):

    re = {'record':record}

    v_org = IF.get_mrc(record['subtomogram'])
    v_org = N.array(v_org, dtype=N.float)

    
    if ('background' in op) and op['background']['for_smoothing']:
        # use given background values instead of taking means
        # IMPORTANT: inside a crowded cytoplasm, the mean intensity of a subtomogram is often larger than the background intensity, in such case the whole subtomogram tend to stand out, therefore we should not pad useing background in such case. On the other hand, we should use backgound values outside the mask
        baseline = op['background']['value'][record['tomogram_id']]['mean']
    else:
        #print 'v_org', v_org.mean(), v_org.std(),
        baseline= v_org.mean()

    v = v_org - baseline

    if not op['density_positive']:     v = -v

    # ----------------------------------------------------------
    # first perform padding and low pass filtering then perform segmentation
    if 'smooth' in op:
        vg = FG.smooth(v, sigma=op['smooth']['sigma'])
    else:
        vg = v

    phi = N.sign(vg - N.abs(vg).mean())
    phi = SA.segment(vg, phi, op['ac_segmentation']['smooth_weight'], op['ac_segmentation']['image_weight'] / vg.var(), print_progress=op['ac_segmentation']['print_progress'])
    if phi is None:     return
    if (phi > 0).sum() == 0:        return
    if (phi < 0).sum() == 0:        return

    if vg[phi > 0].mean() < vg[phi < 0].mean():    phi = -phi      # assume bright regions in subtomograms corresponding to high electron densities

    re['v_org'] = v_org
    re['vg'] = vg
    re['phi'] = phi

    #--------------------------------------------------------------
    # find all connected components that overlaps with a spherical mask of radius sigma, and keep only the largest component
    if 'connect_region_overlap_radius' in op:
        mcr_overlap_sphere = MU.sphere_mask(phi.shape, radius=op['connect_region_overlap_radius'])
    else:
        mcr_overlap_sphere = MU.sphere_mask(phi.shape)

    m_struc = (phi > 0)     # identified structural regions
    if 'erosion' in op:     m_struc = SNM.binary_erosion(m_struc, iterations=op['erosion']['iterations'])       # remove small connections by erotion
    m_struc = N.array(m_struc, dtype=N.int32, order='F')

    re['m_struc'] = m_struc


 
    mcr = core.connected_regions(  m_struc  )

    mcr_largest_lbl = {'num':0}
    for l in range(1, mcr['max_lbl'] + 1):
        mcr_lbl_num_t = N.logical_and((mcr['lbl'] == l), mcr_overlap_sphere).sum()          # we assume the target complex of interest must have the largest overlap with mcr_overlap_sphere
        if mcr_lbl_num_t > mcr_largest_lbl['num']:
            mcr_largest_lbl['num'] = mcr_lbl_num_t
            mcr_largest_lbl['lbl'] = l

    if 'lbl' not in mcr_largest_lbl:    return None

    re['m_struc_target'] =         (mcr['lbl'] == mcr_largest_lbl['lbl'])           # only the target complex region defined by the largest connected component with phi > 0

    # watershed segmentation on phi to find regions occupied by different connected components
    mcr_sws = SW.segment(vol_map=phi, vol_lbl=mcr['lbl'])

    
    m =         N.logical_and(  (mcr_sws['vol_seg_lbl'] == mcr_largest_lbl['lbl']), (phi > (phi.max()*op['ac_segmentation']['phi_propotion_cutoff']))   )          # since phi is signed distance, the cutoff gives the distance to the segment

    if (m > 0.5).sum() == 0:
        print 'no structural region detected'
        return None

    if op['within_sphere_check']:
        # test if the selected structural region (the largest connected component) is totally inside the max sphere of the subtomogram, if not, ignore this segment. In such case, we assume that the complex of interest is totall inside such sphere and rotation will not clip it. This should be useful for eliminating those membrane segments
        m_c = (mcr['lbl'] == mcr_largest_lbl['lbl'])
        m_sp = MU.sphere_mask(m_c.shape)
        m_c_sp = N.logical_and(m_c, (m_sp==0))

        re['m_c'] = m_c
        re['m_c_sp'] = m_c_sp

        if m_c_sp.sum() > 0:
            print 'target complex too large', m_c_sp.sum()
            return None




    # -------------------------------------------------------------
    if ('subtomogram' in op) and ('size' in op['subtomogram']):
        raise       Exception('Need to save resized mask as well')

        siz = N.array(op['subtomogram']['size'])

        # put the mask and subtomogram into a predefined smaller volume
        mw =    N.where(m)

        bdy = N.zeros( (3,2) )        # boundary positions for each dimension
        for i in range(3):
            bdy[i,0] = mw[i].min()
            bdy[i,1] = mw[i].max()

        bdy_c = bdy.mean(axis=1).flatten()        # center location
        bdy_siz =       (bdy[:,1] - bdy[:,0]).flatten()

        if N.any(bdy_siz >= siz):
            print 'particle too large:', record['subtomogram'], bdy_siz, siz
            return None            # in this case, the specified subtomogram size cannot hold the entire segment


        m = CV.cut_from_whole_map(whole_map=N.array(m, dtype=float), c=bdy_c, siz=siz)
        if m is None:
            print 'Cannot cut out smaller subtomogram'
            return None

        v_org = CV.cut_from_whole_map(whole_map=v_org, c=bdy_c, siz=siz)
        if v_org is None:
            print 'Cannot cut out smaller subtomogram'
            return None

    re['m'] = N.array(m, dtype=N.float, order='F')


    v_target = N.array(v_org, order='F')

    if ('background' in op) and op['background']['for_masking']:
        # substract background intensity within mask, and set zero outside mask, this is useful for calculating PCA type of pose normalizations
        # in such case, we can guarantee to segment structural region against background region.
        v_target__default_val = 0.0
        v_target -= op['background']['value'][record['tomogram_id']]['mean']
        v_target[m <= 0.5] = 0.0
    else:
        v_target__default_val = v_org[m > 0.5].mean()
        v_target[m <= 0.5] = v_target__default_val


    re['v_target'] = v_target
    re['v_target__default_val'] = v_target__default_val

    return re


# save subtomograms and other information returned from mask()
def save_data(r, op):
    record = r['record']
    out_path = os.path.join(os.path.abspath(op['out_dir']), record['subtomogram'][len(op['common_path']):])

    if os.path.isfile(out_path):
        # ignore those already segmented subtomograms
        record['subtomogram'] = out_path
        return record

    out_path_dir = os.path.dirname(out_path)
    if not os.path.isdir(out_path_dir):
        try:
            os.makedirs(out_path_dir)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(out_path_dir):
                pass
            else: raise

    if 'v_target' in r:    IF.put_mrc(r['v_target'], out_path)

    if op['debug']:
        if 'v_org' in r:       IF.put_mrc(r['v_org'], out_path+'-org.mrc')
        if 'vg' in r:          IF.put_mrc(r['vg'], out_path+'-pg.mrc')
        if 'phi' in r:         IF.put_mrc(r['phi'], out_path+'-ac-phi.mrc')

        if 'm_struc' in r:     IF.put_mrc(r['m_struc'], out_path+'-struc.mrc')

        if 'm_c' in r:         IF.put_mrc(m_c, out_path+'-mc.mrc')
        if 'm_c_sp' in r:      IF.put_mrc(m_c_sp, out_path+'-mc-sp.mrc')

        if 'm' in r:           IF.put_mrc(r['m'], out_path+'-ac-seg.mrc')


    record = r['record']
    record['subtomogram'] = out_path
    print 'generated', record['subtomogram'],
    if 'background' in op:        print 'using bg',
    print


    return record



def main(op_file):

    with open(op_file) as f:    op = json.load(f)

    if 'background' in op:
        with open(op['background']['value_file']) as f:        bgv = json.load(f)
        del op['background']['value_file']

        bgv = {_['tomogram_id']:_ for _ in bgv}

        op['background']['value'] = bgv


    with open(op['data_config_in']) as f:           data = json.load(f)


    op['common_path'] = os.path.commonprefix([ _['subtomogram'] for _ in data ])

    if False:
        # just for debugging
        for i,r in enumerate(data):
            r['sequential_id'] = i
            mask(record=r, op=op)
     
    # parallel processing
    pool = Pool()
    pool_apply = []
    for i,r in enumerate(data):
        r['sequential_id'] = i
        #mask(record=r, op=op)
        pool_apply.append(  pool.apply_async( func=mask, kwds={'record':r, 'op':op} )  )

    data_new = []
    for pa in pool_apply:
        r = pa.get(99999999)
        if r is None:   continue

        rs = save_data(r, op)
        data_new.append(rs)
        print '\r', r['sequential_id'], '      ',

    del pool
    del pool_apply

    print 'successfully segmented subtomograms:', len(data_new)

    with open(op['data_config_out'], 'w') as f:           json.dump(data_new, f, indent=2)



if __name__ == '__main__':
    main(sys.argv[1])

