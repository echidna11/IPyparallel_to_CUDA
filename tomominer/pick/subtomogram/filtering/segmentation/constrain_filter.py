#!/usr/bin/env python


'''

we want to keep only those subtomograms that clearly defines a particle inside

given a subtomogram, 

option 1: first gaussian smooth. Then perform level set segmentation to get structural regions. 
option 2: use anistropic diffusion and gaussian mixture for segmentation, see segmentation/ani_dif_int_clust.py

then find connected components, and perform watershed on phi. Then the watershed segment that covers the center of the subtomogram contains the structure of interest. We require the corresponding connected region to be the largest  connected region. We further require the connected component of interest do not touch the boundary of the subtomogram, or within the max sphere of the subtomogram

the output of the program is a list of filtered subtomograms

'''

import os
import sys
import json
import numpy as N
from multiprocessing.pool import Pool as Pool


import tomominer.io.file as IF
import tomominer.segmentation.watershed as SW
import tomominer.segmentation.active_contour.chan_vese.segment as SACS
import tomominer.segmentation.ani_dif_int_clust as SADIC
import tomominer.filter.gaussian as FG
import tomominer.segmentation.util as SU
import tomominer.model.util as MU

import tomominer.core as core

def check(d, op):
    v = IF.read_mrc(d['subtomogram'])['value']

    if 'debug' in op:       
        subtomogram_file_name = os.path.splitext(os.path.basename(d['subtomogram']))[0]
    else:
        subtomogram_file_name = None


    if not op['density positive']:    v = -v
    if 'gaussian sigma' in op:      v = FG.smooth(v, sigma=op['gaussian sigma'])

    if 'debug' in op:   IF.put_mrc(v, os.path.join(op['debug']['out dir'], '%s-0-v.mrc'%(subtomogram_file_name,)), overwrite=True)

    if 'level set segmentation' in op:
        assert      'anistropic diffusion segmentation' not in op

        # level set segmentation 
        phi = SACS.segment_with_postprocessing(v, op['level set segmentation'])
        if phi is None:     return {'d':d, 'success':False, 'err': 'phi is None'}

        struc = phi > 0             # get structural region

    elif 'anistropic diffusion segmentation' in op:
        assert      'level set segmentation' not in op

        # anistropic diffusion segmentation
        sre = SADIC.do_segmentation(v, op['anistropic diffusion segmentation'])
        struc = (sre['lbl'] == sre['lbl'].max())          # get structural region
        phi = struc.astype(N.float)
        phi -= 0.5

        #phi = core.ac_distance_transform_3d(phi.astype(N.uint8))         # distance transform to convert structural region into a level set
        phi = SACS.ac_reinit(phi)           # distance transform to convert structural region into a level set

        if 'debug' in op:   IF.put_mrc(sre['vf'], os.path.join(op['debug']['out dir'], '%s-1-ani.mrc'%(subtomogram_file_name,)), overwrite=True)
        if 'debug' in op:   IF.put_mrc(sre['vf_mean'], os.path.join(op['debug']['out dir'], '%s-2-ani-mean.mrc'%(subtomogram_file_name,)), overwrite=True)
        if 'debug' in op:   IF.put_mrc(sre['lbl'], os.path.join(op['debug']['out dir'], '%s-3-lbl.mrc'%(subtomogram_file_name,)), overwrite=True)


    if 'debug' in op:   IF.put_mrc(phi, os.path.join(op['debug']['out dir'], '%s-7-phi.mrc'%(subtomogram_file_name,)), overwrite=True)
    if 'debug' in op:   IF.put_mrc(struc, os.path.join(op['debug']['out dir'], '%s-8-struc.mrc'%(subtomogram_file_name,)), overwrite=True)



    # get connected components
    cc = SU.connected_components(struc)
    if 'debug' in op:   IF.put_mrc(cc['label'], os.path.join(op['debug']['out dir'], '%s-5-cc.mrc'%(subtomogram_file_name,)), overwrite=True)

    if cc['num'] > 1:

        # find the component with largest size
        target_lbl = None
        target_structure_size = -1

        for l in range(1, cc['num']+1):
            l_size = (cc['label'] == l).sum()
            if l_size <= target_structure_size:     continue
            target_structure_size = l_size
            target_lbl = l

        assert      target_lbl is not None
        if ("minimum volume proportion" in op) and (target_structure_size < (float(v.size) * op['minimum volume proportion'])):
            #print target_structure_size
            return {'d':d, 'success':False, 'err': 'target_structure_size < minimum volume proportion'}         # use this constrain to require the minimum size of the target complex

        target_structure =      cc['label'] == target_lbl

        if op['center occupancy check']:
        # check if the largest component occupies center region (after watershed expansion). Note that such constrain may filter out the case that a large sphere contains a small ball inside. However, in this case the small ball should be regareded as an isolated complex?
            sws = SW.segment(vol_map=phi, vol_lbl=cc['label'])

            # get the watershed label at center
            mid_co = N.round(N.array(v.shape) / 2)
            if sws['vol_seg_lbl'][mid_co[0], mid_co[1], mid_co[2]] != target_lbl :     return {'d':d, 'success':False}


    elif cc['num'] == 1:
        target_lbl = 1
        target_structure =      cc['label'] == target_lbl

    else:
        return {'d':d, 'success':False, 'err': 'cc_num == 0'}

    if 'debug' in op:   IF.put_mrc(target_structure, os.path.join(op['debug']['out dir'], '%s-6-ts.mrc'%(subtomogram_file_name,)), overwrite=True)


    if op['within max sphere check']:
        # check if target_structure is within the max sphere of the subtomogram
        sm = MU.sphere_mask(shape=v.shape, radius=((N.min(v.shape)/2.0)-1))
        if N.logical_and(target_structure, (sm == 0)).sum() > 0:        return {'d':d, 'success':False, 'err': 'outside max sphere'}

    if op['notouch subtomogram boundary check']:
        # check if target_structure touches the boundary of the subtomogram
        bm = N.zeros(target_structure.shape, dtype=N.int)        # mask of boundary region
        bm[0,:,:] = 1       ;       bm[-1,:,:] = 1       ;       bm[:,0,:] = 1       ;       bm[:,-1,:] = 1       ;       bm[:,:,0] = 1       ;       bm[:,:,-1] = 1
        if N.logical_and(target_structure, (bm == 1)).sum() > 0:        return {'d':d, 'success':False, 'err': 'touches subtomogram boundary'}



    return {'d':d, 'success':True, 'err': ''}




def main():
    with open('pick_subtomogram_filtering_segmentation_constrain_filter__op.json') as f:        op = json.load(f)

    with open(op['input data json file']) as f:     dj = json.load(f)



    pool = Pool()
    pre = []
    for d in dj:
        #c = check(d=d, op=op['check'])
        pre.append( pool.apply_async(func=check, kwds={'d':d, 'op':op['check']}) )

    common_prefix = os.path.commonprefix([_['subtomogram'] for _ in dj])

    djt = []
    for i, r in enumerate(pre):
        c = r.get(999999)
        print '\r', len(djt), c['d']['subtomogram'][len(common_prefix):], c['err'], '                            ',            ;              sys.stdout.flush()
        if not c['success']:    continue
        djt.append(c['d'])

    with open(op['output data json file'], 'w') as f:       json.dump(djt, f, indent=2)



if __name__ == '__main__':
    main()



