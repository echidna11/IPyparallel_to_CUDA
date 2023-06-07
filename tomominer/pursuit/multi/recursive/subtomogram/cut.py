#!/usr/bin/env python


'''

Given original subtomograms, and the corresponding cluster average, first segment cluster average. And decide minimum bounding sphere that can contain the segment. Then for each subtomogram, rotate the minimum bounding sphere according to the subtomogram's alignment. Then cut out subvolume inside subtomogram. In the same time, need to prepare a missing wedge mask (make sure the FFT center of the new mask coenside with the original one!!!).


~/ln/tomominer/tomominer/pursuit/multi/recursive/subtomogram/cut.py

'''


import os
import sys
import json
import copy
import cPickle as pickle
import numpy as N
import pymatlab

import tomominer.io.file as IF
import tomominer.simulation.tomogram.template_bounding_sphere as STT
import tomominer.geometry.rotate as GR
import tomominer.geometry.point_cloud.util as GPU
import tomominer.geometry.ang_loc as GA
import tomominer.image.vol.util as IVU
import tomominer.model.util as MU
import tomominer.pursuit.multi.util as PMU
import tomominer.segmentation.active_contour.chan_vese.segment as SACS
import tomominer.segmentation.util as SU



def process_one_template(t, op, session, dj):
    v = IF.read_mrc(t['subtomogram'])['value']
    v_org = N.copy(v)
    if not op['density_positive']:      v = -v

    phi = SACS.segment_with_postprocessing(v, op['template segmentation'])
    if phi is None:     return

    s = phi > (phi.max() * op['template segmentation']['phi_propotion_cutoff'])

    # in the experimental data, the contrast may be low, and there could be more than one connected components, in such case, we select the largest
    cc = SU.connected_components(s)
    assert  cc['num'] > 0


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
        s =     cc['label'] == target_lbl



    vt = N.copy(v_org)
    vt[N.logical_not(s)] = float('nan')
    template_file = os.path.abspath(os.path.join(op['out dir'], 'clus-vol-avg-seg.mrc'))
    if not os.path.isdir(os.path.dirname(template_file)):       os.makedirs(os.path.dirname(template_file))
    IF.put_mrc(vt, template_file, overwrite=True)


    s_t = N.ones(s.shape)
    s_t[N.logical_not(s)] = float('nan')


    # ------------------------------------------
    # calculate bounding sphere of s
    b = STT.get_bounding_sphere(s_t, session=session)
    c = N.reshape(b['c'], (1,3))
    rad = N.ceil(b['r'])
    rad = N.min(    (N.min(v.shape) / 2, rad)       )        # prevent the radius to be larger than the original radius
    print 'bounding sphere radius', rad

    # ------------------------------------------
    # assert that all items in dj has the same subtomogram mask file, generate new mask file, make sure that the FFT center coinside
    mask_file_org = set(_['mask'] for _ in dj)
    assert      len(mask_file_org) == 1
    mask_file_org = list(mask_file_org)[0]

    mask_org = IF.read_mrc(mask_file_org)['value']
    mask_org_c = IVU.fft_mid_co(mask_org.shape)
    mask_org_c_val = mask_org[mask_org_c[0], mask_org_c[1], mask_org_c[2]]
    mask_org[mask_org_c[0], mask_org_c[1], mask_org_c[2]] = float('nan')            # use this marker to make sure the cutted mask has the same FFT origion

    mask = mask_org[(mask_org_c[0]-rad):(mask_org_c[0]+rad), (mask_org_c[1]-rad):(mask_org_c[1]+rad), (mask_org_c[2]-rad):(mask_org_c[2]+rad)]
    mask_c = IVU.fft_mid_co(mask.shape)
    assert      N.isnan(mask[mask_c[0], mask_c[1], mask_c[2]])
    mask[mask_c[0], mask_c[1], mask_c[2]] = mask_org_c_val
    mask[MU.sphere_mask(shape=mask.shape, center=IVU.fft_mid_co(mask.shape)) == 0] = 0       # further adding a spherical mask

    mask_file = os.path.abspath(os.path.join(op['out dir'], 'wedge-mask.mrc'))
    if not os.path.isdir(os.path.dirname(mask_file)):       os.makedirs(os.path.dirname(mask_file))
    IF.put_mrc(mask, mask_file, overwrite=True)
    

    # ------------------------------------------
    # generate individual cutted subtomograms
    if 'template guided subtomogram segmentation' in op:
        op['template guided subtomogram segmentation']['density_positive'] = op['density_positive']
        assert      'smooth_weight' not in op['template guided subtomogram segmentation']       ;            op['template guided subtomogram segmentation']['smooth_weight'] = op['template segmentation']['smooth_weight']
        assert      'image_weight' not in op['template guided subtomogram segmentation']        ;           op['template guided subtomogram segmentation']['image_weight'] = op['template segmentation']['image_weight']
        
        print 'template guided subtomogram segmentation', op['template guided subtomogram segmentation']


    djn = []
    for djt in dj:
        print '\r', djt['subtomogram'][len(op['path common prefix']):],         ;       sys.stdout.flush()

        v = IF.read_mrc(djt['subtomogram'])['value']

        # rotate the bounding sphere center according to alignment
        ang, loc = GA.reverse_transform_ang_loc(ang=N.array(djt['angle']), loc_r=N.array(djt['loc']))

        if 'template guided subtomogram segmentation' in op:
            # optional step: mask the subtomogram using template guided segmentation in pursuit.multi.util.template_guided_segmentation()
            sr = GR.rotate(s, angle=ang, loc_r=loc)

            v_tgs = PMU.template_guided_segmentation(v=v, m=sr, op=op['template guided subtomogram segmentation'])

            if v_tgs is not None:
                v = v_tgs
                v[N.logical_not(N.isfinite(v))] = v[N.isfinite(v)].mean()

        # get the rotated minimum bounding sphere center
        cr = GPU.rotate_translate(c, c0=N.array(v.shape)/2.0, angle=ang, loc_r=loc)
        cr = cr.flatten()
        cr = N.round(cr)


        # check if the generated subvolume will go outside the original subtomogram, if so, discard
        cr_s = cr - rad
        if N.any(cr_s < 0):   
            print 'outside left'         ;       sys.stdout.flush()
            continue

        cr_e = cr + rad
        if N.any(cr_e >= N.array(v.shape)):   
            print 'outside right'         ;       sys.stdout.flush()
            continue

        vt = v[cr_s[0]:cr_e[0], cr_s[1]:cr_e[1], cr_s[2]:cr_e[2]]

        # save the new subvolume
        subtomogram_file = os.path.abspath(os.path.join(op['out dir'], djt['subtomogram'][len(op['path common prefix']):]))
        if not os.path.isdir(os.path.dirname(subtomogram_file)):       os.makedirs(os.path.dirname(subtomogram_file))
        IF.put_mrc(vt, subtomogram_file, overwrite=True)

        djt = copy.deepcopy(djt)
        djt['original'] = {'subtomogram':djt['subtomogram'], 'mask':djt['mask']}
        djt['subtomogram'] = subtomogram_file
        djt['mask'] = mask_file

        if 'template' in djt:       del djt['template']
        if 'score' in djt:          del djt['score']
        if 'loc' in djt:            del djt['loc']
        if 'angle' in djt:          del djt['angle']

        djn.append(djt)



    with open(os.path.join(op['out dir'], 'data_config.json'), 'w') as f:       json.dump(djn, f, indent=2)
    print 'extracted', len(djn), 'subtomograms'




def main():
    
    with open('pursuit_multi_recursive_subtomogram_cut__op.json') as f:     op = json.load(f)

    if 'template guided subtomogram segmentation' in op:        op['out dir'] = op['out dir'] + '-seg'
    
    if not os.path.isdir(op['out dir']):    os.makedirs(op['out dir'])

    with open(op['data json file']) as f:      dj_full = json.load(f)

    if op['use true class label']:
        true_label = {_['subtomogram']:_['cluster_label'] for _ in dj_full if ('cluster_label' in _)}
    else:
        true_label = None

    with open(os.path.join(op['pass dir'], 'cluster_average_select.pickle'), 'rb') as f:      cas = pickle.load(f)
    stcf = cas['selected_templates_common_frame']
    st = cas['selected_templates']
    tki = cas['tk_info']


    session = pymatlab.session_factory(options='-nodisplay')

    session.run( 'addpath(\'%s\')'%(op['matlab_code_path']) )
    session.run( 'addpath(\'%s\')'%(op['matlab_bounding_code_path']) )


    # find the path common prefix
    djs = []
    for c in st:        djs.extend(tki[st[c]['subtomogram']]['data_json'])

    op['path common prefix'] = os.path.commonprefix([_['subtomogram'] for _ in djs])


    for c in st:
        if ('selected clusters' in op) and (c not in set(op['selected clusters'])):     continue

        print 'processing cluster', c       ;       sys.stdout.flush()

        opt = copy.deepcopy(op)
        opt['out dir'] = os.path.join( op['out dir'], 'clus-%03d'%(c,), 'subtomograms' )

        # add true class label information, if avaliable
        dj = copy.deepcopy(tki[st[c]['subtomogram']]['data_json'])
        for d in dj:
            if 'cluster_label' in d:    del d['cluster_label']
            if (true_label is not None) and (d['subtomogram'] in true_label):      d['cluster_label'] = true_label[d['subtomogram']]

        print len(dj), 'subtomograms'
        process_one_template(t=stcf[c], dj=dj, op=opt, session=session)

    if 'selected clusters' not in op:
        # export list of the rest subtomograms
        out_dir_t = os.path.join( op['out dir'], 'clus-rest', 'subtomograms' )
        if not os.path.isdir(out_dir_t):     os.makedirs(out_dir_t)

        djs_s = set([_['subtomogram'] for _ in djs])
        with open(os.path.join(out_dir_t, 'data_config.json'), 'w') as f:       json.dump([_ for _ in dj_full if _['subtomogram'] not in djs_s], f, indent=2)



if __name__ == '__main__':
    main()

