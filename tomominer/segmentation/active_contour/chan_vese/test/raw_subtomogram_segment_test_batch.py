#!/usr/bin/env python


# test what you can get by directly applying chan-vese segmentation on raw subtomograms



if __name__ == '__main__':

    import json
    with open('raw_subtomogram_segment_test_batch__op.json') as f:   op = json.load(f)

    import os
    out_dir = os.path.abspath(op['out_dir'])
    if not os.path.isdir(out_dir):      os.makedirs(out_dir)

    with open(op['data_json_file']) as f:   dj = json.load(f)


    import numpy as N
    import tomominer.io.file as IF
    import tomominer.segmentation.active_contour.chan_vese.segment as SACS
    import tomominer.filter.gaussian as FG

    for i, d in enumerate(dj):

        v = IF.read_mrc(d['subtomogram'])['value']

        if not op['density_positive']:      v = -v

        v_phi = SACS.segment_with_poseprocessing(v, op['ac_segmentation'])
        if v_phi is None:   continue
        IF.put_mrc(v_phi, os.path.join(out_dir, '%05d-2-phi.mrc'%(i,)), overwrite=True)

        IF.put_mrc((v_phi>0), os.path.join(out_dir, '%05d-3-phi-seg.mrc'%(i,)), overwrite=True)

        vg = FG.smooth(v, sigma=op['smooth_sigma'])
        IF.put_mrc(vg, os.path.join(out_dir, '%05d-0-vg.mrc'%(i,)), overwrite=True)

        vg_phi = segment(vg, op['ac_segmentation'])
        if vg_phi is None:  continue
        IF.put_mrc(vg_phi, os.path.join(out_dir, '%05d-1-vg-phi.mrc'%(i,)), overwrite=True)

        print '\r', i, '                  ',


