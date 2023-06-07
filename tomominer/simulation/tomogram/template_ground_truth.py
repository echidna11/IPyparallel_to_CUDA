#!/usr/bin/env python

# use backprojection reconstruction to generate ground truth template

import os
import sys
import json
import pickle
import numpy as N
import tomominer.io.file as IOF
import tomominer.image.vol.util as CV

if __name__ == '__main__':
    with open('template_ground_truth__config.json') as f:    op = json.load(f)
    with open(op['templete_select_file']) as f:    maps = pickle.load(f)
    with open(op['model_tomogram_subtomogram__config__file']) as f:    mts_op = json.load(f)
    with open(mts_op['back_projection_reconstruction_config_file']) as f:    bp_op = json.load(f)
    if not os.path.isdir(op['out_dir']):     os.makedirs(op['out_dir'])

    # get the size of subtomograms
    mask = IOF.get_mrc(op['subtomogram_mask_file'])
    size = N.array(mask.shape);     print 'subtomogram size', size

    bp_op['random_seed'] = 0
    bp_op['model']['missing_wedge_angle'] = 0
    bp_op['model']['SNR'] = float('NaN')
    
    import oct2py
    session = oct2py.Oct2Py() #pymatlab.session_factory(options='-nodesktop -nodisplay')

    import tomominer.simulation.back_projection_reconstruction__octave as BCR
    file_paths = {}
    for pid in maps:
        m = CV.resize_center(v=maps[pid]['map'], s=size, cval=0.0)
        vb = BCR.do_reconstruction(session=session, dm=m, op=bp_op)
        vb = N.array(vb, dtype=float, order='F')

        if not mts_op['density_positive']:  vb = -vb
        file_paths[pid] = os.path.abspath(     os.path.join(op['out_dir'], '%s.mrc'%(pid))     )
        IOF.put_mrc(vb, file_paths[pid], overwrite=True)

    # generate a JSON file that records file pathes
    ci = []
    for pid in file_paths:
        ci_t = {'cluster_label':maps[pid]['id'], 'pdb_id':pid, 'ground_truth':file_paths[pid]}
        ci.append(ci_t)

    with open(os.path.join(op['out_dir'], 'cluster_info.json'), 'w') as f:      json.dump(ci, f, indent=2)

