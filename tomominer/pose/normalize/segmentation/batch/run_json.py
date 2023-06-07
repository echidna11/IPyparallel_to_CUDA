#!/usr/bin/env python

import os
import sys
import json
from multiprocessing.pool import Pool
import numpy as N

import tomominer.segmentation.mask.batch_data_json as SMB
import tomominer.pose.normalize.util as PNU
import tomominer.geometry.rotate as GR
import tomominer.io.file as IF



'''
do pose normalization through level set and segmentation: 1) gaussian filtering, 2) chan vese 3) removal of neighbor structures 4) calculate center of mass of phi on phi >0 region, 5) apply PCA to phi only on phi >0 region. 6)  use PCA for pose normalization (translate and rotate subtomograms), then perform k-means clustering to separate subtomograms. The cluster number k can be manually chosen from inspecting the cluster averages, and hierarchical clustering of averages.

~/ln/tomominer/tomominer/pose/normalize/segmentation/batch/run_json.py

'''


def normalize(record, op):
    if os.path.isfile(record['segmentation']['pose']['subtomogram']):   return {'record':record}

    mr = SMB.mask(record=record, op=op['segmentation'])
    if mr is None:      return

    phi = N.zeros(mr['phi'].shape)
    phi[mr['m_struc_target']] = mr['phi'][mr['m_struc_target']]


    # calculate center of mass, and pca
    c = PNU.center_mass(phi)

    mid_co = N.array(phi.shape) / 2
    if N.sqrt(N.square(c - mid_co).sum()) > (N.min(phi.shape) * op['center_mass_max_displacement_proportion']):    return        # fail if the center of mass is too far away to the subtomogram center, in such case pose normalization may push the structure outside the subtomogram

    pca_re = PNU.pca(v=phi, c=c)
    rm = pca_re['v']

    record['segmentation']['pose']['c'] = c.tolist()
    record['segmentation']['pose']['rm'] = rm.tolist()
    record['segmentation']['pose']['w'] = pca_re['w'].tolist()



    phi_pn =  GR.rotate(phi, rm=rm, c1=c, default_val=0)
    v_target_pn = GR.rotate(mr['v_target'], rm=rm, c1=c, default_val=mr['v_target__default_val'])        # generate segmented and pose normalized subtomogram



    return  {'mr':mr, 'phi':phi, 'phi_pn':phi_pn, 'v_target_pn':v_target_pn, 'record':record}


def save_data(r, op):

    record = r['record']

    if 'v_target_pn' not in r:  return record      # this indicates record has already been processed

    if os.path.isfile(record['segmentation']['pose']['subtomogram']):
        # ignore those already segmented subtomograms
        return record

    out_path_dir = os.path.dirname(record['segmentation']['pose']['subtomogram'])
    if not os.path.isdir(out_path_dir):
        try:
            os.makedirs(out_path_dir)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(out_path_dir):
                pass
            else: raise

    #IF.put_mrc(r['phi_pn'], record['segmentation']['pose']['subtomogram']+'-phi-pn.mrc')       # for debugging
    IF.put_mrc(r['mr']['v_target'], record['segmentation']['subtomogram'])        # record the segmented subtomogram, in case if there is a need to use it for classification
    IF.put_mrc(r['v_target_pn'], record['segmentation']['pose']['subtomogram'])

    print 'generated', record['segmentation']['pose']['subtomogram'],
    if 'background' in op['segmentation']:        print 'using bg',
    print

    return record


def main(op, pool=None):


    if 'background' in op['segmentation']:      del op['segmentation']['background']

    with open(op['data_config_in']) as f:           data = json.load(f)

    op['common_path'] = os.path.commonprefix([ _['subtomogram'] for _ in data ])

    
    # parallel processing
    if pool is None:    pool = Pool()

    pool_apply = []
    for i,r in enumerate(data):
        out_path_root = os.path.join(os.path.abspath(op['out_dir']), r['subtomogram'][len(op['common_path']):])
        out_path_root = os.path.splitext(out_path_root)[0]

        assert 'segmentation' not in r
        r['segmentation'] = {}
        r['segmentation']['subtomogram'] = out_path_root + '-seg.mrc'
        r['segmentation']['pose'] = {}
        r['segmentation']['pose']['subtomogram'] = out_path_root + '-seg-pn.mrc'

        #normalize(record=r, op=op)
        pool_apply.append(  pool.apply_async( func=normalize, kwds={'record':r, 'op':op} )  )

    data_new = []
    for pa in pool_apply:
        r = pa.get(99999999)
        if r is None:   continue

        rs = save_data(r, op)
        data_new.append(rs)

    del pool
    del pool_apply

    print 'successfully segmented and pose normalized subtomograms:', len(data_new)

    with open(op['data_config_out'], 'w') as f:           json.dump(data_new, f, indent=2)




if __name__ == '__main__':
    with open(sys.argv[1]) as f:    op = json.load(f)
    with open(op['segmentation_op_file']) as f:           op['segmentation'] = json.load(f)
    del op['segmentation_op_file']

    main(op)

