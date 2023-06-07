
'''
generate pose normalized subtomgorams through following steps
first approximate subtomograms as a level set, then perform pose normalization on only the positive part of the level set

the functions here are used to support

~/ln/frequent_structure/data/out/method/pose-norm-level-set/scripts/pose_norm_level_set.py




~/ln/tomominer/tomominer/pose/normalize/segmentation/batch/run_json_1.py
'''



import os, sys, json, copy
from multiprocessing.pool import Pool
import numpy as N

import tomominer.pose.normalize.util as PNU
import tomominer.geometry.rotate as GR
import tomominer.io.file as IF
import tomominer.filter.gaussian as FG
import tomominer.segmentation.active_contour.chan_vese.segment as SA

import traceback




def level_set(record, op):
    v_org = IF.get_mrc(record['subtomogram'])
    
    v = v_org - v_org.mean()

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

    return {'phi':phi, 'v_org': v_org, 'v': v, 'vg':vg}



def normalize(record, op):
    if os.path.isfile(record['pose']['subtomogram']):   return {'record':record}

    ls = level_set(record=record, op=op['segmentation'])
    if ls is None:      return

    phi = N.zeros(ls['phi'].shape)
    phi[ls['phi'] > 0] = ls['phi'][ls['phi'] > 0]


    # calculate center of mass, and pca
    c = PNU.center_mass(phi)

    mid_co = (N.array(phi.shape)-1) / 2.0        # center of volume          # IMPORTANT: following python convension, in index starts from 0 to size-1!!! So (siz-1)/2 is real symmetry center of the volume
    if N.sqrt(N.square(c - mid_co).sum()) > (N.min(phi.shape) * op['center_mass_max_displacement_proportion']):    return        # fail if the center of mass is too far away to the subtomogram center, in such case pose normalization may push the structure outside the subtomogram

    rm = PNU.pca(v=phi, c=c)['v']
    loc_r__pn = rm.T.dot(mid_co - c)

    record['pose']['c'] = c.tolist()
    record['pose']['loc_r'] = loc_r__pn.tolist()
    record['pose']['rm'] = rm.tolist()

    if True:
        phi_pn =  GR.rotate(phi, rm=rm, c1=c, default_val=0)
        v_org_pn = GR.rotate_pad_mean(ls['v_org'], rm=rm, c1=c)        # generate segmented and pose normalized subtomogram
    else:
        # perform pose normalization according to loc_r
        phi_pn =  GR.rotate(phi, rm=rm, loc_r=loc_r__pn, default_val=0)
        v_org_pn = GR.rotate_pad_mean(ls['v_org'], rm=rm, loc_r=loc_r__pn)

    return  {'ls':ls, 'phi':phi, 'phi_pn':phi_pn, 'v_org_pn':v_org_pn, 'record':record}

def do_normalize(record, op):

    try:
        r = normalize(record, op)
    except :
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))

    return r


def save_data(r, op):

    record = r['record']

    if 'v_org_pn' not in r:  return record      # this indicates record has already been processed

    if os.path.isfile(record['pose']['subtomogram']):
        # ignore those already segmented subtomograms
        return record

    out_path_dir = os.path.dirname(record['pose']['subtomogram'])
    if not os.path.isdir(out_path_dir):
        try:
            os.makedirs(out_path_dir)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(out_path_dir):
                pass
            else: raise

    #IF.put_mrc(r['phi_pn'], record['pose']['subtomogram']+'-phi-pn.mrc')       # for debugging
    IF.put_mrc(r['v_org_pn'], record['pose']['subtomogram'])

    print 'generated', record['pose']['subtomogram'],
    print

    return record


def main(op, pool=None, n_chunk=1000):

    op['data_config out'] = os.path.abspath(op['data_config out'])
    if not os.path.isdir(os.path.dirname(op['data_config out'])):       os.makedirs(os.path.dirname(op['data_config out']))

    with open(op['data_config in']) as f:           data = json.load(f)

    # convert relative path to absolute path, if needed
    for d in data:
        if not os.path.isabs(d['subtomogram']):     d['subtomogram'] = os.path.abspath(os.path.join(os.path.dirname(op['data_config in']), d['subtomogram']))
        if not os.path.isabs(d['mask']):     d['mask'] = os.path.abspath(os.path.join(os.path.dirname(op['data_config in']), d['mask']))


    op['common_path'] = os.path.commonprefix([ _['subtomogram'] for _ in data ])

    if os.path.isfile(op['data_config out']):
        with open(op['data_config out']) as f:           data_new = json.load(f)
        subtomograms_processed = set([_['subtomogram'] for _ in data_new])
        print 'loaded', len(subtomograms_processed), 'processed subtomograms'
    else:
        data_new = []
        subtomograms_processed = set()

    if 'multiprocessing' not in op:     op['multiprocessing'] = False

    if op['multiprocessing']:
        # parallel processing
        if pool is None:    pool = Pool()

        pool_apply = []
        for i,r in enumerate(data):
            if r['subtomogram'] in subtomograms_processed:      continue

            out_path_root = os.path.join(os.path.abspath(op['out dir']), r['subtomogram'][len(op['common_path']):])
            out_path_root = os.path.splitext(out_path_root)[0]

            assert 'segmentation' not in r
            r['pose'] = {}
            r['pose']['subtomogram'] = out_path_root + '-seg-pn.mrc'
            r['i'] = i

            #normalize(record=r, op=op)
            pool_apply.append(  pool.apply_async( func=do_normalize, kwds={'record':r, 'op':op} )  )

        for pa in pool_apply:
            r = pa.get(99999999)
            if r is None:   continue
            
            print '\rprocessing subtomogram', r['record']['i'], '   ', 'successfully processed', len(data_new)
            sys.stdout.flush()
            del r['record']['i']

            if False:
                rs = save_data(r, op)
                data_new.append(rs)
            else:
                # do not save subtomograms, to reduce storage consumption
                # this saves rotation parameters, center of mass for each subtomogram
                data_new.append(r['record'])

            if (len(data_new) % n_chunk) == 0:
                # keep saving data, in case the process terminates for some problem....
                with open(op['data_config out'], 'w') as f:           json.dump(data_new, f, indent=2)


        del pool
        del pool_apply

    else:
        for i,r in enumerate(data):
            if r['subtomogram'] in subtomograms_processed:      continue

            out_path_root = os.path.join(os.path.abspath(op['out dir']), r['subtomogram'][len(op['common_path']):])
            out_path_root = os.path.splitext(out_path_root)[0]
            assert 'segmentation' not in r
            r['pose'] = {}
            r['pose']['subtomogram'] = out_path_root + '-seg-pn.mrc'
            nr = normalize(record=r, op=op)
            if nr is None:  continue
            data_new.append(nr['record'])
            print '\rprocessing subtomogram', i, '                    ', 'successfully processed', len(data_new), '                             ',         ;           sys.stdout.flush()


            if (len(data_new) % n_chunk) == 0:
                # keep saving data, in case the process terminates for some problem....
                with open(op['data_config out'], 'w') as f:           json.dump(data_new, f, indent=2)


    print 'successfully pose normalized subtomograms in total:', len(data_new)

    with open(op['data_config out'], 'w') as f:           json.dump(data_new, f, indent=2)




'''
modified from 
~/ln/tomominer/tomominer/pose/normalize/segmentation/batch/run_json.py
'''

