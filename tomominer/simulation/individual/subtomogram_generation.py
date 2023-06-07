#!/usr/bin/env python

'''
load density maps, and generate subtomograms of multiple classes, adapted from    DataProcessingImputation.simulate_several_class_subtomograms_commands()

use multiprocessing for parallel generation


todo: first prepare all dj files for ground truth and instances, and save them, 
then load existing dj file, and generate corresponding instances if the instance file does not exist

'''


import os, sys, shutil, json, pickle, copy, uuid, multiprocessing, tempfile
import pymatlab
from multiprocessing.pool import Pool
import numpy as N
import tomominer.image.vol.util as CV
import tomominer.image.vol.wedge.util as IVWU
import tomominer.io.file as IF
import tomominer.io.path_util as IP
import tomominer.simulation.back_projection_reconstruction as BPR
import tomominer.geometry.rotate as GR
import tomominer.geometry.ang_loc as AAL
import tomominer.model.util as MU


def ground_truth_dj_prepare(cluster_label, out_dir):
    dj = []
    for pdb_id in cluster_label.keys():
        fn = os.path.join(out_dir, '%d.mrc'%(cluster_label[pdb_id]))
        dj.append( {'cluster_label':cluster_label[pdb_id], 'pdb_id':pdb_id, 'ground_truth':fn} )
    return dj


def ground_truth_generate(maps, dj, op):
    op = copy.deepcopy(op)
    op['model']['missing_wedge_angle'] = 0
    op['model']['SNR'] = float('nan')

    session = None

    for d in dj:
        fn = d['ground_truth']
        pdb_id = d['pdb_id']

        if os.path.isfile(fn):      continue

        if session is None:
            session = pymatlab.session_factory(options='-nodesktop -nodisplay')
        mb = BPR.do_reconstruction(session, maps[pdb_id], op)
        if not os.path.isdir(os.path.dirname(fn)):      os.makedirs(os.path.dirname(fn))
        IF.put_mrc(mb, fn)


def subtomogram_dj_prepare(instance_num, shape, loc_proportion, out_dir, max_file_num_per_dir=100):
    loc_max = N.array(shape, dtype=float) * loc_proportion

    dj = []
    for i in xrange(instance_num):
        # form a path name that seperates the files using digits
        f_dir = os.path.join(out_dir, IP.id_digit_seperation(i=i, max_item_per_dir=max_file_num_per_dir))

        fn = os.path.join(f_dir, '%d.mrc'%(i))

        angle = AAL.random_rotation_angle_zyz()
        loc_r = (N.random.random(3)-0.5)*loc_max

        dj.append( {'id':i, 'uuid':str(uuid.uuid4()), 'subtomogram':fn, 'model':{'angle':angle.tolist(), 'loc':loc_r.tolist(), 'bp_rand_seed':N.random.randint(N.iinfo(N.int32).max)} } )

    return dj


def subtomogram_generate(maps, dj, op):

    session = pymatlab.session_factory(options='-nodesktop -nodisplay')
    session.run('clear all;')
    session.run(BPR.param_prepare(op))

    for i, d in enumerate(dj):
        # generate and write the subtomograms
        angle = N.array(d['model']['angle'])
        loc_r = N.array(d['model']['loc'])

        vr = GR.rotate(maps[d['cluster_label']], angle=angle, loc_r=loc_r, default_val=0.0)

        f, fn = tempfile.mkstemp();     os.close(f)
        IF.put_mrc(vr, fn, overwrite=True)
        session.putvalue('vr_fn', fn)
        session.run('im = tom_mrcread(vr_fn);')
        session.run('m = im.Value;')
        os.remove(fn)

        session.run('RandStream.setGlobalStream(RandStream(\'mt19937ar\',\'seed\', %d));'%(d['model']['bp_rand_seed']))
        session.run('mbp = GenerateSimulationMap.backprojection_reconstruction(reconstruction_param, m, reconstruction_param.model.SNR);')
        if ('inverse_intensity' in op) and op['inverse_intensity']:        session.run('mbp = -mbp')

        session.putvalue('bp_fn', fn)
        session.run('tom_mrcwrite(mbp, \'name\', bp_fn);')

        vr_bp = IF.read_mrc_vol(fn)
        if N.isfinite(vr_bp).sum() != vr_bp.size:
            print 'bp reconstruction fail for', d['subtomogram']            ;            sys.stdout.flush()
            os.remove(fn)
            continue

        f_dir = os.path.dirname(d['subtomogram'])
        if not os.path.isdir(f_dir):
            try:
                os.makedirs(f_dir)
            except:
                pass

        assert  not os.path.exists(d['subtomogram'])
        shutil.copyfile(fn, d['subtomogram'])
        os.remove(fn)

        print '\r',     float(i+1) / len(dj),         ;               sys.stdout.flush()

    #return {'maps': maps}

def all_generate(op):
    spacing = op['spacing']
    resolution = op['resolution']

    out_dir = os.path.join(op['out_dir'], str(op['backprojection_reconstruction']['model']['missing_wedge_angle']), str(op['backprojection_reconstruction']['model']['SNR']))
    print 'out_dir:', out_dir
    if not os.path.isdir(out_dir):            os.makedirs(out_dir)
    
    # Prepare Ground Truth
    gt_dj_file = os.path.join(out_dir, 'cluster_info.json')
    if not os.path.isfile(gt_dj_file):
        cluster_label = {pdb_id:l for l, pdb_id in enumerate(op['pdb_ids'])}
        gt_dj = ground_truth_dj_prepare(cluster_label=cluster_label, out_dir=os.path.join(out_dir, 'ground_truth'))
        with open( gt_dj_file, 'w' ) as f:        json.dump(gt_dj, f, indent=2)
        print 'saved', gt_dj_file
    else:
        with open(gt_dj_file) as f:     gt_dj = json.load(f)
        cluster_label = {_['pdb_id']:_['cluster_label'] for _ in gt_dj}
        print 'loaded', gt_dj_file

    del gt_dj_file

    print '# load maps'
    with open(op['map_file']) as f:     maps_org = pickle.load(f)

    maps = {}
    for l, pdb_id in enumerate(op['pdb_ids']):
        if op['density_positive']:
            maps[pdb_id] = maps_org[pdb_id][spacing][resolution]['map']
        else:
            maps[pdb_id] = -maps_org[pdb_id][spacing][resolution]['map']

    print '# resize the maps to same size'
    maps = CV.resize_center_batch_dict(maps, enlarge_factor=op['enlarge_factor'], cval=0.0)
    size = N.array(maps[pdb_id].shape)
    print 'map size', size

    print '# generate noise free ground truth'
    ground_truth_generate(maps, dj=gt_dj, op=op['backprojection_reconstruction'])

    print '# genetate and save wedge mask'
    mask_file = os.path.join(out_dir, 'wedge-mask.mrc')
    mask = IVWU.wedge_mask(size, op['backprojection_reconstruction']['model']['missing_wedge_angle']) * MU.sphere_mask(size);     assert      N.all(N.isfinite(mask))
    IF.put_mrc(mask, mask_file)

    # load exising data_json file, if any, in order to keep the previously generated information
    dj_file = os.path.join(out_dir, 'data_config.json')
    if not os.path.isfile(dj_file):
        dj = []

        for l, pdb_id in enumerate(op['pdb_ids']):
            djt = subtomogram_dj_prepare(instance_num=op['instance_num'], shape=maps[pdb_id].shape, loc_proportion=op['loc_proportion'], out_dir=os.path.join(out_dir, 'subtomograms', str(l)))

            for d in djt:
                d['mask'] = mask_file
                d['cluster_label'] = l

            dj.extend(djt)
            del djt

        with open(dj_file, 'w') as f:       json.dump(dj, f, indent=2)
        print 'saved', dj_file
    else:
        with open(dj_file) as f:       dj = json.load(f)
        print 'loaded', dj_file


    # we only generate subtomograms that do not exist previously
    dj_new = []
    for d in dj:
        if os.path.isfile(d['subtomogram']):    continue
        dj_new.append(d)

    dj = dj_new

    print '# parallel generate instances', len(dj)
    if len(dj) == 0:        return

    pool = Pool(processes=op['worker_num'])
    pool_inds = list(range(len(dj)))
    chunk_size = int(N.ceil(len(dj) / multiprocessing.cpu_count()))
    chunk_size = N.max((chunk_size, 2))

    maps_label = {cluster_label[_]:maps[_] for _ in maps}       # just re-indexing density maps using cluster labels instead of pdb_ids, for the convinence of subtomogram generation
    
    pool_results = []
    while pool_inds:
        pool_inds_t = pool_inds[:chunk_size]

        bp_op = copy.deepcopy(op['backprojection_reconstruction'])
        #subtomogram_generate(maps=maps_label, dj=[dj[_] for _ in pool_inds_t], op=bp_op)
        paar = pool.apply_async(   func=subtomogram_generate, kwds={'maps':maps_label, 'dj':[dj[_] for _ in pool_inds_t], 'op':bp_op}    )
        pool_results.append(paar)
        pool_inds = pool_inds[chunk_size:]

    for r in pool_results:
        r.get(999999)

    print "Done"


def main():
    op_file = 'subtomogram_generation__op.json'
    with open(op_file) as f:    op = json.load(f)
    op['out_dir'] = os.path.abspath(op['out_dir'])
    if not os.path.isdir(op['out_dir']):        os.mkdir(op['out_dir'])
    shutil.copyfile(op_file, os.path.join(op['out_dir'], op_file))

    if 'worker_num_max' in op:
        op['worker_num'] = min(op['worker_num_max'], multiprocessing.cpu_count())
    else:
        op['worker_num'] = multiprocessing.cpu_count()

    for snr in op['backprojection_reconstruction']['model']['SNR_s']:
        for wedge in op['backprojection_reconstruction']['model']['missing_wedge_angle_s']:
            op_t = copy.deepcopy(op)
            del op_t['backprojection_reconstruction']['model']['SNR_s']
            del op_t['backprojection_reconstruction']['model']['missing_wedge_angle_s']
            op_t['backprojection_reconstruction']['model']['SNR'] = snr
            op_t['backprojection_reconstruction']['model']['missing_wedge_angle'] = wedge
            all_generate(op=op_t)


if __name__ == "__main__":
    main()
    

