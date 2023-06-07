#!/usr/bin/env python


'''
generate crowded tomograms containing many hetrogenious complexes, 
    0) select templates                             (template_select.py)
    1) segment templates, then estimate radius      (template_bounding_sphere.py)
    1) generate models through sphere packing       (model_generation.py)
    2) generate density maps according to models
    3) generate tomograms 

then perform peak detection and extract subtomograms, 
record ground truth lables of the extracted subtomograms

in addition, generate ground truth templates (template_ground_truth.py)

perform the above in parallel and obtain large amount of subtomograms

to prepare classify config file, use json/list/concatenate.py to concat all info.json file.  For example ~/ln/tomominer/tomominer/json/list/concatenate.py `find ../../subtomograms/ -name 'info.json'` > data_config.json

'''

import os
import sys
import pickle
import json
import shutil
import copy
import numpy as N
import tomominer.io.file as IOF

# function below is pipeline for processing a single model (particluar SNR of particular model number), i.e.:
# 1. generate density map
# 2. perform backprojection reconstruction
# 3. detect peaks
# 4. extract and save subtomograms
# 5. record the subtomogram's true label if any

def model_processing(op):

    assert      op['model_id'] is not None
    with open(op['templete_select_file']) as f:     maps = pickle.load(f)
    with open(op['model_generation_config_file']) as f:     model_op = json.load(f)
    with open(op['model_file']) as f:       models = json.load(f)
    with open(op['bounding_sphere_file']) as f:     bounding_spheres = pickle.load(f)
    bounding_spheres = bounding_spheres['bounding_spheres']
    # total tomogram volume parameters: box size
    box  = model_op['packing']['param']['box']
    map_size = N.array( [box['x'], box['y'], box['z']] )
    # out directory is structured like this: /subtomograms/<totoal number of particles>/<x, y, z parameters of tomogram volume>
    # CREATING DIRECTORIES
    # processing directory
    processing_dir_root = os.path.join(op['out_dir'], str(model_op['copy_number']['total']), '%d-%d-%d'%(model_op['packing']['param']['box']['x'], model_op['packing']['param']['box']['y'], model_op['packing']['param']['box']['z']))
    
    # Output directory of extracted subtomograms
    # the extension of out dir starts with whether peak detection (p) or ground truth (g) location is used
    out_dir = os.path.abspath( os.path.join(processing_dir_root, op['model_id'], str(op['bp']['model']['missing_wedge_angle']), str(op['bp']['model']['SNR']), ('peak-dog' if 'peak' in op else 'peak-true'), str(op['subtomogram']['size_ratio'] if 'size_ratio' in op['subtomogram'] else op['subtomogram']['size']), 'subtomograms'))
    
    # Tomogram directory (With noise)
    tomogram_dir = os.path.join(processing_dir_root, op['model_id'], str(op['bp']['model']['missing_wedge_angle']), str(op['bp']['model']['SNR']), 'tomogram')
    if not os.path.exists(tomogram_dir):      os.makedirs(tomogram_dir)
    # Tomorgram file
    tomogram_file = os.path.join(tomogram_dir, '%s-bp.npy'%(op['model_id'],))
    
    # density map directory (Original)
    density_map_dir = os.path.join(processing_dir_root, 'density-maps')
    if not os.path.exists(density_map_dir):      os.makedirs(density_map_dir)
    # density map file
    density_map_file = os.path.join(density_map_dir, '%s.npy'%(op['model_id']))
    
    #To skip loading of tomorgams and debsity files, directly check whether info file is available, because it is generated at the end of subtomogram extraction
    if op['regenerate'] is "False":
        if 'peak' in op:
            if os.path.isfile(os.path.join(out_dir, "info.json")):
                print "Particle picking and subtomogram extraction already done for model ", op['model_id'], " SNR = ", str(op['bp']['model']['SNR'])
                print "If you want to regenerate particles and subtomograms, make regenerate variable in program to be True"
                return
        elif os.path.isfile(os.path.join(out_dir, "info.json")):
            print "Ground truth avialable for model ", op['model_id'], " SNR = ", str(op['bp']['model']['SNR'])
            return
                
    # If we want to make sure that subtomograms are RE-GENERATED,make regenerate = "True" in json file
    if os.path.isdir(out_dir):     shutil.rmtree(out_dir)
    os.makedirs(out_dir)
        
    # If density file is not present then create one. We need not create the map file again and again for different SNR values. So use existing one.
    if not os.path.isfile(density_map_file):
        print '# generate density map according to model'
        import tomominer.simulation.tomogram.model_to_density_map as MOD
        v = MOD.get_map(size=map_size, model=models[op['model_id']]['instances'], bounding_spheres=bounding_spheres, maps=maps)

        if op['density_positive'] == False:     v = -v

        # Write we save density map in npy (numpy) format, and convert to mrc format later (to visualize and apply template matching).
        # The program used to convert .npy to .mrc need to be editted so that it also includes header information while writing mrc file.

        #IOF.put_mrc(N.array(v, dtype=float, order='F'), density_map_file, overwrite=True)
        N.save(density_map_file, v)
    else:
        print 'loading', density_map_file
        v = N.load(density_map_file)

    # If tomogram file is not available then create one. This one is with Noise.
    if not os.path.isfile(tomogram_file):
        print '# construct simulation tomograms'
        import oct2py
        session = oct2py.Oct2Py()

        import tomominer.simulation.back_projection_reconstruction__octave as BCR
        op['bp']['random_seed'] = models[op['model_id']]['random_seed']
        vt = BCR.do_reconstruction(session=session, dm=v, op=op['bp'], verbose=True)

        del session     # close matlab/octave to save storage

        #IOF.put_mrc( N.array(vt, dtype=float, order='F'), tomogram_file, overwrite=True )
        # File saved in .npy format
        N.save(tomogram_file, vt)
    else:
        print 'loading', tomogram_file
        vt = N.load(tomogram_file)

 
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    # make a copy of all configuration and model files, for backup purpose
    config_dir = os.path.join(out_dir, 'config')
    if not os.path.isdir(config_dir):
        os.makedirs(config_dir)
        for f in os.listdir('./'):
            if not os.path.isfile(f):       continue
            if f.endswith('.json') or f.endswith('.pickle'):
                shutil.copyfile(f, os.path.join(config_dir, f))


    import tomominer.pick.particle_picking_dog as PPPD
    import tomominer.pick.particle_picking_dog__filter as PPPDF
    import tomominer.pick.particle_picking_dog__permutated_peaks as PPPDP
    import tomominer.pick.particle_picking_dog__permutated_peaks__filter as PPDPPF

    # Because 'peak' in config file, program will peform dog-particle-picking on simulated tomogram
    if 'peak' in op:
        print '# perform DoG peak detection'
        # performs particle picking in single job (no partitioning, if in future we have bigger simulated tomograms, then we can use peak_partition function instead of peak, but then we also need to take care of overpartitions, duplicated particles, i.e. filtering is required)
        peak_pp = PPPD.peak(v=vt, s1=op['peak']['sigma1'], s2 = op['peak']['sigma1']*op['peak']['sigma_k'], find_maxima=op['density_positive'])

        if op['density_positive']:
            peak_pp = sorted(peak_pp, key=lambda _:(-_['val']))     # order according to peak value
        else:
            peak_pp = sorted(peak_pp, key=lambda _:(_['val']))     # order according to peak value

        # assign id to particles
        for i in range(len(peak_pp)):       peak_pp[i]['id'] = i

        filter_op = {}
        if 'top_num' in op['peak']:     filter_op['top_num'] = op['peak']['top_num']
        # filtering peaks on the basis of inter-peak distance, peaks are retained if distance between peaks is greater than sigma1
        peak_ppf = PPPDF.do_filter(pp=peak_pp, peak_dist_min=op['peak']['sigma1'], op=filter_op)
        assert      len(peak_ppf) > 0

        # permute tomogram and pick particles, vol_size is size of tomogram that will be permuted
        # perm_peak returns the list of values of peaks picked in permuted tomogram
        perm_p_val = PPPDP.perm_peak(v=vt, vol_size=op['peak']['permute']['vol_size'], s1=op['peak']['sigma1'], s2=op['peak']['sigma1']*op['peak']['sigma_k'], find_maxima=op['density_positive'])
        peak_ppf = PPDPPF.filter_peaks(pp=peak_ppf, v=perm_p_val, std_threshold_factor=op['peak']['permute']['filter']['std_threshold_factor'], find_maxima=op['density_positive'])
        assert      len(peak_ppf) > 0

        peak_re = {'x':N.array([_['x'] for _ in peak_ppf]), 'p':[_['val'] for _ in peak_ppf]}

        # todo: generate a configure file in /tmp for visual inspection using particle_picking_dog__display__imod.py
        
        print '# correspond particle picking with ground truth'
        # peak_true_label returns pdb_id for each picked peak by checking in which bounding sphere it falls (i.e. predicts label depending on distance from true locations)
        true_lbl = PPPD.peak_true_label(xp=peak_re['x'], xt_rec=models[op['model_id']]['instances'], rad=bounding_spheres)

        print '# use particle picking\'s DoG to set subtomogram size'

        # if size_ratio is given in config file then use sigma1(radius)*2*size_ratio as length of subtomogram, else use size directly if mentioned in config file
        if 'size_ratio' in op['subtomogram']:
            assert      'size' not in op['subtomogram']
            # size should be 3 dimentional with each length = diameter*ratio
            # here 2*sigma1 is diameter and 3 outside square brackets means 3 dimentional (short cut :) )
            op['subtomogram']['size'] = [  int(N.ceil(2.0 * op['peak']['sigma1'] * op['subtomogram']['size_ratio']))  ] * 3
            print '# use sigma1 subtomogram size', op['subtomogram']['size']
 
        elif 'size' in op['subtomogram']:
            assert      'size_ratio' not in op['subtomogram']
            op['subtomogram']['size'] = [    op['subtomogram']['size']    ] * 3

        # extract subtomograms
        print '# extract subtomograms',

        import tomominer.pick.subtomogram_extraction as PSE
        
        subtomogram_paths = PSE.subtomogram_extraction(v=vt, x=peak_re['x'], siz=N.array(op['subtomogram']['size']), out_dir=out_dir)
        print len(subtomogram_paths), 'subtomograms extracted'

    # if 'peak' not in config file then no particle picking, just extract ground truth
    else:
        
        print '# obtain peaks from ground truth, then extract subtomograms'
        instances = models[op['model_id']]['instances']
        peak_re = {}
        peak_re['x'] = [_['x'] for _ in instances];     peak_re['x'] = N.array(peak_re['x'])
        true_lbl = [_['pdb_id'] for _ in instances]
        peak_re['r'] = [bounding_spheres[_]['r'] for _ in true_lbl]
        peak_re['angle'] = [_['angle'] for _ in instances]

        # here we exract subtomograms of varying sizes, depending on diameter of complex. Finally all the subtomograms will be made of same size by appending zeros/ noise on all sides
        
        if 'size_ratio' in op['subtomogram']:
            assert      'size' not in op['subtomogram']
            # save diameter of biggest particle
            diameter_max = max( bounding_spheres[_]['r'] for _ in bounding_spheres ) * 2
            op['subtomogram']['size'] = [  int(N.ceil(diameter_max * op['subtomogram']['size_ratio']))  ] * 3
            #print '# use max diameter of ground truth complexes to set subtomogram size', op['subtomogram']['size']
        
        elif 'size' in op['subtomogram']:
            assert      'size_ratio' not in op['subtomogram']
            op['subtomogram']['size'] = [    op['subtomogram']['size']    ] * 3
        
        print '# extract subtomograms', ('with ground truth masked according to subtomogram__ground_truth_mask_ratio' if ('ground truth mask ratio' in op['subtomogram']) else '')

        import tomominer.pick.subtomogram_extraction_for_variable_size as PSEVS
        
        subtomogram_paths = PSEVS.subtomogram_extraction(v=vt, x=peak_re['x'], out_dir=out_dir, radius = N.array(peak_re['r']), op=op)
        print len(subtomogram_paths), 'subtomograms extracted'



    print '# generate and store missing wedge mask (combined with spherical mask) file'
    import tomominer.image.vol.wedge.util as IVMU
    # generate wegde mask depending on missing wedge angle and subtomogram size
    mask = IVMU.wedge_mask(op['subtomogram']['size'], op['bp']['model']['missing_wedge_angle'])


    mask_file = os.path.join(out_dir, 'wedge-mask.mrc')
    IOF.put_mrc( N.array(mask, dtype=float, order='F'), mask_file, overwrite=True)

    # creating info.json file
    print '# record information'
    info = []
    for i in range(len(true_lbl)):
        if i not in subtomogram_paths:      continue

        rec = {}
        rec['peak'] = {}
        rec['peak']['id'] = i
        rec['peak']['loc'] = [_ for _ in peak_re['x'][i]]
        if 'p' in peak_re:      rec['peak']['val'] = peak_re['p'][i]
        if 'angle' in peak_re:  rec['peak']['angle'] = peak_re['angle'][i]      # record true angle
        rec['subtomogram'] = subtomogram_paths[i]
        rec['mask'] = mask_file
        rec['tomogram_id'] = op['model_id']

        if true_lbl[i] is not None:
            rec['pdb_id'] = true_lbl[i]
            rec['cluster_label'] = maps[rec['pdb_id']]['id']

        info.append(rec)


    with open(os.path.join(out_dir, 'info.json'), 'w') as f:        json.dump(info, f, indent=2)
    del v
    del vt


def main():

    #op_file = sys.argv[1]
    op_file = 'model_tomogram_subtomogram__config.json'
    
    with open(op_file) as f:    op = json.load(f)

    with open(op['model_generation_config_file']) as f:     model_op = json.load(f)

    with open(op['back_projection_reconstruction_config_file']) as f:     bp_op = json.load(f)

    # if multiprocessing then this else done sequentially for each model number and further for each SNR.
    # So its better to do multiprocessing
    if op['multi_processing_process_num'] > 0:
        # parallel processing
        import multiprocessing
        from multiprocessing.pool import Pool

        pool = Pool(processes=min(op['multi_processing_process_num'], multiprocessing.cpu_count()))
        pool_re = []
        for model_id in range(model_op['model_number']):
            if ('selected_models' in op) and (model_id not in set(op['selected_models'])):      continue

            for snr in bp_op['model']['SNR_s']:
                opt = copy.deepcopy(op)
                opt['model_id'] = str(model_id)
                opt['bp'] = copy.deepcopy(bp_op)
                del opt['bp']['model']['SNR_s']
                opt['bp']['model']['SNR'] = snr

                #model_processing(op=op)
                pool_re.append( pool.apply_async(func=model_processing, kwds={'op':opt} )  )

        for i, r in enumerate(pool_re):
            r.wait()

        re = [_.get(999999)    for _ in pool_re]

    else:

        for model_id in range(model_op['model_number']):
            if ('selected_models' in op) and (model_id not in set(op['selected_models'])):      continue

            for snr in bp_op['model']['SNR_s']:
                opt = copy.deepcopy(op)
                opt['model_id'] = str(model_id)
                opt['bp'] = copy.deepcopy(bp_op)
                del opt['bp']['model']['SNR_s']
                opt['bp']['model']['SNR'] = snr

                model_processing(op=opt)



if __name__ == '__main__':
    
    main()

