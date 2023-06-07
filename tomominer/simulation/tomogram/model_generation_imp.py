#!/usr/bin/env python

# given bounding sphere results, pack the spheres and models containing location and orientations
# the packing is done using Zach's packing method

import os
import sys
import json
import pickle
import copy
import random
import multiprocessing
from multiprocessing.pool import Pool
import numpy as N
import tomominer.geometry.pack.sphere.imp.imp as GPSII
import tomominer.statistics.util as SU


def generate_model(op):
    op = copy.deepcopy(op)
    # ---------------------------------------------
    # set radius
    param = op['packing']['param']
    param['molecule'] = {}
    param['molecule']['types'] = {}

    # use maps integer id, for consistancy with template_ground_truth.py
    bs = op['bounding_spheres']
    bs = {bs[_]['id']:bs[_] for _ in bs}        # just change index from pdb id to integer id
    for i in bs:
        param['molecule']['types'][i] = {'id':i, 'mass':1.0, 'radius':bs[i]['r'], 'diffusion_rate':0.5}

    # -------------------------------------------
    # generate and assign instance numbers
    pdb_ids={_:bs[_]['pdb_id'] for _ in bs}

    prop = {}
    if 'reference_prop' in op:
        print 'use given frequency'
        for i in pdb_ids:       prop[i] = op['reference_prop'][pdb_ids[i]]

    else:
        for i in pdb_ids:        prop[i] = float( N.random.random() )

        if "high_freq" in op:
            # we can manually define a number of high frequency complexes
            minimum_freq = float(op['high_freq']['minimum_freq'])      ;       assert minimum_freq <= 1
            for i in pdb_ids:
                if pdb_ids[i] in set(op['high_freq']['pdb_ids']):
                    prop[i] = minimum_freq + (float(N.random.random()) * (1 - minimum_freq))

            del minimum_freq


    pymol_dir = os.path.join(os.getcwd(), 'pymol')
    if not os.path.isdir(pymol_dir):        os.makedirs(pymol_dir)

    models = {}
    model_i = 0

    pool = Pool()
    while True: 
        # parallel processing
        pool_apply = []
        for i in range(min(multiprocessing.cpu_count(), op['model_number'])):
            opt = copy.deepcopy(op)
            opt['model id'] = i
            opt['pymol file'] = os.path.join(pymol_dir, '%d.pym'%(i,))
            opt['random_seed'] = random.randint(0, 2**32-1)           # seed for python, numpy and matlab. Note python and numpy uses different random seeds, needs to set seperately
            pool_apply.append(pool.apply_async(func=generate_model_single, kwds={'op':opt, 'prop':prop, 'pdb_ids':pdb_ids}))
        for pa in pool_apply:
            pre = pa.get(99999)
    
            print 'model', model_i, ' packaged particle number ', len(pre['conf'])
            if len(pre['conf']) == 0:   continue
            models[model_i] = {}
            models[model_i]['instances'] = pre['conf']
            models[model_i]['pack score'] = pre['pack score']
            models[model_i]['pack temprature'] = pre['pack temprature']
            models[model_i]['pack inside box num'] = pre['pack inside box num']
            models[model_i]['random_seed'] = pre['random_seed']
            model_i += 1
            if model_i >= op['model_number']:   break

        if model_i >= op['model_number']:   break

    del pool

    with open('model_generation_imp__out.json', 'w') as f:       json.dump(models, f, indent=2)

    # write out statistics for inspection
    instance_count = {}
    for model_i in models:
        count = {}
        for pdb_id in pdb_ids.itervalues():            count[pdb_id] = len(     [_ for _ in models[model_i]['instances'] if (_['pdb_id']) == pdb_id]     )
        instance_count[model_i] = count

    prop_pdb = {}
    for i in prop:  prop_pdb[pdb_ids[i]] = prop[i]

    with open('model_generation_imp__stat.json', 'w') as f:       json.dump({'instance_count':instance_count, 'proportion':prop_pdb}, f, indent=2)




def generate_model_single(op, prop, pdb_ids):

    # seed for python, numpy and matlab. Note python and numpy uses different random seeds, needs to set seperately
    random.seed(op['random_seed'])
    N.random.seed(op['random_seed'])


    # generate number of instances according to proportion

    prop_seq = SU.proportion_sample(size=op['copy_number']['total'], prop=prop);      prop_seq = N.array(prop_seq, dtype=N.int)

    types = op['packing']['param']['molecule']['types']
    
    conf = []
    for i in range(len(prop_seq)):        conf.append( { 'id':i, 'type':prop_seq[i], 'mass':types[prop_seq[i]]['mass'],  'radius':float(types[prop_seq[i]]['radius']), 'diffusion_rate':types[prop_seq[i]]['diffusion_rate'] }  )
    
    # perform packing
    pre = GPSII.do_packing(conf=conf, op=op['packing']['param'], pymol_file_name=op['pymol file'])

    # filter packing result and fill up rotation angle
    pre_filtered = []
    for c in pre['conf']:   
        if c['type'] not in pdb_ids: continue
             
        c['pdb_id'] = pdb_ids[ c['type'] ]
        c['x'] = [int(_) for _ in c['x']]       # use an integer location so that it is more accurate for paste instances and cut subtomgrams

        ang = 2 * N.pi * N.random.random(3)
        c['angle'] = [_ for _ in ang]

        pre_filtered.append(c)

    return {'conf':pre_filtered, 'pack score':pre['score'], 'pack temprature':pre['temprature'], 'pack inside box num':pre['inside box num'], 'random_seed':op['random_seed']}




def main():

    # load options
    #op_file = sys.argv[1]
    op_file = 'model_generation_imp__op.json'
    with open(op_file) as f:    op = json.load(f)

    # load particle radius for bounding shperes
    with open(op['bounding_sphere_file']) as f:   bounding_spheres = pickle.load(f)
    bounding_spheres = bounding_spheres['bounding_spheres']

    op['bounding_spheres'] = bounding_spheres

    
    if 'reference_prop' in op:
        with open(op['reference_prop']) as f:       rp = json.load(f)['proportion']           # the reference proportion is usually inside model_generation_imp__stat.json
        op['reference_prop'] = rp

    generate_model(op=op)



if __name__ == '__main__':
    main()

