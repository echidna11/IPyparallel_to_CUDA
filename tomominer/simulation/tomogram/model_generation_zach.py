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

import tomominer.geometry.pack.zach_sphere as GPZS
import tomominer.statistics.util as SU


def generate_model(op):

    op = copy.deepcopy(op)


    raise   Exception('the following part is not correct, need to replace with maps integer id, for consistancy with template_ground_truth.py')
    pdb_ids_v = list( op['bounding_spheres'].keys() )
    pdb_ids = { (_+1) : pdb_ids_v[_] for _ in range(len(pdb_ids_v)) }       # type id must start from 1, so ...


    # ---------------------------------------------
    # set radius
    param = op['packing']['param']

    assert      'types' not in param['molecule']
    param['molecule']['types'] = []

    for i in pdb_ids:
        param['molecule']['types'].append( {'id':i, 'mass':1.0, 'radius':op['bounding_spheres'][ pdb_ids[i] ]['r'] * op['bounding_sphere_radius_factor'], 'diffusion_rate':0.5} )



    # -------------------------------------------
    # generate and assign instance numbers
    conf = op['packing']['conf']
    if 'molecule' not in conf:  conf['molecule'] = {} 
    assert      'random' not in conf['molecule']

    prop = {}
    for i in pdb_ids:        prop[i] = float( N.random.random() )

    if "high_freq" in op:
        # we can manually define a number of high frequency complexes
        minimum_freq = float(op['high_freq']['minimum_freq'])      ;       assert minimum_freq <= 1
        for i in pdb_ids:
            if pdb_ids[i] in set(op['high_freq']['pdb_ids']):
                prop[i] = minimum_freq + (float(N.random.random()) * (1 - minimum_freq))

        del minimum_freq


    models = {}
    model_i = 0

    pool = Pool()
    while True: 
        # parallel processing

        pool_apply = []
        for i in range(multiprocessing.cpu_count()):
            pool_apply.append(         pool.apply_async(func=generate_model_single, kwds={'op':copy.deepcopy(op), 'prop':prop, 'pdb_ids':pdb_ids} )  )

        for pa in pool_apply:
            pre_filtered = pa.get(99999)

            print 'model', model_i, ' packaged particle number ', len(pre_filtered)
            if len(pre_filtered) == 0:      continue

            models[model_i] = {}
            models[model_i]['instances'] = pre_filtered
            models[model_i]['random_seed'] = random.randint(0, 2**32-1)           # seed for python, numpy and matlab. Note python and numpy uses different random seeds, needs to set seperately

            model_i += 1

            if model_i >= op['model_number']:   break

        if model_i >= op['model_number']:   break

    del pool

    with open('model_generation__out.json', 'w') as f:       json.dump(models, f, indent=2)

    # write out statistics for inspection
    instance_count = {}
    for model_i in models:
        count = {}
        for pdb_id in pdb_ids.itervalues():            count[pdb_id] = len(     [_ for _ in models[model_i]['instances'] if (_['pdb_id']) == pdb_id]     )
        instance_count[model_i] = count

    prop_pdb = {}
    for i in prop:  prop_pdb[pdb_ids[i]] = prop[i]

    with open('model_generation__stat.json', 'w') as f:       json.dump({'instance_count':instance_count, 'proportion':prop_pdb}, f, indent=2)




def generate_model_single(op, prop, pdb_ids):
    conf = op['packing']['conf']


    # generate number of instances according to proportion

    prop_seq = SU.proportion_sample(size=op['copy_number']['total'], prop=prop);      prop_seq = N.array(prop_seq, dtype=N.int)
    
    conf['molecule']['random'] = []
    for i in pdb_ids:        conf['molecule']['random'].append( { 'type':i, 'number':int(N.sum(prop_seq == i)), 'ignore_overlap':1 }  )



    # perform packing
    pre = GPZS.do_packing(op['packing'])


    # filter packing result and fill up rotation angle
    pre_filtered = []
    for c in pre['coordinates']:   
        if c['type'] not in pdb_ids: continue
             
        c['pdb_id'] = pdb_ids[ c['type'] ]
        c['x'] = [int(_) for _ in c['x']]       # use an integer location so that it is more accurate for paste instances and cut subtomgrams

        ang = 2 * N.pi * N.random.random(3)
        c['angle'] = [_ for _ in ang]

        pre_filtered.append(c)

    return pre_filtered



if __name__ == '__main__':


    # load options
    #op_file = sys.argv[1]
    op_file = 'model_generation_zach__op.json'
    with open(op_file) as f:    op = json.load(f)

    assert      op['copy_number']['total'] <= op['packing']['param']['n_max']

    # load particle radius for bounding shperes
    with open(op['bounding_sphere_file']) as f:   bounding_spheres = pickle.load(f)
    bounding_spheres = bounding_spheres['bounding_spheres']

    op['bounding_spheres'] = bounding_spheres

    op['packing']['param']['vol_high'] = op['packing']['param']['box']      # just change from imp convension to Zach's program
    del op['packing']['param']['vol_high']

    generate_model(op=op)



