#!/usr/bin/env python


import numpy as N
import tomominer.geometry.rotate as GR
import tomominer.image.vol.util as CV


'''
convert a model to a density map
parameters:   size: container size / map size     model : center location and rotation angles of particles            bounding_sphere : bounding_sphere information       maps : density map of individual complexes
'''

def get_map(size, model, bounding_spheres, maps):

    v = N.zeros(shape=size, dtype=N.float)

    for i in range(len(model)):
        mdl = model[i]
        pid = mdl['pdb_id']
        vr = GR.rotate( v=maps[ pid ]['map'], angle=N.array(mdl['angle']), c1=N.array(bounding_spheres[pid]['c']), default_val=0.0 )
        vr[vr < 0.0] = 0.0

        CV.add_to_whole_map(whole_map=v, vol=vr, c=N.array(mdl['x']))


    return v


if __name__ == '__main__':

    import os
    import sys
    import json
    import pickle
    import scipy.io as SIO



    op_file = sys.argv[1]
   
    with open(op_file) as f:    op = json.load(f)


    with open(op['templete_select_file']) as f:     maps = pickle.load(f)

    with open(op['model_generation_config_file']) as f:     model_op = json.load(f)
    with open(op['model_file']) as f:       models = json.load(f)

    with open(op['bounding_sphere_file']) as f:     bounding_spheres = pickle.load(f)
    bounding_spheres = bounding_spheres['bounding_spheres']


    vol_high = model_op['packing']['param']['vol_high']
    map_size = N.array( [vol_high['x'], vol_high['y'], vol_high['z']] )


    import scipy.io as SIO
    def model_to_map_local(map_size, model_id, model, bounding_spheres, maps, out_dir):
        v = get_map(size=map_size, model=models[str(model_id)]['instances'], bounding_spheres=bounding_spheres, maps=maps)
        out_file = os.path.join(out_dir, 'model-%s.mat'%(model_id))
        SIO.savemat(out_file, {'v':v})

        return out_file


    from multiprocessing.pool import Pool

    pool = Pool()
    pool_re = []
    for model_id in models:
        pool_re.append( pool.apply_async(func=model_to_map_local, kwds={'map_size':map_size, 'model_id':model_id, 'model':models[model_id]['instances'], 'bounding_spheres':bounding_spheres, 'maps':maps, 'out_dir':op['out_dir']}) )

    for r in pool_re:       r.get()




