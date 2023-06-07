#!/usr/bin/env python




"""
In a model file, remove particles whose bounding sphere is outside the volume boundary
"""

'''
~/ln/tomominer/tomominer/simulation/tomogram/analysis/filter_boundary_particles.py
'''


import json, pickle

import numpy as N


def main():
    with open('filter_boundary_particles__op.json') as f:         op = json.load(f)

    with open(op['model op file']) as f:       model_op = json.load(f)
    with open(op['model file']) as f:       models = json.load(f)
    with open(op['bounding sphere file'], 'rb') as f:       bounding_spheres = pickle.load(f)['bounding_spheres']

    box  = model_op['packing']['param']['box']
    map_size = N.array( [box['x'], box['y'], box['z']] )


    models_new = {}
    for model_id in models:
        model_t = []

        for m in models[model_id]['instances']:
            x = m['x']
            pid = m['pdb_id']
            r = bounding_spheres[pid]['r']
            if N.any( (x-r) <= 0 ):     continue
            if N.any( (x+r) >= map_size ):  continue
            model_t.append(m)

        models_new[model_id] = {'instances': model_t}

    with open(op['model file out'], 'w') as f:      json.dump(models_new, f, indent=2)


if __name__ == '__main__':
    main()

