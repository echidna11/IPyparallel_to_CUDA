#!/usr/bin/env python



'''
keep only those particles in a model that are suficiently close to the particles of detected patterns, so it is easier for us to perform visual comparison



~/ln/tomominer/tomominer/simulation/tomogram/analysis/filter_reference_particles.py

'''

import os, json, pickle
import numpy as N

from scipy.spatial.distance import cdist


def main():
    with open('filter_reference_particles__op.json') as f:         op = json.load(f)

    with open(op['data']) as f:         dj = json.load(f)

    with open(op['model']) as f:        models = json.load(f)
    models = {int(_):models[_] for _ in models}

    with open(op['bounding sphere'], 'rb') as f:       bounding_spheres = pickle.load(f)['bounding_spheres']


    models_new = {}
    for model_id in models:
        model_t = []

        p = [_['peak']['loc'] for _ in dj if (int(_['tomogram_id']) == model_id)]
        p = N.array(p)

        for m in models[model_id]['instances']:
            x = N.array(m['x']).reshape((1,3))

            pid = m['pdb_id']
            r = bounding_spheres[pid]['r']

            if cdist(p, x).min() > r:   continue

            model_t.append(m)
            print '\r', model_id, len(model_t), '        ',

        models_new[model_id] = {'instances': model_t}

    with open(op['model out'], 'w') as f:      json.dump(models_new, f, indent=2)


if __name__ == '__main__':
    main()

