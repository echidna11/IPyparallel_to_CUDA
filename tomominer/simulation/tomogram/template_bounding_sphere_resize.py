#!/usr/bin/env python



'''
resize bounding sphere so that the density map can be more compact
'''

import json, pickle, copy

if __name__ == "__main__":
    with open('template_bounding_sphere_resize__op.json') as f:     op = json.load(f)

    with open(op['input file'], 'rb') as f:         bs = pickle.load(f)['bounding_spheres']

    for p, b in bs.iteritems():     b['r'] *= op['ratio']

    with open(op['output file'], 'wb') as f:         pickle.dump({'bounding_spheres': bs}, f)

   
