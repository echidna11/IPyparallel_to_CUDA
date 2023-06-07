#!/usr/bin/env python


'''
filter subtomograms according to max translation



~/ln/tomominer/tomominer/template/search/filter/translation_max.py
'''



import os, json
import numpy as N

import tomominer.io.file as IV

if __name__ == '__main__':
    with open('template_search_filter_translation_max__op.py') as f:     op = json.load(f)
    
    with open(op['input data json file']) as f:     dj = json.load(f)
    
    v = IV.read_mrc_vol(dj[0]['subtomogram'])
    rad = N.max(v.shape) / 2.0

    if 'max radius proportion' in op:
        assert 'max radius' not in op
        max_rad = rad * op['max radius proportion']
    elif 'max radius' in op:
        assert 'max radius proportion' not in op
        max_rad = op['max radius']
    else:
        raise Exception('max radius')



    djn = []
    for d in dj:
        r = N.sqrt(N.square(d['loc']).sum())
        if r > max_rad:     continue
        djn.append(d)


    with open(op['output data json file']) as f:    json.dump(djn, f, indent=2)



