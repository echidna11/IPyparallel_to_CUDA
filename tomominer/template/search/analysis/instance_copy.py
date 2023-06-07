#!/usr/bin/env python


'''
copy rigid transformed instances to a folder, according to decreasing order of alignment scores, 

~/ln/tomominer/tomominer/template/search/analysis/instance_copy.py
'''


import os, sys, json

import numpy as N

import tomominer.io.file as IF
import tomominer.geometry.rotate as GR
import tomominer.filter.gaussian as FG


if __name__ == '__main__':

    with open('template_search_analysis_instance_copy__op.json') as f:      op = json.load(f)

    if not os.path.isdir(op['out_dir']):      os.makedirs(op['out_dir'])

    with open(op['data_file']) as f:    dj = json.load(f)

    dj = sorted(dj, key = lambda _:_['score'], reverse=True)

    for i, d in enumerate(dj):
        v = IF.read_mrc_vol(d['subtomogram'])

        if 'gaussian_sigma' in op:      v = FG.smooth(v, op['gaussian_sigma'])

        v = GR.rotate(v, angle=N.array(d['angle']), loc_r=N.array(d['loc']))

        IF.put_mrc(v, os.path.join(op['out_dir'], '%06d-%04d.mrc'%(i,d['tomogram_id'])))    

        print '\r', i, '            ',            ;       sys.stdout.flush()



