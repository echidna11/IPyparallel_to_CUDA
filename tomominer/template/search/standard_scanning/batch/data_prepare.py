#!/usr/bin/env python



'''
convert tomograms from mrc format to npy format, partition if needed


~/ln/tomominer/tomominer/template/search/standard_scanning/batch/data_prepare.py
'''

import os, sys, json
import cPickle as pickle
import numpy as N

import tomominer.io.file as IF

def main():
    with open('config_prepare__op.json') as f:    op = json.load(f)
    with open(op['config_stat_out_file']) as f:     st_out = json.load(f)


    print 'uploading template'
    v = IF.read_mrc_vol(op['template'])
    N.save(st_out['template_tmp'], v)


    with open(st_out['jobs_file'], 'rb') as f:        bl = pickle.load(f)['partitions']

    v = IF.read_mrc_vol(op['tomogram'])

    # load large map, then store partition files
    print 'uploading partitions'
    vp = None
    for i, b in enumerate(bl):
        bp = b['base']
        if not os.path.isdir(b['dir']):            os.makedirs(b['dir'])
        if not os.path.isfile(b['map_file']):
            vp = v[bp[0,0]:bp[0,1], bp[1,0]:bp[1,1], bp[2,0]:bp[2,1]]
            N.save(b['map_file'], vp)
        print '\r %5d / %5d    '%(i, len(bl)),      ;       sys.stdout.flush()
    


if __name__ == '__main__':
    main()

