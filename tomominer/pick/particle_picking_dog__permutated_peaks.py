#!/usr/bin/env python



'''
randomly sample values inside a tomogram and form smaller volumes, then apply the same DoG detect peaks inside such volumes, and collect peak values for later use to filter out peaks in particle_picking_dog__out.json before do distance filtering.

'''

import numpy as N

import tomominer.io.file as IOF
import tomominer.pick.particle_picking_dog as PP

def perm_peak(v, vol_size, s1, s2, find_maxima):
    v = N.random.choice(v.flatten(), size=tuple(vol_size)) 

    ps = PP.peak(v=v.astype(N.float32), s1=s1, s2=s2, find_maxima=find_maxima)

    return [_['val'] for _ in ps]


if __name__ == '__main__':

    import json
    with open('particle_picking_dog__permutated_peaks__op.json') as f:      op = json.load(f)

    with open('particle_picking_dog__config.json') as f:    pp_op = json.load(f)

    with open(pp_op['tomogram_info_file']) as f:                   tom_info = json.load(f)

    
    peaks = {} 
    for tom in tom_info:
        if ('selected_tomogram_ids' in pp_op) and (tom['id'] not in set(pp_op['selected_tomogram_ids'])):     continue

        print 'loading tomogram', tom['id'],
        v = IOF.read_mrc(tom['vol_file'])['value']
        print 'done'

        sigma1 = pp_op['sigma1'] / tom['voxel_spacing']
        peaks[tom['id']] = perm_peak(v=v, vol_size=op['vol_size'], s1=sigma1, s2=sigma1*pp_op['sigma_k'], find_maxima=pp_op['find_maxima'])        # convert sigma1 from nm unit to voxel unit, according to voxel spacing


    import cPickle as pickle
    with open('particle_picking_dog__permutated_peaks__out.pickle', 'wb') as f:     pickle.dump(peaks, f, protocol=-1)

