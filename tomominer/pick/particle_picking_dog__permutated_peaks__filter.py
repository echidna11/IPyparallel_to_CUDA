#!/usr/bin/env python



'''
particle_picking_dog__permutated_peaks__filter.py

given particle_picking_dog__permutated_peaks.py result, calculate a threshold, and use it to filter out peaks

'''

import numpy as N

# parameters: pp list of peaks from original data,    v list of permuted peak values
def filter_peaks(pp, v, std_threshold_factor, find_maxima):
    v = N.array(v)

    if find_maxima:
        t = v.mean() + std_threshold_factor*v.std()
        return [_ for _ in pp if (_['val'] > t)]
    else:
        t = v.mean() - std_threshold_factor*v.std()
        return [_ for _ in pp if (_['val'] < t)]
    


if __name__ == '__main__':
    
    import json
    
    with open('particle_picking_dog__permutated_peaks__filter__op.json') as f:      op = json.load(f)

    print 'std_threshold_factor', op['std_threshold_factor']

    with open('particle_picking_dog__config.json') as f:    pp_op = json.load(f)

    with open(op['real_peak_file']) as f: pp = json.load(f)

    import cPickle as pickle
    with open(op['permutated_peak_file'], 'rb') as f:   perm_pp = pickle.load(f)


    pp_f = []
    for tom_i in perm_pp:
        # calculate a threshold
        pp_f.extend(        filter_peaks(pp=[_ for _ in pp if (_['tomogram_id'] == tom_i)], v=perm_pp[tom_i], std_threshold_factor=float(op['std_threshold_factor']), find_maxima=pp_op['find_maxima'])      )

    print 'peak number reduced from', len(pp), 'to', len(pp_f)
    with open('particle_picking_dog__permutated_peaks__filter__out.json', 'w') as f:     json.dump(pp_f, f, indent=2)

