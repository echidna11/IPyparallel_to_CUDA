#!/usr/bin/env python

# find the Difference of Gauss function that best approximates a template through parameter scanning


import json
import numpy as N

import tomominer.io.file as IF
import tomominer.geometry.rotate as GR
import tomominer.model.util as MU
import tomominer.filter.convolve as FC


if __name__ == '__main__':
    
    with(open('template_gaussian_approximation__op.json')) as f:     op = json.load(f)

    t = IF.read_mrc(op['template_file'])['value'].astype(N.float)

    if not op["intensity_positive"]:    t = -t

    if 'enlarge_factor' in op:      t = GR.rotate(t, siz2=N.round(N.array(t.shape)*op['enlarge_factor']), default_val=t.mean())        # sometimes when the template volume is tightly bounded to the complex enlarge it a bit

    
    sig1_s = N.arange(start=op['sigma1']['start'], stop=op['sigma1']['end'], step=op['sigma1']['step'])
    best = {}
    for s1 in sig1_s:
        for k in op['k_factors']:
            g = MU.difference_of_gauss_function(t.shape, s1, s1 * k)
            c = FC.pearson_correlation_simple(t, g)
            
            c_max = c.max()

            if ('c' in best) and (c_max <= best['c']):      continue

            best['c'] = c_max
            best['s1'] = s1
            best['k'] = k

            print c_max, s1, k


