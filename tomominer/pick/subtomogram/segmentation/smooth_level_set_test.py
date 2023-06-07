#!/usr/bin/env python



'''
objective: to see which gaussian smoothing scale is most suitable to be combined with level set segmentation

given each gaussian smoothing, perform a level set segmentation, then calculate optimal score that uses the intensity difference of the original subtomogram (instead of the smoothed map)
'''

import json
import copy
import numpy as N


import tomominer.io.file as IF
import tomominer.filter.gaussian as FG
import tomominer.segmentation.active_contour.chan_vese.segment as SACS

def smooth_segmentation_score(v, op):
    vg = FG.smooth(v, sigma=op['sigma'])
    phi = SACS.segment_with_postprocessing(vg, op=op)

    s_intensity = 0.0

    ind =   phi <= 0
    s_intensity += N.square(v[ind] - v[ind].mean()).sum()

    ind =   phi > 0
    s_intensity += N.square(v[ind] - v[ind].mean()).sum()

    s_boundary = (phi == 0).sum()

    return      op['smooth_weight'] * s_boundary + (op['image_weight'] / v.var()) * s_intensity





def main():
    with open('pick_subtomogram_segmentation_smooth_level_set_test__op.json') as f:     op = json.load(f)

    v = IF.read_mrc(op['subtomogram file'])['value']

    sigma_s = N.arange(start=op['sigma']['start'], stop=op['sigma']['stop'], step=op['sigma']['step'])

    score = {}
    for sigma in sigma_s:
        ops = copy.deepcopy(op['smooth segmentation'])
        ops['sigma'] = sigma

        score[sigma] = smooth_segmentation_score(v, op=ops)
        print sigma, score[sigma]
        


if __name__ == '__main__':
    main()

