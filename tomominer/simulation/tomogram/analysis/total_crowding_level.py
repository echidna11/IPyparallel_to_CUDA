#!/usr/bin/env python


'''
calculate the total crowding level (volume occupancy)

~/ln/tomominer/tomominer/simulation/tomogram/analysis/total_crowding_level.py
'''




import os, json, pickle, copy

import numpy as N


import tomominer.geometry.point_cloud.util as GPU


def main():
    with open('total_crowding_level__op.json') as f:         op = json.load(f)

    if 'contour' in op:
        print 'using contour level', op['contour'], 'for calculation'
    else:
        print 'using segmentation for calculation'

    with open(op['map file'], 'rb') as f:       maps = pickle.load(f)
    with open(op['model op file']) as f:       model_op = json.load(f)
    with open(op['model file']) as f:       models = json.load(f)
    models = {int(_):models[_] for _ in models}


    with open(op['bounding sphere file'], 'rb') as f:       bs = pickle.load(f)

    box  = model_op['packing']['param']['box']
    map_size = N.array( [box['x'], box['y'], box['z']] )

    total_vol = 0.0
    total_seg_vol = 0.0
    for model_id in models:
        print 'processing model', model_id
        for pid in maps:
            if 'contour' in op:
                v = maps[pid]['map']
                seg = (v >= op['contour'])
            else:
                seg = N.isfinite(bs['segments'][pid])

            if ('convex hull' in op) and op['convex hull']:                seg = GPU.mask_convex_hull(seg)
            
            seg_vol = float(seg.sum())
            print 'pdb id', pid, 'volume', seg_vol

            model_t = [_ for _ in models[model_id]['instances'] if (_['pdb_id'] == pid)]
            total_seg_vol += len(model_t) * seg_vol

        total_vol += N.prod(map_size)

    stat = {'total_vol':total_vol, 'total_seg_vol':total_seg_vol, 'crowding':total_seg_vol/total_vol}
    with open(op['stat out'], 'w') as f:        json.dump(stat, f, indent=2)

    print 'volume occupancy', "{0:.2f}".format(stat['crowding']*100), "%"


if __name__ == '__main__':
    main()




'''
related code

~/ln/tomominer/tomominer/simulation/tomogram/analysis/class_isosurface_plot.py
'''




'''
# export maps for inspection of surface


map_file = '../../template_select.pickle'

import pickle
with open(map_file, 'rb') as f:       maps = pickle.load(f)

import os
import tomominer.io.file as IF
for pid in maps:    IF.put_mrc(maps[pid]['map'], os.path.join('/tmp/m', pid+'.mrc'))


'''

