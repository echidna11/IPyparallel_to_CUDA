#!/usr/bin/env python



# perform level set segmentation, then make density map to display connected components

import os
import json
import numpy as N

import tomominer.io.file as IF
import tomominer.filter.gaussian as FG
import tomominer.pick.segmentation.level_set.util as PSLU
import tomominer.segmentation.util as SU


def main():

    with open('pick_segmentation_level_set_visual_check__op.json') as f:     op = json.load(f)
    with open(op['segmentation op file']) as f:     seg_op = json.load(f)

    if not os.path.isfile(op['clip']['tomogram file']):
        print 'loading', op['input tomogram file']

        v = IF.read_mrc(op['input tomogram file'])['value']
        cs = op['clip']['start']
        ce = op['clip']['end']
        v = v[cs[0]:ce[0], cs[1]:ce[1], cs[2]:ce[2]]
        
        IF.put_mrc(v, op['clip']['tomogram file'], overwrite=True)

    else:
        print 'loading', op['clip']['tomogram file']
        v = IF.read_mrc(op['clip']['tomogram file'])['value']

    if not seg_op['intensity positive']:    v = -v

    if 'smooth sigma' in seg_op:
        print 'smoothing', seg_op['smooth sigma']
        v = FG.smooth(v, seg_op['smooth sigma'])

    IF.put_mrc(v, op['out clip file'], overwrite=True)

    print 'segmenting'
    phi = PSLU.segment(v=v, op=seg_op)
    IF.put_mrc(phi, op['out phi file'], overwrite=True)
 
    seg = phi>0
    IF.put_mrc(seg, op['out segmentation file'], overwrite=True)


    print 'finding connected components'
    cc = SU.connected_components(seg)
    IF.put_mrc(cc['label'], op['out component file'], overwrite=True)


if __name__ == '__main__':
    main()


