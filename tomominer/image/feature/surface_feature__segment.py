#!/usr/bin/env python

# given surface feature calculated from surface_feature.py, set a cutoff using mean value, and use cutoff to keep the largest connected components


if __name__ == '__main__':

    import json
    with open('surface_feature__segment__op.json') as f:    op = json.load(f)
    with open('surface_feature__op.json') as f:     ss_op = json.load(f)


    import tomominer.io.file as IF
    r = IF.read_mrc(ss_op['out_file'])['value']

    m =        r > (r.mean() + r.std()*op['threshold_std_ratio'])
    IF.put_mrc(m, '/tmp/m.mrc', overwrite=True)

    import time
    cur_tim = time.time()    
    import tomominer.segmentation.connected_regions as SC
    cr = SC.connected_regions(m)
    print 'connected_regions()', time.time() - cur_tim, 'sec'


    print 'calculating and store segment sizes'
    siz = [0] * (cr['max_lbl']+1)
    for l in range(1, cr['max_lbl']+1):        
        siz[l] = (cr['lbl'] == l).sum()
        #print l, siz[l]

    import numpy as N
    siz_i = N.argsort(- N.array(siz) )

    for l in siz_i:        
        if siz[l] < op['segment_size_min']:     continue
        print l, siz[l], '\t',

    s = N.zeros(m.shape)
    for l in range(len(siz)):
        if siz[l] < op['segment_size_min']:     continue
        s[cr['lbl'] == l] = 1

    IF.put_mrc(s, op['out_file'], overwrite=True)

