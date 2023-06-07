#!/usr/bin/env python


# given a tomogram, guess the tilt angle range
import numpy as N
from scipy.stats.stats import pearsonr

import tomominer.image.vol.wedge.util as W

def wedge_mask_cor(v_abs, ops):

    for op in ops:
        m = W.wedge_mask(size=v_abs.shape, ang1=op['ang1'], ang2=op['ang2'], direction=op['direction'])
        m = m.astype(N.float)

        op['cor'] = float(pearsonr(v_abs.flatten(), m.flatten())[0])

    return ops


if __name__ == '__main__':

    import warnings
    warnings.simplefilter("once")

    print 'WARNING: need to manually inspect the output best angles, in order to determine proper angles. Because sometimes the best angles are too small, in such case we want to choose less best but larger angles. Therefore we need to manually edit tomogram_prepare__out__wedge.json'
    
    import json
    with open('tilt_angle_inference__op.json') as f:        op = json.load(f)
    with open(op['tomogram_info_file']) as f:     tom_info = json.load(f)

   
    if 'processes' in op:
        n_proc = op['processes']
    else:
        import multiprocessing
        n_proc = multiprocessing.cpu_count()

    print 'using', n_proc, 'cpus'

    from multiprocessing.pool import Pool
    pool = Pool(processes=n_proc)
    pool_results = []

    import os, sys, gc
    from numpy.fft import fftshift, fftn
    import tomominer.io.file as IF
    for tom in tom_info:
        if 'selected_tomogram_ids' in op:
            if tom['id'] not in set(op['selected_tomogram_ids']):        continue

        print 'tomogram', tom['id'], tom['vol_file'],
        v = IF.read_mrc(tom['vol_file'], show_progress=True)['value']
        print 'size', v.shape

        if 'vol_clip' in op:
            vol_clip = op['vol_clip']
            v = v[vol_clip[0][0]:vol_clip[0][1], vol_clip[1][0]:vol_clip[1][1], vol_clip[2][0]:vol_clip[2][1]]
            print 'clipping', vol_clip, 'the clip must as thick (and cube) as possible, and contain as many image figures as possible'

        v = v.astype(N.float)
        v = N.log(  N.abs(      fftshift(    fftn(v)  )   ) )       # it is very important to use log scale!!

        gc.collect()        # clean up memory before forking

        
        wedge_angle_scan_range = op['wedge_angle_scan_range']
        tasks = []
        for ang1 in range(-wedge_angle_scan_range[1], -wedge_angle_scan_range[0]+1):
            for ang2 in range(wedge_angle_scan_range[0], wedge_angle_scan_range[1]+1):
                for direction in range(3):
                    tasks.append({'ang1':ang1, 'ang2':ang2, 'direction':direction})

        n_chunk = 20
        while tasks:
            #wedge_mask_cor(v_abs=v, ops=tasks[:n_chunk])
            pool_results.append(    pool.apply_async(   func=wedge_mask_cor, kwds={'v_abs':v, 'ops':tasks[:n_chunk]}    )    )
            tasks = tasks[n_chunk:]


        best = None
        for re in pool_results:
            for r in re.get(9999999):
                print '\r', r['ang1'], r['ang2'], r['direction'], r['cor'], '         ',            ;       sys.stdout.flush()

                if best is None:
                    best = r
                    continue

                if r['cor'] > best['cor']:
                    best = r
                    print

        assert  best is not None    
        tom['wedge'] = best

        # write to temp dir for visual inspection
        IF.put_mrc(v, os.path.join('/tmp', '%d-abs.mrc'%(tom['id'])), overwrite=True)
        IF.put_mrc(W.wedge_mask(size=v.shape, ang1=best['ang1'], ang2=best['ang2'], direction=best['direction']), os.path.join('/tmp', '%d-mask.mrc'%(tom['id'])), overwrite=True)


        with open(op['output_file'], 'w') as f:        json.dump(tom_info, f, indent=2)



