#!/usr/bin/env python



# functions for particle picking using Difference of Gaussian

import os
import sys
import gc as GC
import multiprocessing
from multiprocessing.pool import Pool as Pool
import time
import numpy as N
import scipy.spatial.distance as SPD


import tomominer.filter.gaussian as FG
import tomominer.filter.local_extrema as FL
import tomominer.io.file as IOF



# peak detection
# parameters:       v: vol map      s1: sigma1      s2: sigma2
# return:       x: coordinates      p: peak values      vg: dog smoothed vol map        find_maxima: find local maxima or minima
def peak(v, s1, s2, find_maxima=True, top_num=None):

    p = peak__partition__single_job(v=v, s1=s1, s2=s2, find_maxima=find_maxima, top_num=top_num)

    return p['ps']



# partition the volume then detect peaks for each partition, note that this will result in redundant peaks!! Clean up must be done afterwards!!
def peak__partition(v, s1, s2, find_maxima=True, partition_op=None, multiprocessing_pool=None, top_num=None):

    if partition_op is None:
        # in this case, just generate a single partition
        siz_max = max(v.shape)
        partition_op = {   'nonoverlap_width': siz_max*2,    'overlap_width': siz_max*2    }
    
    import tomominer.image.vol.partition as IVP
    b = IVP.gen_bases(v.shape, nonoverlap_width=partition_op['nonoverlap_width'], overlap_width=partition_op['overlap_width'])
    print 'partition num', b.shape

    ps = []

    if multiprocessing_pool is not None:
        pool = multiprocessing_pool
        pool_re = []
        for i0 in range(b.shape[0]):
            for i1 in range(b.shape[1]):
                for i2 in range(b.shape[2]):
                    bp = N.squeeze(b[i0,i1,i2,:,:])
                    pool_re.append(      pool.apply_async(func=peak__partition__single_job, kwds={'v':v[bp[0,0]:bp[0,1], bp[1,0]:bp[1,1], bp[2,0]:bp[2,1]], 's1':s1, 's2':s2, 'base':bp, 'find_maxima':find_maxima, 'partition_id':(i0,i1,i2), 'save_vg':(partition_op['save_vg'] if 'save_vg' in partition_op else False), 'top_num':top_num } )        )

        for pool_re_t in pool_re:
            ppsj = pool_re_t.get(9999999)
            ps.extend(  ppsj['ps']  )
            print '\r', ppsj['partition_id'], '                     '     ;       sys.stdout.flush()


    else:

        for i0 in range(b.shape[0]):
            for i1 in range(b.shape[1]):
                for i2 in range(b.shape[2]):
                    bp = N.squeeze(b[i0,i1,i2,:,:])
                    ppsj = peak__partition__single_job(v=v[bp[0,0]:bp[0,1], bp[1,0]:bp[1,1], bp[2,0]:bp[2,1]], s1=s1, s2=s2, base=bp, find_maxima=find_maxima, partition_id=(i0,i1,i2), save_vg=(partition_op['save_vg'] if 'save_vg' in partition_op else False), top_num=top_num)
                    ps.extend(      ppsj['ps']   )
                    print '\r', ppsj['partition_id'], '                     '    ;       sys.stdout.flush()


    # order peaks in ps according to values
    if find_maxima:
        ps = sorted(ps, key=lambda _:(-_['val']))
    else:
        ps = sorted(ps, key=lambda _:_['val'])

    return ps



def peak__partition__single_job(v, s1, s2, base=None, find_maxima=None, partition_id=None, save_vg=False, top_num=None):
    assert  find_maxima is not None

    vg = FG.dog_smooth__large_map(v, s1=s1, s2=s2)

    if save_vg:                IOF.put_mrc(v, '/tmp/%d-%d-%d--v.mrc'%(partition_id[0], partition_id[1], partition_id[2]), overwrite=True)        # save the smoothed partition for inspection

    del v



    if find_maxima:
        #print 'local_maxima()'         ;       sys.stdout.flush()
        x = FL.local_maxima(vg)
    else:
        #print 'local_minima()'         ;       sys.stdout.flush()
        x = FL.local_minima(vg)


    p = vg[x]
    x = N.array(x).T

    if base is not None:
        assert      base.shape[0] == x.shape[1]
        assert      base.shape[1] == 2

        for dim_i in range(x.shape[1]):        x[:,dim_i] += base[dim_i,0]


    ps = []
    for i in range(len(p)):        ps.append(      {'val':float(p[i]), 'x':[_ for _ in x[i,:]]}    )

    if save_vg:              IOF.put_mrc(vg, '/tmp/%d-%d-%d--vg.mrc'%(partition_id[0], partition_id[1], partition_id[2]), overwrite=True)        # save the smoothed partition for inspection

    del vg;         GC.collect()

    # My modification for seelcting top 10000 particles in each partition
    if top_num is not None:
        if find_maxima:
            ps = sorted(ps, key=lambda _:(-_['val']))
        else:
            ps = sorted(ps, key=lambda _:(_['val']))
        ps = ps[:top_num]
    # My modification ends
    
    return {'ps':ps, 'partition_id':partition_id}





# for simulation data, if we know the true center locations of complexes, we can compare with detected peaks to assign ground truth class labels (if there is a single detected peak within certain radius of the ground truth)
# parameters:       xp : location of peaks        xt_rec : records of true location and labels following the output of model_generation.py         rad : radius, following output of template_bounding_sphere.py
# return: class labels
def peak_true_label(xp, xt_rec, rad, rad_ratio=1.0):


    xt = N.zeros( (len(xt_rec), xp.shape[1]) )

    for i, xtr in enumerate(xt_rec):
        xt[i,:] = xtr['x']


    d = SPD.cdist( xp, xt )     # distance matrix between predicted and true locations

    tf = N.zeros(d.shape, dtype=int)       # indicator of whether the distance is within certain threshold
    for i, xtr in enumerate(xt_rec):
        pid = xtr['pdb_id']

        ind =   d[:,i] < (  float( rad[ pid ][ 'r' ] ) * rad_ratio )

        tf[ind, i] = 1


    # defining true label by finding ONE to ONE correspondance of predicted peaks and true locations, defined by bounding sphere radius
    lbl = [None] * len(xp)
    for i in range(len(xp)):
        if tf[i,:].sum() != 1:      continue

        ind = N.where(tf[i,:] == 1)[0][0] # return the column where tf[i,:]==1
        # double check
        if tf[:,ind].sum() != 1:    continue

        lbl[i] = xt_rec[ind]['pdb_id']

    return lbl

'''
load a configure file, then batch processing by loading a tomogram then perform particle picking
finally output all coordinates into a json file
'''
if __name__ == '__main__':
    
    import json
    with open('particle_picking_dog__config.json') as f:    op = json.load(f)
    if 'partition' not in op:       op['partition'] = None

    if 'multiprocessing_process_num' not in op:         op['multiprocessing_process_num'] = 0
    if op['multiprocessing_process_num'] > 0:
        multiprocessing_pool = Pool(processes=min(op['multiprocessing_process_num'], multiprocessing.cpu_count()))
    else:
        multiprocessing_pool = None


    with open(op['tomogram_info_file']) as f:                   tom_info = json.load(f)



    peaks = []
    for tom in tom_info:
        if ('selected_tomogram_ids' in op) and (tom['id'] not in set(op['selected_tomogram_ids'])):     continue

        print 'loading tomogram', tom['id'],
        mrc = IOF.read_mrc(tom['vol_file'])
        print 'done'

        sigma1 = op['sigma1'] / tom['voxel_spacing']       # convert sigma1 from nm unit to voxel unit, according to voxel spacing
        ps = peak__partition(v=mrc['value'].astype(N.float32), s1=sigma1, s2=sigma1*op['sigma_k'], find_maxima=op['find_maxima'], partition_op=op['partition'], multiprocessing_pool=multiprocessing_pool, top_num=op['top_num'])
        del mrc;        GC.collect()

        for peak_i in range(len(ps)):
            peak_t = ps[peak_i]
            peak_t['tomogram_id'] = tom['id']
            peak_t['id'] = peak_i

            peaks.append(peak_t)

        with open('particle_picking_dog__out.json', 'w') as f:      json.dump(peaks, f, indent=2)       # save peaks after EVERY iteration, so that you can see these peaks early





