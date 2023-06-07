

import copy
import numpy as N

import tomominer.pursuit.multi.util as PMU
import tomominer.statistics.ssnr as SS

import tomominer.pursuit.multi.util as PMU
import tomominer.io.file as IV
import tomominer.align.refine.gradient_refine as AFGF


'''
given template alignment and best match information, sequencially include subtomograms of best alignment and sequentially calculate FSC scores. Then look at local maxima of the FSC sequence, and choose the one that has best correlation with the template
see idea 1.1 of:                https://docs.google.com/document/d/1NmYhEczJ8A--KI5PStf8hyuPM5QPfRLz7RftOfkn8_k/edit
parameters:       dj: data_json       
return:           selected data_json

code modified from tomominer.pursuit.multi.util.cluster_formation_alignment_fsc__by_local_maxima()

'''


def ssnr_sequential(self, dj, op):
    print 'ssnr_sequential()', op


    # order subtomograms according to decrease of alignment scores
    dj = sorted( copy.deepcopy(dj), key=lambda _ : float(_['score']), reverse=True )

    # parallel calculate sequential ssnr
    ssnr_s = SS.ssnr_sequential_parallel(self=self, data_json_dict={0:dj}, op=op)

    return ssnr_s[0]


def ssnr_sequential__global_maxima(ssnr_s, dj, op):

    op = copy.deepcopy(op)

    # retrive FSC score sequence
    fsc_s = N.array([N.sum(_) for _ in ssnr_s['fsc']])

    import scipy.ndimage.filters as SDF
    if 'gaussian_smooth_sigma' in op:
        fsc_s = SDF.gaussian_filter1d(fsc_s, op['gaussian_smooth_sigma'])           # sometimes there are too many local maxima, can smooth a bit to reduce


    # find local maxima and collect corresponding data_json

    if ('min_expansion_size' in op):        fsc_s[:op['min_expansion_size']] = -1.0        # remove those values that corresponds to small sets

    i = N.argmax(fsc_s)
    return {'i':i, 'data_json':copy.deepcopy(dj[:(i+1)]), 'fsc':fsc_s[i]}



def ssnr_sequential__local_maxima(ssnr_s, dj, op):

    op = copy.deepcopy(op)

    # retrive FSC score sequence
    fsc_s = N.array([N.sum(_) for _ in ssnr_s['fsc']])

    import scipy.ndimage.filters as SDF
    if 'gaussian_smooth_sigma' in op:
        fsc_s = SDF.gaussian_filter1d(fsc_s, op['gaussian_smooth_sigma'])           # sometimes there are too many local maxima, can smooth a bit to reduce


    # find local maxima and collect corresponding data_json
    dj_lm = {}

    import tomominer.filter.local_extrema as FL
    inds = FL.local_maxima(fsc_s)
    inds = inds[0]

    if ('min_expansion_size' in op):        inds = inds[inds > op['min_expansion_size']]        # remove those local maxima that corresponds to small sets
    if 'min_fs' in op:      inds = [_ for _ in inds if (fsc_s[k][_] > op['min_fsc'])]         # in such case we require the refined set average to have a higher FSC score than the original average


    if len(inds) == 0:  return dj_lm

    if op['debug']:        print 'number of local maxima', len(inds)
    
    for i in inds:          dj_lm[i] = {'i':i, 'data_json':copy.deepcopy(dj[:(i+1)]), 'fsc':fsc_s[i]}

    return dj_lm

def averaging(self, dj_lm, op):
    cav = PMU.cluster_averaging_vols(self=self, clusters={_:dj_lm[_]['data_json'] for _ in dj_lm}, op=op)
    cav = cav['cluster_avg_dict']

    for i in dj_lm:        
        if i not in cav:        continue
        dj_lm[i]['average'] = cav[i]



def correlation(dj_lm, tk):
    tv = IV.read_mrc_vol(str(tk['subtomogram']))
    tm = IV.read_mrc_vol(str(tk['mask']))
    tfv = AFGF.foster_transform_single_vol(v=tv, m=tm)
    del tv, tm

    # calculate correlation between original templates and new averages
    for i in dj_lm:
        if 'average' not in dj_lm[i]:   continue
        vf_maximum = AFGF.foster_transform_single_vol(v=dj_lm[i]['average']['vol'], m=dj_lm[i]['average']['mask'])
        dj_lm[i]['score'] = N.sum(tfv * vf_maximum)        # correlation score


def get_best_correlated_maxima(dj_lm):
    # select those achieve maximum correlations, record results
    best = None
    for i in dj_lm:
        if 'score' not in dj_lm[i]:  continue
        if (best is None) or (best['score'] < dj_lm[i]['score']):                best = dj_lm[i]


    return best
    


def ssnr_sequential_expansion__local_maxima(self, dj, op):
    print 'ssnr_sequential_expansion__local_maxima()'

    op = copy.deepcopy(op)

    ssnr_op = copy.deepcopy(op['ssnr'])
    ssnr_op['n_chunk'] = op['n_chunk']
    ssnr_s = ssnr_sequential(self, dj, op=ssnr_op)

    dj_lm = ssnr_sequential__local_maxima(ssnr_s=ssnr_s, dj=dj, op=op)

    avg_op = copy.deepcopy(op['averaging'])
    avg_op['n_chunk'] = op['n_chunk']
    averaging(self=self, dj_lm=dj_lm, op=avg_op)
    correlation(dj_lm=dj_lm, tk=op['template']['key'])
    best = get_best_correlated_maxima(dj_lm=dj_lm)

    return best



def ssnr_sequential_expansion__global_maxia(self, dj, op):
    print 'ssnr_sequential_expansion__global_maixmum()'

    op = copy.deepcopy(op)

    ssnr_op = copy.deepcopy(op['ssnr'])
    ssnr_op['n_chunk'] = op['n_chunk']
    ssnr_s = ssnr_sequential(self, dj, op=ssnr_op)

    return ssnr_sequential__global_maxima(ssnr_s=ssnr_s, dj=dj, op=op)

