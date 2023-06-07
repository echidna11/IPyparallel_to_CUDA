import json, os, copy, sys, random, gc, time
import numpy as N
import numpy.fft as NF
import tomominer.io.file as IF
import tomominer.geometry.rotate as GR
import scipy.ndimage as SN
from scipy.stats import pearsonr
from scipy.ndimage.filters import gaussian_filter as SNFG
import warnings, math
from random import randrange
from tomominer.io.cache import Cache
warnings.filterwarnings("ignore")

def min_max(x):
    x = (x-x.min())/(x.max()-x.min())
    return x

def mask_min_max(x, t):
    m = N.mean(x)
    sd = N.std(x)
    tmp = m - t*sd
    mask = x.copy()
    mask[mask<=tmp] = -1.0
    mask[mask>tmp] = -0.0
    mask = -mask
    return mask

def mask_normalize(x, t):
    m = N.mean(x)
    sd = N.std(x)
    tmp = m - t*sd
    mask = x.copy()
    mask[mask>tmp] = 0
    mask[mask<=tmp] = 1
    return mask

def normalize(x):
    if x.std()==0:
        return x
    else:
        return (x-x.mean())/x.std()

def get_significant_points(v):
    """
    Retrieve all points with a density greater than one standard deviation above the mean.
    Return:
        An array of 4-tuple (indices of the voxels in x,y,z format and density value) 
    """
    sig_points = []
    boo = v > (v.mean() + v.std())
    for z in range(v.shape[0]):
        for y in range(v.shape[1]):
            for x in range(v.shape[2]):
                if boo[z][y][x]:
                    sig_points.append(N.array([z,y,x,v[z][y][x]]))
    return N.array(sig_points)

def get_random_significant_pairs(v, amount):
    """
    Return an array of tuple pairs of significant points randomly chosen from 'get_significant_points' function.
    For use with the DCCF and DLSF scoring functions.
    
    Arguments:
        *amount*
            number of significant point pairs to return.
    """
    sig_points = get_significant_points(v)
    sig_pairs = []
    size = len(sig_points)
    if amount <= size*size:
        random_pairs = []
        for r in range(amount):
            tmp = (randrange(size), randrange(size))
            if tmp not in random_pairs:
                fst = sig_points[tmp[0]]
                snd = sig_points[tmp[1]]
                new_value = N.array([fst[0], fst[1], fst[2], snd[0], snd[1], snd[2], fst[3]-snd[3]])
                sig_pairs.append(new_value)
    else:
        for r in range(amount):
            fst = sig_points[randrange(size)]
            snd = sig_points[randrange(size)]
            new_value = N.array([fst[0], fst[1], fst[2], snd[0], snd[1], snd[2], fst[3]-snd[3]])
            sig_pairs.append(new_value)
    return N.array(sig_pairs)

def MI(v1, v2, layers1=20, layers2=20, mask_array=None, normalised=False):
    # based on mask_array MI calculated on complete map (All 1 mask), Overlap region (AND on masks)
    
    if mask_array is None:
        mask_array = N.ones(v1.shape, dtype=int)
    else:
        mask_sum = N.sum(mask_array)
        if mask_sum == 0:
            return 0.0

    v1_copy = v1.copy()
    v2_copy = v2.copy()
    v1_copy = v1_copy*mask_array
    v2_copy = v2_copy*mask_array
    # for us the formule for number of layers gives the same value 20, as default value. :)
    # sturges rule provides a way of calculating number of bins : 1+math.log(number of points)
    layers1 = int(1 + math.log(v1.size, 2))
    layers2 = int(1 + math.log(v2.size, 2))
    layers1 = max(layers1,8)
    layers2 = max(layers2,8)

    P, _, _ = N.histogram2d(v1_copy.ravel(), v2_copy.ravel(), bins=(layers1, layers2))
    P /= P.sum()
    p1 = P.sum(axis=0)
    p2 = P.sum(axis=1)
    p1 = p1[p1 > 0]
    p2 = p2[p2 > 0]
    P = P[P > 0]
    Hx_ = (-p1 * N.log2(p1)).sum()
    Hy_ = (-p2 * N.log2(p2)).sum()
    Hxy_ = (-P * N.log2(P)).sum()
    del P, p1, p2, v1_copy, v2_copy
    gc.collect()
    
    if normalised:
        if Hxy_ == 0.0:
            return 0.0
        else:
            return (Hx_ + Hy_)/Hxy_
    else:
        return Hx_ + Hy_ - Hxy_

def cross_correlation(x, y):
    #if x and y are normalised: i.e. mean = 0 and std = 1 or x and y are centered: i.e. mean=0
    #then cross_correlation = pearson_correlation
    return N.sum(x*y)/N.sqrt(N.sum(N.square(x))*N.sum(N.square(y)))

## Where we calculate scores for clusters
def calculate_scores_from_clusters(self, d1, d2, gf, return_key=True): # doing single pair in each task
    # Subtomogram 1
    vm = self.cache.get_mrc(d1['mask'])
    v  = self.cache.get_mrc(d1['subtomogram'])
    ang = N.array(d1["angle"], dtype=N.float)
    loc = N.array(d1["loc"], dtype=N.float)
    vr_1 = None
    vm_1 = None
    vr_1 = GR.rotate_pad_mean(v, angle=ang, loc_r=loc);           assert N.all(N.isfinite(vr_1))
    vm_1 = GR.rotate_mask(vm, angle=ang);           assert N.all(N.isfinite(vm_1))

    # Subtomogram 2
    vm = self.cache.get_mrc(d2['mask'])
    v  = self.cache.get_mrc(d2['subtomogram'])
    ang = N.array(d2["angle"], dtype=N.float)
    loc = N.array(d2["loc"], dtype=N.float)
    vr_2 = None
    vm_2 = None
    vr_2 = GR.rotate_pad_mean(v, angle=ang, loc_r=loc);           assert N.all(N.isfinite(vr_2))
    vm_2 = GR.rotate_mask(vm, angle=ang);           assert N.all(N.isfinite(vm_2))

    mask_cutoff = 0.5
    vm_1[vm_1 < mask_cutoff] = 0.0
    vm_1[vm_1 >= mask_cutoff] = 1.0
    vm_2[vm_2 < mask_cutoff] = 0.0
    vm_2[vm_2 >= mask_cutoff] = 1.0
    masks_logical_and = N.logical_and(vm_1, vm_2)
    masks_logical_and_flag = False
    if masks_logical_and.sum() < 2:
        masks_logical_and_flag = True
    else:
        masks_logical_and = masks_logical_and.flatten()
        masks_logical_and = N.where(masks_logical_and==True)[0]
    
    # Now we have rotated subtomograms. So calculate all the scores.
    scores = {}
    scores['pearson_correlation'] = {}
    scores['estimated-SNR'] = {}
    scores['fourier_pcorr'] = {}
    scores['fourier_pcorr_mw'] = {}
    scores['constrained_pcorr_on_normalised'] = {}
    scores['constrained_pcorr_on_gaussian_filtered'] = {}
    scores['overlap_score'] = {}
    scores['overlap_pcorr_on_normalised'] = {}
    scores['overlap_pcorr_on_gaussian_filtered'] = {}
    scores['lsf'] = {}
    scores['constained_lsf'] = {}
    scores['overlap_lsf'] = {}
    scores['dlsf'] = {}
    scores['mutual_information'] = {}
    scores['constrained_mutual_information'] = {}
    scores['overlap_mutual_information'] = {}
    scores['mutual_information_norm'] = {}
    scores['constrained_mutual_information_norm'] = {}
    scores['overlap_mutual_information_norm'] = {}

    # # MinMax_Normalization
    # vr_1_norm = min_max(vr_1)
    # vr_2_norm = min_max(vr_2)

    # ZeroMean UnitVariance Normalization
    vr_1_norm = normalize(vr_1)
    vr_2_norm = normalize(vr_2)

    for gaussian_filter in range(gf):
        # Gaussian Filter
        vr_1_gf = SNFG(vr_1_norm, gaussian_filter)
        vr_2_gf = SNFG(vr_2_norm, gaussian_filter)
        
        # Create Masks
        threshold_i = 1.5

        # # When using min_max normalization
        # vr_1_mask = mask_min_max(vr_1_gf, threshold_i)
        # vr_2_mask = mask_min_max(vr_2_gf, threshold_i)

        # When using normal distribution normalization
        vr_1_mask = mask_normalize(vr_1_gf, threshold_i)
        vr_2_mask = mask_normalize(vr_2_gf, threshold_i)
        
        # Pearson Correlation
        pcorr = pearsonr(vr_1_gf.flatten(), vr_2_gf.flatten())
        if math.isnan(pcorr[0]):
            scores['pearson_correlation'][str(gaussian_filter)] = 0.0
        else:
            scores['pearson_correlation'][str(gaussian_filter)] = pcorr[0]
        
        # Estimated-SNR
        if pcorr[0]==1:
            scores['estimated-SNR'][str(gaussian_filter)] = 1000.0
        elif math.isnan(pcorr[0]):
            scores['estimated-SNR'][str(gaussian_filter)] = 0.0
        else:
            scores['estimated-SNR'][str(gaussian_filter)] = pcorr[0]/(1-pcorr[0])

        # pcorr in fourier domain without missing wedge
        vr_1_f = (NF.fftshift(NF.fftn(vr_1_gf.copy()))).real
        vr_2_f = (NF.fftshift(NF.fftn(vr_2_gf.copy()))).real
        pcorr = pearsonr(vr_1_f.flatten(), vr_2_f.flatten())
        if math.isnan(pcorr[0]):
            scores['fourier_pcorr'][str(gaussian_filter)] = 0.0
        else:
            scores['fourier_pcorr'][str(gaussian_filter)] = pcorr[0]

        # pcorr in fourier domain with missing wedge
        if masks_logical_and_flag:
            scores['fourier_pcorr_mw'][str(gaussian_filter)] = 0.0
        else:
            pcorr = pearsonr(vr_1_f.flatten()[masks_logical_and], vr_2_f.flatten()[masks_logical_and])
            if math.isnan(pcorr[0]):
                scores['fourier_pcorr_mw'][str(gaussian_filter)] = 0.0
            else:
                scores['fourier_pcorr_mw'][str(gaussian_filter)] = pcorr[0]

        del vr_1_f, vr_2_f
        gc.collect()

        # Constrained pcorr on normalised: Use mask from guassian filtered map but calculate pcorr on original map
        # Constrained pcorr on normalised guassian_filtered: Use mask from guassian filtered map but calculate pcorr on guassian filtered again
        tmp_or = N.logical_or(vr_1_mask, vr_2_mask)
        if tmp_or.sum() < 2:
            scores['constrained_pcorr_on_normalised'][str(gaussian_filter)] = 0.0
            scores['constrained_pcorr_on_gaussian_filtered'][str(gaussian_filter)] = 0.0
        else:
            tmp_or = tmp_or.flatten()
            tmp_or = N.where(tmp_or==True)[0]
            # pcorr between original subtomograms in the region where any one of the subtomogram is 1 after thresholding
            constrained_pcorr = pearsonr(vr_1_norm.flatten()[tmp_or], vr_2_norm.flatten()[tmp_or])
            if math.isnan(constrained_pcorr[0]):
                scores['constrained_pcorr_on_normalised'][str(gaussian_filter)] = 0.0
            else:
                scores['constrained_pcorr_on_normalised'][str(gaussian_filter)] = constrained_pcorr[0]
            constrained_pcorr = pearsonr(vr_1_gf.flatten()[tmp_or], vr_2_gf.flatten()[tmp_or])
            if math.isnan(constrained_pcorr[0]):
                scores['constrained_pcorr_on_gaussian_filtered'][str(gaussian_filter)] = 0.0
            else:
                scores['constrained_pcorr_on_gaussian_filtered'][str(gaussian_filter)] = constrained_pcorr[0]
        
        # Overlap score
        tmp_and = N.logical_and(vr_1_mask, vr_2_mask)
        scores['overlap_score'][str(gaussian_filter)] = float(tmp_and.sum())/min(vr_1_mask.sum(), vr_2_mask.sum())
        
        # Overlap pcorr on normalised
        # Overlap pcorr on guassian filtered
        if tmp_and.sum() < 2:
            scores['overlap_pcorr_on_normalised'][str(gaussian_filter)] = 0.0
            scores['overlap_pcorr_on_gaussian_filtered'][str(gaussian_filter)] = 0.0
        else:
            tmp_and = tmp_and.flatten()
            tmp_and = N.where(tmp_and==True)[0]
            # Pearson correlation between original subtomograms in the region of overlap
            overlap_score_pcorr = pearsonr(vr_1_norm.flatten()[tmp_and], vr_2_norm.flatten()[tmp_and])
            if math.isnan(overlap_score_pcorr[0]):
                scores['overlap_pcorr_on_normalised'][str(gaussian_filter)] = 0.0
            else:
                scores['overlap_pcorr_on_normalised'][str(gaussian_filter)] = overlap_score_pcorr[0]
            overlap_score_pcorr = pearsonr(vr_1_gf.flatten()[tmp_and], vr_2_gf.flatten()[tmp_and])
            if math.isnan(overlap_score_pcorr[0]):
                scores['overlap_pcorr_on_gaussian_filtered'][str(gaussian_filter)] = 0.0
            else:
                scores['overlap_pcorr_on_gaussian_filtered'][str(gaussian_filter)] = overlap_score_pcorr[0]

        # LSF
        scores['lsf'][str(gaussian_filter)] = ((vr_1_gf - vr_2_gf)**2).mean()
        # constrained LSF
        scores['constained_lsf'][str(gaussian_filter)] = ((vr_1_gf.flatten()[tmp_or] - vr_2_gf.flatten()[tmp_or])**2).mean()
        # overlap LSF
        scores['overlap_lsf'][str(gaussian_filter)] = ((vr_1_gf.flatten()[tmp_and] - vr_2_gf.flatten()[tmp_and])**2).mean()
        
        # DLSF
        vr_1_sig_pairs = get_random_significant_pairs(vr_1_gf, 10000)
        tmp_score = 0.0
        for p in vr_1_sig_pairs:
            z1 = int(p[0])
            y1 = int(p[1])
            x1 = int(p[2])
            z2 = int(p[3])
            y2 = int(p[4])
            x2 = int(p[5])
            dens = p[6]
            prot_dens = vr_2_gf[z1][y1][x1] - vr_2_gf[z2][y2][x2]
            tmp_score += (dens-prot_dens)**2
        scores['dlsf'][str(gaussian_filter)] = tmp_score/vr_1_gf.size
        del vr_1_sig_pairs, tmp_score, prot_dens, z1, z2, y1, y2, x1, x2, dens
        gc.collect()
        
        # Mutual Information
        tmp_and = N.logical_and(vr_1_mask, vr_2_mask)
        tmp_or = N.logical_or(vr_1_mask, vr_2_mask)
        scores['mutual_information'][str(gaussian_filter)] = MI(vr_1_gf, vr_2_gf, mask_array=None, normalised=False)
        scores['constrained_mutual_information'][str(gaussian_filter)] = MI(vr_1_gf, vr_2_gf, mask_array=tmp_or, normalised=False)
        scores['overlap_mutual_information'][str(gaussian_filter)] = MI(vr_1_gf, vr_2_gf, mask_array=tmp_and, normalised=False)
        scores['mutual_information_norm'][str(gaussian_filter)] = MI(vr_1_gf, vr_2_gf, mask_array=None, normalised=True)
        scores['constrained_mutual_information_norm'][str(gaussian_filter)] = MI(vr_1_gf, vr_2_gf, mask_array=tmp_or, normalised=True)
        scores['overlap_mutual_information_norm'][str(gaussian_filter)] = MI(vr_1_gf, vr_2_gf, mask_array=tmp_and, normalised=True)

    del vr_1, vm_1, vr_2, vm_2
    gc.collect()

    if return_key:
        re_key = self.cache.save_tmp_data(scores, fn_id=self.task.task_id)
        assert re_key is not None
        return {'key':re_key}
    else:
        return scores
