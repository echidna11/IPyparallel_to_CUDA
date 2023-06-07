import json, os, pickle, psutil, copy, uuid
import sys, shutil, math, random, gc, time
import numpy as N
import numpy.fft as NF
import tomominer.io.file as IF
from scipy.ndimage import affine_transform
import tomominer.image.vol.util as GV
from scipy.stats import pearsonr
from scipy.ndimage.filters import gaussian_filter as SNFG
import warnings, math
# from random import randrange
# warnings.filterwarnings("ignore")

def eulerAnglesToRotationMatrix(theta):
    R_x = N.array([[1,         0,                  0                   ],
                    [0,         math.cos(theta[0]), -math.sin(theta[0]) ],
                    [0,         math.sin(theta[0]), math.cos(theta[0])  ]])

    R_y = N.array([[math.cos(theta[1]),    0,      math.sin(theta[1])  ],
                    [0,                     1,      0                   ],
                    [-math.sin(theta[1]),   0,      math.cos(theta[1])  ]])

    R_z = N.array([[math.cos(theta[2]),    -math.sin(theta[2]),    0],
                    [math.sin(theta[2]),    math.cos(theta[2]),     0],
                    [0,                     0,                      1]])

    R = N.dot(R_z, N.dot( R_y, R_x ))
    del R_x, R_y, R_z
    return R

def var__local(self, data_json, gaussian_filter, labels=None, return_key=True, mask_cutoff=0.5):
    if labels is None:        labels = [0] * len(data_json)
    sum_v = {}
    prod_sum_v = {}
    mask_sum = {}
    for i, r in enumerate(data_json):
        m = IF.get_mrc(r['mask'])
        size_m = N.array(m.shape, dtype = N.float)
        c_m = (size_m -1 ) / 2.0
        v = IF.get_mrc(r['subtomogram'])
        size_v = N.array(v.shape, dtype = N.float)
        c_v = (size_v -1 ) / 2.0
        angle = N.array(r['model']['angle'])
        misalign_angle =  N.array(r['model']['misalign_angle'])
        rm1 = eulerAnglesToRotationMatrix(angle)
        rm1 = rm1.transpose() # reverse the rotation angles
        rm2 = eulerAnglesToRotationMatrix(misalign_angle)
        # We need to rotate back the subtomogram to original position because angle was used to rotate map for simulation. Also we need to rotate further for misalignment. Similar with mask.
        # Rotation matrix for both map and mask is same
        rm = N.matmul(rm1, rm2)
        # Rotate subtomogram
        c = -rm.dot(c_v) + c_v
        vr = affine_transform(input=v, matrix=rm, offset=c, output_shape=size_v.astype(N.int), cval=float('NaN'))
        vr[N.logical_not(N.isfinite(vr))] = v.mean()
        # Rotate mask
        c = -rm.dot(c_m) + c_m
        vm = affine_transform(input=m, matrix=rm, offset=c, output_shape=size_m.astype(N.int), cval=float('NaN'))
        vm[N.logical_not(N.isfinite(vm))] = 0.0

        
        if gaussian_filter != 0:
            vr = SNFG(vr, gaussian_filter)
        
        vr = NF.fftshift( NF.fftn(vr) )

        vr[vm < mask_cutoff] = 0.0
        if labels[i] not in sum_v:
            sum_v[labels[i]] = vr
        else:
            sum_v[labels[i]] += vr

        if labels[i] not in prod_sum_v:
            prod_sum_v[labels[i]] = vr * N.conj(vr)
        else:
            prod_sum_v[labels[i]] += vr * N.conj(vr)

        if labels[i] not in mask_sum:
            mask_sum[labels[i]] = N.zeros(vm.shape, dtype=N.int)
        
        mask_sum[labels[i]][vm >= mask_cutoff] += 1
        del m, v, rm1, rm2, rm, vr, vm, size_m, c_m, size_v, c_v, c, angle, misalign_angle
        gc.collect()
    re = {'sum':sum_v, 'prod_sum':prod_sum_v, 'mask_sum':mask_sum}
    if return_key:
        re_key = self.cache.save_tmp_data({}, fn_id=self.task.task_id)
        assert re_key is not None
        return {'key':re_key}
    else:
        return re
        
def var__global(self, clusters, op):
    data_json_copy = []
    labels_copy = []
    for l in clusters:
        for d in clusters[l]:
            data_json_copy.append(d)
            labels_copy.append(l)

    tasks = []
    n_chunk = int( N.ceil(len(data_json_copy) / self.runner.work_queue.get_worker_number()))
    n_chunk = N.max((n_chunk, 10))
    while data_json_copy:
        data_json_copy_part = data_json_copy[:n_chunk]
        labels_copy_part = labels_copy[:n_chunk]

        tasks.append(self.runner.task(module='tomominer.simulation.calculate_scores', method='var__local', kwargs={'data_json':data_json_copy_part, 'labels':labels_copy_part, 'gaussian_filter': op['gaussian_filter'] if 'gaussian_filter' in op else 0}))

        data_json_copy = data_json_copy[n_chunk:]
        labels_copy = labels_copy[n_chunk:]


    del data_json_copy
    del labels_copy
    gc.collect()
    sum_global = {}
    prod_sum = {}
    mask_sum = {}
    for res in self.runner.run__except(tasks):
        with open(res.result['key']) as f:     re = pickle.load(f)
        os.remove(res.result['key'])
        for l in re['sum']:
            if l not in sum_global:
                sum_global[l] = re['sum'][l]
            else:
                sum_global[l] += re['sum'][l]
            
            if l not in prod_sum:
                prod_sum[l] = re['prod_sum'][l]   
            else:
                prod_sum[l] += re['prod_sum'][l]

            if l not in mask_sum:
                mask_sum[l] = re['mask_sum'][l]
            else:
                mask_sum[l] += re['mask_sum'][l]
    return {'sum':sum_global, 'prod_sum':prod_sum, 'mask_sum':mask_sum}

# get a volume containing radius
def ssnr__get_rad(siz):
    grid = GV.grid_displacement_to_center(siz, GV.fft_mid_co(siz))
    rad = GV.grid_distance_to_center(grid)
    return rad

def ssnr_to_fsc(ssnr):
    fsc = ssnr / (2.0 + ssnr)
    return fsc

# get index within certain frequency band
def ssnr_rad_ind(rad, r, band_width_radius):
    ###
    return ( abs(rad - r) <= band_width_radius )

def ssnr__given_stat(sum_v, prod_sum, mask_sum, rad=None, op=None):
    op = copy.deepcopy(op)

    if 'band_width_radius' not in op:
        op['band_width_radius'] = 1.0

    if 'mask_sum_threshold' not in op:
        op['mask_sum_threshold'] = 2.0
    else:
        op['mask_sum_threshold'] = N.max((op['mask_sum_threshold'], 2.0))

    siz = N.array(sum_v.shape)
    subtomogram_num = mask_sum.max()
    avg = N.zeros(sum_v.shape, dtype=N.complex) + N.nan
    ind = mask_sum > 0
    avg[ind] = sum_v[ind]  /  mask_sum[ind];    avg_abs_sq = N.real(  avg * N.conj( avg ) )
    del ind

    var = N.zeros(sum_v.shape, dtype=N.complex) + N.nan
    ind = mask_sum >= op['mask_sum_threshold']
    var[ind] = ( prod_sum[ind] - mask_sum[ind]*(avg[ind]*N.conj(avg[ind])) ) / ( mask_sum[ind] - 1 );   var = N.real(var)
    del ind

    if rad is None:     rad = ssnr__get_rad(siz)
    vol_rad = int( N.floor( N.min(siz) / 2.0 ) + 1)
    ssnr = N.zeros(vol_rad) + N.nan     # this is the SSNR of the AVERAGE image

    # the interpolation can also be performed using scipy.ndimage.interpolation.map_coordinates()
    for r in range(vol_rad):
        ind = ssnr_rad_ind(rad=rad, r=r, band_width_radius=op['band_width_radius']) # in order to use it as an index or mask, must convert to a bool array, not integer array!
        ind[mask_sum < op['mask_sum_threshold']] = False
        ind[N.logical_not(N.isfinite(avg))] = False
        ind[N.logical_not(N.isfinite(var))] = False
        if op['method'] == 1:
            if var[ind].sum() > 0:
                ssnr[r] = (mask_sum[ind] * avg_abs_sq[ind]).sum() / var[ind].sum()
            else:
                ssnr[r] = 0.0
        del ind

    assert N.all(N.isfinite(ssnr))
    fsc = ssnr_to_fsc(ssnr)
    return {'ssnr':ssnr, 'fsc':fsc}

def ssnr_parallel(self, clusters, op):
    c = var__global(self=self, clusters=clusters, op=op)
    fs = {}
    if 'ssnr' not in fs:
        fs['ssnr'] = {}
    if 'fsc' not in fs:
        fs['fsc'] = {}
    if 'fsc_sum' not in fs:
        fs['fsc_sum'] = {}
    for l in c['sum']:
        ssnr_re = ssnr__given_stat(sum_v=c['sum'][l], prod_sum=c['prod_sum'][l], mask_sum=c['mask_sum'][l], op=op['ssnr'])
        if l not in fs['ssnr']:
            fs['ssnr'][l] = ssnr_re['ssnr'].tolist()
        if l not in fs['fsc']:
            fs['fsc'][l] = ssnr_re['fsc'].tolist()
        if l not in fs['fsc_sum']:
            fs['fsc_sum'][l] = sum(fs['fsc'][l])
    
    del c
    gc.collect()
    return fs

def calculate_scores(self, d1, d2, misalign_angles, return_key=True):
    time_1 = time.time()
    # Subtomogram 1
    m = IF.get_mrc(d1['mask'])
    size_m = N.array(m.shape, dtype = N.float)
    c_m = (size_m -1 ) / 2.0
    v = IF.get_mrc(d1['subtomogram'])
    size_v = N.array(v.shape, dtype = N.float)
    c_v = (size_v -1 ) / 2.0
    angle = N.array(d1['model']['angle'])
    rm = eulerAnglesToRotationMatrix(angle)
    rm = rm.transpose() # reverse the rotation angles
    # We need to rotate back the subtomogram to original position because angle was used to rotate map for simulation.
    # Rotation matrix for both map and mask is same
    # Rotate subtomogram
    c = -rm.dot(c_v) + c_v
    vr_1 = affine_transform(input=v, matrix=rm, offset=c, output_shape=size_v.astype(N.int), cval=0.0)
    # Rotate mask
    c = -rm.dot(c_m) + c_m
    vm_1 = affine_transform(input=m, matrix=rm, offset=c, output_shape=size_m.astype(N.int), cval=0.0)

    # Subtomogram 2
    m = IF.get_mrc(d2['mask'])
    size_m = N.array(m.shape, dtype = N.float)
    c_m = (size_m -1 ) / 2.0
    v = IF.get_mrc(d2['subtomogram'])
    size_v = N.array(v.shape, dtype = N.float)
    c_v = (size_v -1 ) / 2.0
    angle = N.array(d2['model']['angle'])
    misalign_angle =  N.array(misalign_angles)
    rm1 = eulerAnglesToRotationMatrix(angle)
    rm1 = rm1.transpose() # reverse the rotation angles
    rm2 = eulerAnglesToRotationMatrix(misalign_angle)
    # We need to rotate back the subtomogram to original position because angle was used to rotate map for simulation. Also we need to rotate further for misalignment. Similar with mask.
    # Rotation matrix for both map and mask is same
    rm = N.matmul(rm1, rm2)
    # Rotate subtomogram
    c = -rm.dot(c_v) + c_v
    vr_2 = affine_transform(input=v, matrix=rm, offset=c, output_shape=size_v.astype(N.int), cval=0.0)
    # Rotate mask
    c = -rm.dot(c_m) + c_m
    vm_2 = affine_transform(input=m, matrix=rm, offset=c, output_shape=size_m.astype(N.int), cval=0.0)
    time_2 = time.time()
    # Now we have rotated subtomograms. So calculate all the scores.
    scores = {}
    pcorr = pearsonr(vr_1.flatten(), vr_2.flatten())
    scores['pearson_correlation'] = pcorr
    scores['estimated-SNR'] = pcorr[0]/(1-pcorr[0])
    # Guassian filters
    filter_sigmas = range(6)[1:] # 1 to 5
    for filter_sigma in filter_sigmas:
        vr_1_g = SNFG(vr_1, filter_sigma)
        vr_2_g = SNFG(vr_2, filter_sigma)
        pcorr = pearsonr(vr_1_g.flatten(), vr_2_g.flatten())
        score_name = "pearson_correlation_gaussian_filtered_sigma_"+str(filter_sigma)
        scores[score_name] = pcorr
        #scores['estimated-SNR'] = pcorr[0]/(1-pcorr[0])
    
    #scores['misalign_angles'] = misalign_angles
    #scores['sf1'] = sf1(vr_1, vr_2)
    #fs[misalign]['scores'] = ssnr_parallel(clusters=clusters, op=op, fs=fs[misalign]['scores'], rep=rep)

    del m, size_m, c_m, v, size_v, c_v, angle, misalign_angle, rm1, rm2, rm, c, vr_1, vr_2, vm_1, vm_2, pcorr, filter_sigma, filter_sigmas, vr_1_g, vr_2_g, score_name
    gc.collect()

    if return_key:
        re_key = self.cache.save_tmp_data(scores, fn_id=self.task.task_id)
        assert re_key is not None
        return {'key':re_key}
    else:
        return scores

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

'''
def MI(v1, v2, layers1=20, layers2=20, mask_array=None, normalised=True):
    # based on mask_array MI calculated on complete map (All 1 mask), contoured map (OR on masks), Overlap region (AND on masks)
    
    if mask_array is None:
        mask_array = N.ones(v1.shape, dtype=int)
    else:
        mask_sum = N.sum(mask_array)
        if mask_sum == 0:
            return 0.0

    # for us the formule for number of layers gives the same value 20, as default value. :)
    # sturges rule provides a way of calculating number of bins : 1+math.log(number of points)
    layers1 = int(1 + math.log(v1.size, 2))
    layers2 = int(1 + math.log(v2.size, 2))
    layers1 = max(layers1,8)
    layers2 = max(layers2,8)
    
    # digitize masked map based on layers
    map1_bin = v1.copy()
    map2_bin = v2.copy()
    map1_bin = map1_bin*mask_array
    map2_bin = map2_bin*mask_array
    map1_bin = digitize_map(v=map1_bin, cutoff=v1[mask_array].min(), nbins=layers1, left=True)
    map2_bin = digitize_map(v=map2_bin, cutoff=v2[mask_array].min(), nbins=layers2, left=True)
    
    # make sure the outside region is filled with zeros
    map1_bin = map1_bin*mask_array
    map2_bin = map2_bin*mask_array
    
    # background frequencies from the whole map
    bins1 = []
    for i in range(layers1+2): bins1.append(i)
    bins2 = []
    for i in range(layers2+2): bins2.append(i)
    
    # calculate frequency of bins
    map1_freq = N.histogram(map1_bin, bins1)[0][1:]
    map2_freq = N.histogram(map2_bin, bins2)[0][1:]
    
    if N.sum(map1_freq) == 0:
        #print 'No map overlap (Mutual information score), exiting score calculation..'
        return 0.0
    if N.sum(map2_freq) == 0:
        #print 'No map overlap (Mutual information score), exiting score calculation..'
        return 0.0
    
    list_overlaps = []
    total = 0
    for x in range(layers1):
        mask_array = map1_bin == float(x+1)
        overlap_freq =  N.histogram(map2_bin[mask_array], bins2)[0][1:]
        total += float(N.sum(overlap_freq))
        list_overlaps.append(overlap_freq)

    if total == 0:
        #print 'No map overlap (Mutual information score), exiting score calculation..'
        return 0.0
    
    enter = 0
    Hxy = 0.0
    Hx = 0.0
    Hy = 0.0
    for x in range(layers1):
        # probability of occurrence of x
        p_m1 = map1_freq[x]/float(N.sum(map1_freq))
        if not p_m1 == 0.0:
            Hx += (-p_m1*math.log(p_m1, 2))
        for y in range(layers2):
            enter = 1
            # probability for overlap of bins x and y
            p_comb = list_overlaps[x][y]/total
            if not p_comb == 0:
                Hxy += -p_comb*math.log(p_comb, 2) # joined entropy

            if x == 0:
                # probability of occurrence of y
                p_m2 = map2_freq[y]/float(N.sum(map2_freq))
                if not p_m2 == 0.0:
                    Hy += (-p_m2*math.log(p_m2, 2))
    if enter == 1:
        # normalised MI (Studholme et al.) is used to account for overlap of 'contours'
        if normalised:
            if Hxy == 0.0: return 0.0
            return (Hx + Hy)/Hxy
        return Hx + Hy - Hxy

def digitize_map(v, cutoff, nbins, left=True):
    binMap = v.copy()
    bins = []
    step = (v.max()-float(cutoff))/nbins
    ini = float(cutoff) + (0.000001*step)
    if left:
        ini = float(cutoff) - (0.000001*step)
    bins.append(ini)
    for ii in range(1,nbins+1):
        bins.append(float(cutoff) + ii*step)
    if bins[-1] < v.max():
        bins = bins[:-1]
        bins.append(v.max())

    for z in range(len(v)):
        for y in range(len(v[z])):
            binMap[z][y] = N.digitize(v[z][y], bins)
    return binMap
'''

def cross_correlation(x, y):
    #if x and y are normalised: i.e. mean = 0 and std = 1 or x and y are centered: i.e. mean=0
    #then cross_correlation = pearson_correlation
    return N.sum(x*y)/N.sqrt(N.sum(N.square(x))*N.sum(N.square(y)))

## Where we calculate scores for clusters
def calculate_scores_from_clusters(self, d1, d2, gf, return_key=True): # doing single pair in each task
    # Subtomogram 1
    m = IF.get_mrc(d1['mask'])
    size_m = N.array(m.shape, dtype = N.float)
    c_m = (size_m -1 ) / 2.0
    v = IF.get_mrc(d1['subtomogram'])
    size_v = N.array(v.shape, dtype = N.float)
    c_v = (size_v -1 ) / 2.0
    angle = N.array(d1['model']['angle'])
    misalign_angle =  N.array(d1['model']['misalign_angle'])
    rm1 = eulerAnglesToRotationMatrix(angle)
    rm1 = rm1.transpose() # reverse the rotation angles
    rm2 = eulerAnglesToRotationMatrix(misalign_angle)
    # We need to rotate back the subtomogram to original position because angle was used to rotate map for simulation.
    # Rotation matrix for both map and mask is same
    rm = N.matmul(rm1, rm2)
    # Rotate subtomogram
    c = -rm.dot(c_v) + c_v
    vr_1 = SN.interpolation.affine_transform(input=v, matrix=rm, offset=c, output_shape=size_v.astype(N.int), cval=float('NaN'))
    vr_1[N.logical_not(N.isfinite(vr_1))] = v.mean()
    # Rotate mask
    c = -rm.dot(c_m) + c_m
    vm_1 = SN.interpolation.affine_transform(input=m, matrix=rm, offset=c, output_shape=size_m.astype(N.int), cval=float('NaN'))
    vm_1[N.logical_not(N.isfinite(vm_1))] = 0.0

    # Subtomogram 2
    m = IF.get_mrc(d2['mask'])
    size_m = N.array(m.shape, dtype = N.float)
    c_m = (size_m -1 ) / 2.0
    v = IF.get_mrc(d2['subtomogram'])
    size_v = N.array(v.shape, dtype = N.float)
    c_v = (size_v -1 ) / 2.0
    angle = N.array(d2['model']['angle'])
    misalign_angle =  N.array(d2['model']['misalign_angle'])
    rm1 = eulerAnglesToRotationMatrix(angle)
    rm1 = rm1.transpose() # reverse the rotation angles
    rm2 = eulerAnglesToRotationMatrix(misalign_angle)
    # We need to rotate back the subtomogram to original position because angle was used to rotate map for simulation. Also we need to rotate further for misalignment. Similar with mask.
    # Rotation matrix for both map and mask is same
    rm = N.matmul(rm1, rm2)
    # Rotate subtomogram
    c = -rm.dot(c_v) + c_v
    vr_2 = SN.interpolation.affine_transform(input=v, matrix=rm, offset=c, output_shape=size_v.astype(N.int), cval=float('NaN'))
    vr_2[N.logical_not(N.isfinite(vr_2))] = v.mean()
    # Rotate mask
    c = -rm.dot(c_m) + c_m
    vm_2 = SN.interpolation.affine_transform(input=m, matrix=rm, offset=c, output_shape=size_m.astype(N.int), cval=float('NaN'))
    vm_2[N.logical_not(N.isfinite(vm_2))] = 0.0

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

    del v, size_v, c_v, angle, misalign_angle, rm1, rm2, rm, c, vr_1, vr_2
    gc.collect()
    
    if return_key:
        re_key = self.cache.save_tmp_data(scores, fn_id=self.task.task_id)
        assert re_key is not None
        return {'key':re_key}
    else:
        return scores

'''
def calculate_scores_from_clusters_1(self, pairs, cluster, gf, return_key=True): # doing few pairs in single task
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
    
    for gaussian_filter in range(gf):
        for key in scores:
            scores[key][str(gaussian_filter)] = []

    for i in range(len(pairs)):
        d1 = cluster[pairs[i][0]]
        d2 = cluster[pairs[i][1]]
        
        # Subtomogram 1
        m = IF.get_mrc(d1['mask'])
        size_m = N.array(m.shape, dtype = N.float)
        c_m = (size_m -1 ) / 2.0
        v = IF.get_mrc(d1['subtomogram'])
        size_v = N.array(v.shape, dtype = N.float)
        c_v = (size_v -1 ) / 2.0
        angle = N.array(d1['model']['angle'])
        misalign_angle =  N.array(d1['model']['misalign_angle'])
        rm1 = eulerAnglesToRotationMatrix(angle)
        rm1 = rm1.transpose() # reverse the rotation angles
        rm2 = eulerAnglesToRotationMatrix(misalign_angle)
        # We need to rotate back the subtomogram to original position because angle was used to rotate map for simulation.
        # Rotation matrix for both map and mask is same
        rm = N.matmul(rm1, rm2)
        # Rotate subtomogram
        c = -rm.dot(c_v) + c_v
        vr_1 = affine_transform(input=v, matrix=rm, offset=c, output_shape=size_v.astype(N.int), cval=float('NaN'))
        vr_1[N.logical_not(N.isfinite(vr_1))] = v.mean()
        # Rotate mask
        c = -rm.dot(c_m) + c_m
        vm_1 = affine_transform(input=m, matrix=rm, offset=c, output_shape=size_m.astype(N.int), cval=float('NaN'))
        vm_1[N.logical_not(N.isfinite(vm_1))] = 0.0

        # Subtomogram 2
        m = IF.get_mrc(d2['mask'])
        size_m = N.array(m.shape, dtype = N.float)
        c_m = (size_m -1 ) / 2.0
        v = IF.get_mrc(d2['subtomogram'])
        size_v = N.array(v.shape, dtype = N.float)
        c_v = (size_v -1 ) / 2.0
        angle = N.array(d2['model']['angle'])
        misalign_angle =  N.array(d2['model']['misalign_angle'])
        rm1 = eulerAnglesToRotationMatrix(angle)
        rm1 = rm1.transpose() # reverse the rotation angles
        rm2 = eulerAnglesToRotationMatrix(misalign_angle)
        # We need to rotate back the subtomogram to original position because angle was used to rotate map for simulation. Also we need to rotate further for misalignment. Similar with mask.
        # Rotation matrix for both map and mask is same
        rm = N.matmul(rm1, rm2)
        # Rotate subtomogram
        c = -rm.dot(c_v) + c_v
        vr_2 = affine_transform(input=v, matrix=rm, offset=c, output_shape=size_v.astype(N.int), cval=float('NaN'))
        vr_2[N.logical_not(N.isfinite(vr_2))] = v.mean()
        # Rotate mask
        c = -rm.dot(c_m) + c_m
        vm_2 = affine_transform(input=m, matrix=rm, offset=c, output_shape=size_m.astype(N.int), cval=float('NaN'))
        vm_2[N.logical_not(N.isfinite(vm_2))] = 0.0

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
                scores['pearson_correlation'][str(gaussian_filter)].append(0.0)
            else:
                scores['pearson_correlation'][str(gaussian_filter)].append(pcorr[0])
            
            # Estimated-SNR
            if pcorr[0]==1:
                scores['estimated-SNR'][str(gaussian_filter)].append(1000.0)
            elif math.isnan(pcorr[0]):
                scores['estimated-SNR'][str(gaussian_filter)].append(0.0)
            else:
                scores['estimated-SNR'][str(gaussian_filter)].append(pcorr[0]/(1-pcorr[0]))

            # pcorr in fourier domain without missing wedge
            vr_1_f = (NF.fftshift(NF.fftn(vr_1_gf.copy()))).real
            vr_2_f = (NF.fftshift(NF.fftn(vr_2_gf.copy()))).real
            pcorr = pearsonr(vr_1_f.flatten(), vr_2_f.flatten())
            if math.isnan(pcorr[0]):
                scores['fourier_pcorr'][str(gaussian_filter)].append(0.0)
            else:
                scores['fourier_pcorr'][str(gaussian_filter)].append(pcorr[0])

            # pcorr in fourier domain with missing wedge
            if masks_logical_and_flag:
                scores['fourier_pcorr_mw'][str(gaussian_filter)].append(0.0)
            else:
                pcorr = pearsonr(vr_1_f.flatten()[masks_logical_and], vr_2_f.flatten()[masks_logical_and])
                if math.isnan(pcorr[0]):
                    scores['fourier_pcorr_mw'][str(gaussian_filter)].append(0.0)
                else:
                    scores['fourier_pcorr_mw'][str(gaussian_filter)].append(pcorr[0])

            del vr_1_f, vr_2_f
            gc.collect()

            # Constrained pcorr on normalised: Use mask from guassian filtered map but calculate pcorr on original map
            # Constrained pcorr on normalised guassian_filtered: Use mask from guassian filtered map but calculate pcorr on guassian filtered again
            tmp_or = N.logical_or(vr_1_mask, vr_2_mask)
            if tmp_or.sum() < 2:
                scores['constrained_pcorr_on_normalised'][str(gaussian_filter)].append(0.0)
                scores['constrained_pcorr_on_gaussian_filtered'][str(gaussian_filter)].append(0.0)
            else:
                tmp_or = tmp_or.flatten()
                tmp_or = N.where(tmp_or==True)[0]
                # pcorr between original subtomograms in the region where any one of the subtomogram is 1 after thresholding
                constrained_pcorr = pearsonr(vr_1_norm.flatten()[tmp_or], vr_2_norm.flatten()[tmp_or])
                if math.isnan(constrained_pcorr[0]):
                    scores['constrained_pcorr_on_normalised'][str(gaussian_filter)].append(0.0)
                else:
                    scores['constrained_pcorr_on_normalised'][str(gaussian_filter)].append(constrained_pcorr[0])
                constrained_pcorr = pearsonr(vr_1_gf.flatten()[tmp_or], vr_2_gf.flatten()[tmp_or])
                if math.isnan(constrained_pcorr[0]):
                    scores['constrained_pcorr_on_gaussian_filtered'][str(gaussian_filter)].append(0.0)
                else:
                    scores['constrained_pcorr_on_gaussian_filtered'][str(gaussian_filter)].append(constrained_pcorr[0])
            
            # Overlap score
            tmp_and = N.logical_and(vr_1_mask, vr_2_mask)
            scores['overlap_score'][str(gaussian_filter)].append(float(tmp_and.sum())/min(vr_1_mask.sum(), vr_2_mask.sum()))
            
            # Overlap pcorr on normalised
            # Overlap pcorr on guassian filtered
            if tmp_and.sum() < 2:
                scores['overlap_pcorr_on_normalised'][str(gaussian_filter)].append(0.0)
                scores['overlap_pcorr_on_gaussian_filtered'][str(gaussian_filter)].append(0.0)
            else:
                tmp_and = tmp_and.flatten()
                tmp_and = N.where(tmp_and==True)[0]
                # Pearson correlation between original subtomograms in the region of overlap
                overlap_score_pcorr = pearsonr(vr_1_norm.flatten()[tmp_and], vr_2_norm.flatten()[tmp_and])
                if math.isnan(overlap_score_pcorr[0]):
                    scores['overlap_pcorr_on_normalised'][str(gaussian_filter)].append(0.0)
                else:
                    scores['overlap_pcorr_on_normalised'][str(gaussian_filter)].append(overlap_score_pcorr[0])
                overlap_score_pcorr = pearsonr(vr_1_gf.flatten()[tmp_and], vr_2_gf.flatten()[tmp_and])
                if math.isnan(overlap_score_pcorr[0]):
                    scores['overlap_pcorr_on_gaussian_filtered'][str(gaussian_filter)].append(0.0)
                else:
                    scores['overlap_pcorr_on_gaussian_filtered'][str(gaussian_filter)].append(overlap_score_pcorr[0])

            # LSF
            scores['lsf'][str(gaussian_filter)].append(((vr_1_gf - vr_2_gf)**2).mean())
            # constrained LSF
            scores['constained_lsf'][str(gaussian_filter)].append(((vr_1_gf.flatten()[tmp_or] - vr_2_gf.flatten()[tmp_or])**2).mean())
            # overlap LSF
            scores['overlap_lsf'][str(gaussian_filter)].append(((vr_1_gf.flatten()[tmp_and] - vr_2_gf.flatten()[tmp_and])**2).mean())
            
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
            scores['dlsf'][str(gaussian_filter)].append(tmp_score/vr_1_gf.size)
            del vr_1_sig_pairs, tmp_score, prot_dens, z1, z2, y1, y2, x1, x2, dens
            gc.collect()
            
            # Mutual Information
            tmp_and = N.logical_and(vr_1_mask, vr_2_mask)
            tmp_or = N.logical_or(vr_1_mask, vr_2_mask)
            scores['mutual_information'][str(gaussian_filter)].append(MI(vr_1_gf, vr_2_gf, mask_array=None, normalised=False))
            scores['constrained_mutual_information'][str(gaussian_filter)].append(MI(vr_1_gf, vr_2_gf, mask_array=tmp_or, normalised=False))
            scores['overlap_mutual_information'][str(gaussian_filter)].append(MI(vr_1_gf, vr_2_gf, mask_array=tmp_and, normalised=False))
            scores['mutual_information_norm'][str(gaussian_filter)].append(MI(vr_1_gf, vr_2_gf, mask_array=None, normalised=True))
            scores['constrained_mutual_information_norm'][str(gaussian_filter)].append(MI(vr_1_gf, vr_2_gf, mask_array=tmp_or, normalised=True))
            scores['overlap_mutual_information_norm'][str(gaussian_filter)].append(MI(vr_1_gf, vr_2_gf, mask_array=tmp_and, normalised=True))

        del v, size_v, c_v, angle, misalign_angle, rm1, rm2, rm, c, vr_1, vr_2
        gc.collect()

    if return_key:
        re_key = self.cache.save_tmp_data(scores, fn_id=self.task.task_id)
        assert re_key is not None
        return {'key':re_key}
    else:
        return scores
'''
'''
if True:
    for i in range(22):
        f1 = "../4/20/0.01/30/subtomograms/"+str(i)+"/0/0.mrc"
        v = IF.get_mrc(f1)
        for gf in range(3,4):
            vgf = SNFG(v, gf)
            vgf = min_max(vgf)
            m = N.mean(vgf)
            s = N.std(vgf)
            tmp = m-3*s
            vgf[vgf<=tmp] = -1.0
            vgf[vgf>tmp] = -0.0
            vgf = -vgf
            IF.put_mrc(vgf, str(i)+"_"+str(gf)+".mrc")

    import matplotlib.pyplot as plt
    with open("cluster_labels.json") as f:
        cluster_labels = json.load(f)
    cluster_labels = {cluster_labels[_]:_ for _ in cluster_labels.keys()}
    with open("../4/20/0.01/30/data_config.json") as f:
        dj = json.load(f)
    for i in range(22):
        print "ang\tp\tp_gf\tc_p\tos\tos_p"
        p = []
        p_gf = []
        os = []
        os_p = []
        c_p = []
        dj_tmp = random.sample([_ for _ in dj if _['cluster_label']==i], 2)
        for ang in range(20):
            # Rotate one subtomogram to original position
            v = IF.get_mrc(dj_tmp[0]['subtomogram'])
            size_v = N.array(v.shape, dtype = N.float)
            c_v = (size_v -1 ) / 2.0
            angle = N.array(dj_tmp[0]['model']['angle'])
            rm = eulerAnglesToRotationMatrix(angle)
            rm = rm.transpose()
            c = -rm.dot(c_v) + c_v
            vr_1 = SN.interpolation.affine_transform(input=v, matrix=rm, offset=c, output_shape=size_v.astype(N.int), cval=0.0)
            # Rotate another subtomogram to original position plus misangle
            v = IF.get_mrc(dj_tmp[1]['subtomogram'])
            size_v = N.array(v.shape, dtype = N.float)
            c_v = (size_v -1 ) / 2.0
            angle = N.array(dj_tmp[1]['model']['angle'])
            misalign_angle = N.array([ang*N.pi/100]*3)
            rm1 = eulerAnglesToRotationMatrix(angle)
            rm1 = rm1.transpose()
            rm2 = eulerAnglesToRotationMatrix(misalign_angle)
            rm = N.matmul(rm1, rm2)
            c = -rm.dot(c_v) + c_v
            vr_2 = SN.interpolation.affine_transform(input=v, matrix=rm, offset=c, output_shape=size_v.astype(N.int), cval=0.0)
            # pcorr between original subtomograms
            pcorr = pearsonr(vr_1.flatten(), vr_2.flatten())
            p.append(pcorr[0])
            if ang==0:
                    IF.put_mrc(-min_max(vr_1), str(i)+".mrc")
                    IF.put_mrc(-min_max(vr_2), str(i)+"_r.mrc")
            for gf in range(2,4):
                vr_1_gf = SNFG(vr_1, gf)
                vr_2_gf = SNFG(vr_2, gf)
                if ang==0 or ang==10:
                    IF.put_mrc(-min_max(vr_1_gf), str(i)+"_"+str(gf)+"_"+str(ang)+".mrc")
                    IF.put_mrc(-min_max(vr_2_gf), str(i)+"_"+str(gf)+"_"+str(ang)+"_r.mrc")
                # pcorr between gaussian filtered subtomograms
                pcorr_gf = pearsonr(vr_1_gf.flatten(), vr_2_gf.flatten())
                p_gf.append(pcorr_gf[0])
                
                vr_1_os = copy.deepcopy(vr_1_gf)
                vr_2_os = copy.deepcopy(vr_2_gf)
                vr_1_os = min_max(vr_1_os)
                vr_2_os = min_max(vr_2_os)
                vr_1_os = threshold(vr_1_os)
                vr_2_os = threshold(vr_2_os)
                if ang==0 or ang==10:
                    IF.put_mrc(vr_1_os, str(i)+"_"+str(gf)+"_"+str(ang)+"_os.mrc")
                    IF.put_mrc(vr_2_os, str(i)+"_"+str(gf)+"_"+str(ang)+"_os_r.mrc")
                
                tmp = N.logical_or(vr_1_os, vr_2_os)
                tmp = tmp.flatten()
                tmp = N.where(tmp==1)[0]
                # pcorr between original subtomograms in the region where any one of the subtomogram is 1 after thresholding
                constrained_pcorr = pearsonr(vr_1.flatten()[tmp], vr_2.flatten()[tmp])
                c_p.append(constrained_pcorr[0])

                tmp = N.logical_and(vr_1_os, vr_2_os)
                overlap_score = float(tmp.sum())/(vr_1_os.shape[0]*vr_1_os.shape[1]*vr_1_os.shape[2])
                os.append(overlap_score)
                tmp = tmp.flatten()
                tmp = N.where(tmp==1)[0]
                # Pearson correlation between original subtomograms in the region of overlap
                overlap_score_pcorr = pearsonr(vr_1.flatten()[tmp], vr_2.flatten()[tmp])
                os_p.append(overlap_score_pcorr[0])
                #print ang, "\t", pcorr[0], "\t", pcorr_gf[0], "\t", overlap_score, "\t", overlap_score_pcorr[0]
                print "%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f" % (ang, pcorr[0], pcorr_gf[0], constrained_pcorr[0], overlap_score, overlap_score_pcorr[0])

        plt.plot(p, label="pcorr", color="maroon", linewidth=2)
        plt.plot(p_gf, label="pcorr_gf", color="green", linewidth=2)
        plt.plot(c_p, label="constrained_pcorr", color="pink", linewidth=2)
        plt.plot(os, label="overlap_score", color="orange", linewidth=2)
        plt.plot(os_p, label="overlap_score_pcorr", color="blue", linewidth=2)
        plt.hlines(0, 0, len(p)-1, color='k')
        plt.legend()
        plt.xlim(0, len(p)-1)
        plt.title(cluster_labels[i])
        #plt.show()
        plt.savefig(cluster_labels[i]+".png")
        plt.close()

    for i in range(22):
        print "ang\tp\tp_gf\tc_p\tos\tos_p"
        p = []
        p_gf = []
        os = []
        os_p = []
        c_p = []
        for ang in range(20):
            f1 = "../4/20/0.01/30/subtomograms/"+str(i)+"/0/0.mrc"
            vr_1 = IF.get_mrc(f1)
            size_v = N.array(vr_1.shape, dtype = N.float)
            c_v = (size_v -1 ) / 2.0
            angle = N.array([ang*N.pi/100]*3)
            rm = eulerAnglesToRotationMatrix(angle)
            c = -rm.dot(c_v) + c_v
            vr_2 = SN.interpolation.affine_transform(input=vr_1, matrix=rm, offset=c, output_shape=size_v.astype(N.int), cval=0.0)
            # pcorr between original subtomograms
            pcorr = pearsonr(vr_1.flatten(), vr_2.flatten())
            p.append(pcorr[0])
            for gf in range(2,3):
                vr_1_gf = SNFG(vr_1, gf)
                vr_2_gf = SNFG(vr_2, gf)
                if ang==10:
                    IF.put_mrc(-min_max(vr_1_gf), str(i)+"_"+str(gf)+"_"+str(ang)+".mrc")
                    IF.put_mrc(-min_max(vr_2_gf), str(i)+"_"+str(gf)+"_"+str(ang)+"_r.mrc")
                # pcorr between gaussian filtered subtomograms
                pcorr_gf = pearsonr(vr_1_gf.flatten(), vr_2_gf.flatten())
                p_gf.append(pcorr_gf[0])
                
                vr_1_os = copy.deepcopy(vr_1_gf)
                vr_2_os = copy.deepcopy(vr_2_gf)
                vr_1_os = min_max(vr_1_os)
                vr_2_os = min_max(vr_2_os)
                vr_1_os = threshold(vr_1_os)
                vr_2_os = threshold(vr_2_os)
                if ang==10:
                    IF.put_mrc(vr_1_os, str(i)+"_"+str(gf)+"_"+str(ang)+"_os.mrc")
                    IF.put_mrc(vr_2_os, str(i)+"_"+str(gf)+"_"+str(ang)+"_os_r.mrc")
                
                tmp = N.logical_or(vr_1_os, vr_2_os)
                tmp = tmp.flatten()
                tmp = N.where(tmp==1)[0]
                # pcorr between original subtomograms in the region where any one of the subtomogram is 1 after thresholding
                constrained_pcorr = pearsonr(vr_1.flatten()[tmp], vr_2.flatten()[tmp])
                c_p.append(constrained_pcorr[0])

                tmp = N.logical_and(vr_1_os, vr_2_os)
                overlap_score = float(tmp.sum())/(vr_1_os.shape[0]*vr_1_os.shape[1]*vr_1_os.shape[2])
                os.append(overlap_score)
                tmp = tmp.flatten()
                tmp = N.where(tmp==1)[0]
                # Pearson correlation between original subtomograms in the region of overlap
                overlap_score_pcorr = pearsonr(vr_1.flatten()[tmp], vr_2.flatten()[tmp])
                os_p.append(overlap_score_pcorr[0])
                #print ang, "\t", pcorr[0], "\t", pcorr_gf[0], "\t", overlap_score, "\t", overlap_score_pcorr[0]
                print "%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f" % (ang, pcorr[0], pcorr_gf[0], constrained_pcorr[0], overlap_score, overlap_score_pcorr[0])

        plt.plot(p, label="pcorr", color="maroon", linewidth=2)
        plt.plot(p_gf, label="pcorr_gf", color="green", linewidth=2)
        plt.plot(c_p, label="constrained_pcorr", color="pink", linewidth=2)
        plt.plot(os, label="overlap_score", color="orange", linewidth=2)
        plt.plot(os_p, label="overlap_score_pcorr", color="blue", linewidth=2)
        plt.hlines(0, 0, len(p)-1, color='k')
        plt.legend()
        plt.xlim(0, len(p)-1)
        plt.title(cluster_labels[i])
        plt.savefig(cluster_labels[i]+".png")
        plt.close()

    ### Testing pcorr in fourier domain with and without missing wedge
    # Subtomogram 1
    m = IF.get_mrc(dj[0]['mask'])
    size_m = N.array(m.shape, dtype = N.float)
    c_m = (size_m -1 ) / 2.0
    v = IF.get_mrc(dj[0]['subtomogram'])
    size_v = N.array(v.shape, dtype = N.float)
    c_v = (size_v -1 ) / 2.0
    angle = N.array(dj[0]['model']['angle'])
    rm = eulerAnglesToRotationMatrix(angle)
    rm = rm.transpose() # reverse the rotation angles
    c = -rm.dot(c_v) + c_v
    vr_1 = SN.interpolation.affine_transform(input=v, matrix=rm, offset=c, output_shape=size_v.astype(N.int), cval=float('NaN'))
    vr_1[N.logical_not(N.isfinite(vr_1))] = v.mean()
    # Rotate mask
    c = -rm.dot(c_m) + c_m
    vm_1 = SN.interpolation.affine_transform(input=m, matrix=rm, offset=c, output_shape=size_m.astype(N.int), cval=float('NaN'))
    vm_1[N.logical_not(N.isfinite(vm_1))] = 0.0

    # Subtomogram 2
    m = IF.get_mrc(dj[1]['mask'])
    size_m = N.array(m.shape, dtype = N.float)
    c_m = (size_m -1 ) / 2.0
    v = IF.get_mrc(dj[1]['subtomogram'])
    size_v = N.array(v.shape, dtype = N.float)
    c_v = (size_v -1 ) / 2.0
    angle = N.array(dj[1]['model']['angle'])
    misalign_angle =  N.array([0.5*N.pi, 0.5*N.pi, 0.5*N.pi])
    rm1 = eulerAnglesToRotationMatrix(angle)
    rm1 = rm1.transpose() # reverse the rotation angles
    rm2 = eulerAnglesToRotationMatrix(misalign_angle)
    rm = N.matmul(rm1, rm2)
    c = -rm.dot(c_v) + c_v
    vr_2 = SN.interpolation.affine_transform(input=v, matrix=rm, offset=c, output_shape=size_v.astype(N.int), cval=float('NaN'))
    vr_2[N.logical_not(N.isfinite(vr_2))] = v.mean()
    # Rotate mask
    c = -rm.dot(c_m) + c_m
    vm_2 = SN.interpolation.affine_transform(input=m, matrix=rm, offset=c, output_shape=size_m.astype(N.int), cval=float('NaN'))
    vm_2[N.logical_not(N.isfinite(vm_2))] = 0.0

    IF.put_mrc(vr_1, "vr_1.mrc")
    IF.put_mrc(vr_2, "vr_2.mrc")
    IF.put_mrc(vm_1, "vm_1.mrc")
    IF.put_mrc(vm_2, "vm_2.mrc")

    gaussian_filter = 2
    # Gaussian Filter
    vr_1_gf = SNFG(vr_1, gaussian_filter)
    vr_2_gf = SNFG(vr_2, gaussian_filter)
    IF.put_mrc(vr_1_gf, "vr_1_gf.mrc")
    IF.put_mrc(vr_2_gf, "vr_2_gf.mrc")
    # Pearson Correlation
    pcorr = pearsonr(vr_1_gf.flatten(), vr_2_gf.flatten())
    print "pcorr:", pcorr[0]

    # Estimated-SNR
    print "estimated-SNR:", pcorr[0]/(1-pcorr[0])

    # pcorr in fourier domain without missing wedge
    vr_1_f = (NF.fftshift(NF.fftn(vr_1_gf))).real
    vr_2_f = (NF.fftshift(NF.fftn(vr_2_gf))).real
    IF.put_mrc(vr_1_f, "vr_1_f.mrc")
    IF.put_mrc(vr_2_f, "vr_2_f.mrc")
    pcorr = pearsonr(vr_1_f.flatten(), vr_2_f.flatten())
    print "pcorr_fourier:", pcorr[0]

    # pcorr in fourier domain with missing wedge
    mask_cutoff = 0.5
    vm_1[vm_1 < mask_cutoff] = 0.0
    vm_1[vm_1 >= mask_cutoff] = 1.0
    vm_2[vm_2 < mask_cutoff] = 0.0
    vm_2[vm_2 >= mask_cutoff] = 1.0
    tmp = N.logical_and(vm_1, vm_2)
    IF.put_mrc(vm_1, "vm_1_th.mrc")
    IF.put_mrc(vm_2, "vm_2_th.mrc")
    IF.put_mrc(tmp, "vm_AND.mrc")

    vr_1_f_th = copy.deepcopy(vr_1_f)
    vr_2_f_th = copy.deepcopy(vr_2_f)
    vr_1_f_th[~tmp] = 0.0
    vr_2_f_th[~tmp] = 0.0
    IF.put_mrc(vr_1_f_th, "vr_1_f_th.mrc")
    IF.put_mrc(vr_2_f_th, "vr_2_f_th.mrc")

    tmp = tmp.flatten()
    tmp = N.where(tmp==True)[0]
    pcorr = pearsonr(vr_1_f.flatten()[tmp], vr_2_f.flatten()[tmp])
    print "pcorr_fourier_mw:", pcorr[0]

    vr_1_os = copy.deepcopy(vr_1_gf)
    vr_2_os = copy.deepcopy(vr_2_gf)
    vr_1_os = min_max(vr_1_os)
    vr_2_os = min_max(vr_2_os)
    vr_1_os = threshold(vr_1_os)
    vr_2_os = threshold(vr_2_os)
    IF.put_mrc(vr_1_os, "vr_1_os.mrc")
    IF.put_mrc(vr_2_os, "vr_2_os.mrc")
'''