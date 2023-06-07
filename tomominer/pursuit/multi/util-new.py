'''
at this momemt, put all necessarily modified functions into this file
if the interpolation based classification is a success, after the code become mature, move re-organize all functions into their corresponding folders

~/ln/tomominer/tomominer/pursuit/multi/util.py
'''

import os
import sys
import copy
import time
import traceback
import warnings
import uuid
import cPickle as pickle
import shutil
from collections import defaultdict
import scipy.stats as SCS


#from multiprocessing.pool import ThreadPool as Pool
from multiprocessing.pool import Pool as Pool

import numpy as N
import numpy.fft as NF

from tomominer.io.cache import Cache
from tomominer.common.obj import Object

import tomominer.dimension_reduction.util as DU
import tomominer.image.vol.util as UV
import tomominer.io.file as IV
import tomominer.geometry.ang_loc as AAL
import tomominer.geometry.rotate as GR
import tomominer.align.fast.full as AFF
import tomominer.align.refine.gradient_refine as AFGF
import tomominer.align.util as AU
import tomominer.statistics.ssnr as SS
import tomominer.segmentation.watershed as SW
import tomominer.segmentation.active_contour.chan_vese.segment as SACS
import tomominer.filter.gaussian as FG
import tomominer.model.util as MU

import tomominer.core as core

# aligning two subtomograms without using wedge masks
def align_vols_no_wedge(v1, v2, op=None, logger=None):

    err = None
    try:
        #if logger is not None:      logger.debug("combined_search: v1_key %s, \n m1_key %s,\n v2_key %s,\n m2_key %s", v1_key, m1_key, v2_key, m2_key)
        re = AFF.align_vols_no_wedge(v1=v1, v2=v2, L=op['L'])
        ang = re['angle']
        loc = re['loc']
        score = re['cor']

  
    except Exception as e:      # mxu: just bypass alignment errors

        score = float('nan')
        loc = N.zeros(3)
        ang = N.random.random(3) * (N.pi * 2)      # mxu: randomly assign an angle

        err = e


    return {'angle':ang, 'loc':loc, 'score':score, 'err':err}



def align_vols_with_wedge(v1, m1, v2, m2, op=None, logger=None):

    err = None
    if 'fast_align_and_refine' in op:
        # fast rotation alignment plus parallel stocastic gradient search

        try:
            faar_op = copy.deepcopy(op['fast_align_and_refine'])
            faar_op['fast_max_l'] = op['L']
            faar_op.update({'v1':v1, 'm1':m1, 'v2':v2, 'm2':m2})
            re = AFGF.fast_align_and_refine(faar_op)
            ang = re['ang']
            loc = re['loc_r']
            score = re['cor']

        except Exception as e:      # mxu: just bypass alignment errors
            err = traceback.format_exc()

    else:

        # rotation serach based fast alignment

        #if logger is not None:      logger.debug("combined_search: v1_key %s, \n m1_key %s,\n v2_key %s,\n m2_key %s", v1_key, m1_key, v2_key, m2_key)
        try:
            re = AU.align_vols(v1=v1, m1=m1, v2=v2, m2=m2, L=op['L'])
            ang = re['angle']
            loc = re['loc']
            score = re['score']
        except Exception as e:      # mxu: just bypass alignment errors
            err = traceback.format_exc()


    if err is not None:
        score = float('nan')
        loc = N.zeros(3)
        ang = N.random.random(3) * (N.pi * 2)      # mxu: randomly assign an angle

    return {'angle':ang, 'loc':loc, 'score':score, 'err':err}



def align_keys(self, v1k, v2k, op, v2_out_key=None):

    if self is None:
        self = Object()
        self.cache = Cache()


    if op['with_missing_wedge']:
        v1  = self.cache.get_mrc(v1k['subtomogram'])

        if 'mask' in v1k:   
            m1 = self.cache.get_mrc(v1k['mask'])
        else:   
            m1 = None

        v2  = self.cache.get_mrc(v2k['subtomogram'])

        if 'mask' in v2k:   
            m2 = self.cache.get_mrc(v2k['mask'])
        else: 
            m2 = None

        re = align_vols_with_wedge(v1, m1, v2, m2, op=op)

    else:

        v1  = self.cache.get_mrc(v1k['subtomogram'])
        v2  = self.cache.get_mrc(v2k['subtomogram'])

        re = align_vols_no_wedge(v1=v1, v2=v2, op=op)


    re['v1_key'] = v1k
    re['v2_key'] = v2k
    re['v2_out_key'] = v2_out_key


    # save output to a file
    if v2_out_key is not None:

        v2r = GR.rotate_pad_mean(v2, angle=re['angle'], loc_r=re['loc'])
        IV.put_mrc(v2r, v2_out_key['subtomogram'])

        if op['with_missing_wedge']:
            m2r = GR.rotate_mask(m2, angle=re['angle'])
            IV.put_mrc(m2r, v2_out_key['mask'])

    re['op'] = op

    return re



def pairwise_align_keys__global(self, v1ks, v2ks, op):

    tasks = []
    for i1,v1k in enumerate(v1ks):
        for i2,v2k in enumerate(v2ks):
            op_t = copy.deepcopy(op)
            op_t['i1'] = i1
            op_t['i2'] = i2
            tasks.append(self.runner.task(module='tomominer.pursuit.multi.util', method='align_keys', kwargs={'v1k':v1k, 'v2k':v2k, 'op':op_t}))


    re = []
    for res in self.runner.run__except(tasks):
        re.append(res.result)
        

    return re

def pairwise_align_keys__multiprocessing(self, v1ks, v2ks, op):
    print 'pairwise_align_keys__multiprocessing()'

    pool = Pool()
    pre = []
    for i1,v1k in enumerate(v1ks):
        for i2,v2k in enumerate(v2ks):
            op_t = copy.deepcopy(op)
            op_t['i1'] = i1
            op_t['i2'] = i2

            pre.append( pool.apply_async(func=align_keys, kwds={'self':self, 'v1k':v1k, 'v2k':v2k, 'op':op_t}) )

    re = []
    for i, r in enumerate(pre):
        re.append( r.get(999999) )

        print '\r ', i, float(i) / len(pre),
        sys.stdout.flush()
    
    print 'done'

    return re


def rotate_key(self, vk, vk_out, angle, loc):
    if self is not None:
        v  = self.cache.get_mrc(vk['subtomogram'])
    else:        
        v = IV.get_mrc(vk['subtomogram'])

    vr = GR.rotate_pad_mean(v, angle=angle, loc_r=loc)
    IV.put_mrc(vr, vk_out['subtomogram'])

    if ('mask' in vk) and ('mask' in vk_out):
        if self is not None:        
            m = self.cache.get_mrc(vk['mask'])
        else:       
            m = IV.get_mrc(vk['mask'])

        mr = GR.rotate_mask(m, angle=angle)
        IV.put_mrc(mr, vk_out['mask'])

def copy_keys(vk, vk_out):
    shutil.copyfile(vk['subtomogram'], vk_out['subtomogram'])

    if ('mask' in vk) and ('mask' in vk_out):
        shutil.copyfile(vk['mask'], vk_out['mask'])



# supporse the subtomogram (v) and its mask (vm) is aligned with template t, fill up missing value regions of v using t
def impute_aligned_vols(t, v, vm, normalize=None):
    assert  normalize is not None

    if normalize:
        # since the imputation is used for dimension reduction, we want the dimension reduction to focus on structural variations, not illumination variances. In this case, a pose normalization kind of pre-classification is important to remove golden particles, so that they are not mixed with complexes of similiar structures
        v = (v - v.mean()) / v.std()

        if t is not None:      t = (t- t.mean()) / t.std()

    if t is None:   return v        # if the template is not given, no imputation

    #print 'neighbor_covariance__collect_info():  template given'     # just for debugging

    t_f = NF.fftshift(NF.fftn(t))
    v_f = NF.fftshift(NF.fftn(v))

    v_f[vm == 0] = t_f[vm == 0]     # this is the main imputation step

    v_f_if = N.real( NF.ifftn(NF.ifftshift(v_f)) )
    if normalize:   v_f_if = (v_f_if - v_f_if.mean()) / v_f_if.std()

    if N.all(N.isfinite(v_f_if)):
        return v_f_if
    else:
        print 'warning: imputation failed'
        return v
    

def impute_vols(v, vm, ang, loc, t=None, align_to_template=True, normalize=None):


    ang = N.array(ang, dtype=N.float)      # convert, just in case that the format is incorrect
    loc = N.array(loc, dtype=N.float)      # convert, just in case that the format is incorrect


    v_r = None
    vm_r = None
    t_r = None
    if align_to_template:
        v_r = GR.rotate_pad_mean(v, angle=ang, loc_r=loc)       ;           assert N.all(N.isfinite(v_r))
        vm_r = GR.rotate_mask(vm, angle=ang)                    ;           assert N.all(N.isfinite(vm_r))

        vi = impute_aligned_vols(t=t, v=v_r, vm= (vm_r > 0.5), normalize=normalize)

    else:

        ang_inv, loc_inv = AAL.reverse_transform_ang_loc(ang, loc)
        
        if t is not None:
            assert      N.all(N.isfinite(t))
            t_r = GR.rotate_vol_pad_mean(t, angle=ang_inv, loc_r=loc_inv)       ;       N.all(N.isfinite(t_r))

        vi = impute_aligned_vols(t=t_r, v=v, vm= (vm > 0.5), normalize=normalize)

    return {'vi':vi, 'v_r':v_r, 'vm_r':vm_r, 't_r':t_r}



def impute_vol_keys(vk, ang, loc, normalize, tk=None, align_to_template=True, cache=None):

    v  = cache.get_mrc(vk['subtomogram'])
    if not N.all(N.isfinite(v)):        raise Exception('error loading', vk['subtomogram'])

    vm = cache.get_mrc(vk['mask'])
    if not N.all(N.isfinite(vm)):        raise Exception('error loading', vk['mask'])


    t = None
    if tk is not None:        
        t = IV.get_mrc(tk['subtomogram'])       # we do not cache templates, in order to avoid catching prolems
        if not N.all(N.isfinite(t)):        raise Exception('error loading', tk['subtomogram'])
        #print 'imputation using', tk['subtomogram']

    return impute_vols(v=v, vm=vm, ang=ang, loc=loc, t=t, align_to_template=align_to_template, normalize=normalize)



# calculate the covariance between neighbor voxels, then take average
def neighbor_covariance_avg__parallel(self, data_json, segmentation_tg_op, normalize, n_chunk):

    start_time = time.time()

    data_json_copy = [_ for _ in data_json]
    inds = range(len(data_json_copy))

    tasks = []
    while data_json_copy:
        data_json_copy_part = data_json_copy[:n_chunk]
        inds_t = inds[:n_chunk]

        tasks.append(self.runner.task(module='tomominer.pursuit.multi.util', method='neighbor_covariance__collect_info', kwargs={'data_json':data_json_copy_part, 'segmentation_tg_op':segmentation_tg_op, 'normalize':normalize}))

        data_json_copy = data_json_copy[n_chunk:]
        inds = inds[n_chunk:]


    sum_global = None
    neighbor_prod_sum = None
    for res in self.runner.run__except(tasks):

        with open(res.result) as f:     re = pickle.load(f)
        os.remove(res.result)

        if sum_global is None:
            sum_global = re['sum']
        else:
            sum_global += re['sum']

        assert  N.all(N.isfinite(sum_global))

        
        if neighbor_prod_sum is None:
            neighbor_prod_sum = re['neighbor_prod_sum']    
        else:
            neighbor_prod_sum += re['neighbor_prod_sum']

        assert  N.all(N.isfinite(neighbor_prod_sum))
    
    avg_global = sum_global / len(data_json)
    neighbor_prod_avg = neighbor_prod_sum / len(data_json)


    shift = re['shift']
    cov = N.zeros(neighbor_prod_avg.shape)     # use float16 to reduce data storage
    for i in range(shift.shape[0]):
        cov[:,:,:,i] = neighbor_prod_avg[:,:,:,i] - avg_global * UV.roll(avg_global, shift[i,0], shift[i,1], shift[i,2])

    cov_avg = N.mean(cov, axis=3)

    print 'Calculated neighbor covariance for', len(data_json), 'subtomograms', ' : %2.6f sec' % (time.time() - start_time)

    return cov_avg


# collecting information for calculating the neighbor covariance, calculated at worker side
def neighbor_covariance__collect_info(self, data_json, segmentation_tg_op, normalize):

    sum_local = None
    neighbor_prod_sum = None

    for rec in data_json:
        if 'template' not in rec:  rec['template'] = None

        vri = impute_vol_keys(vk=rec, ang=rec['angle'], loc=rec['loc'], tk=rec['template'], align_to_template=True, normalize=normalize, cache=self.cache) ['vi']

        if (segmentation_tg_op is not None) and (rec['template'] is not None) and ('segmentation' in rec['template']):
            #print 'apply segmentation mask to the interpolated subtomgoram', rec['template']['segmentation']
            phi = IV.read_mrc(rec['template']['segmentation'])['value']

            vri_s = template_guided_segmentation(v=vri, m=(phi>0.5), op=segmentation_tg_op)          # why phi>0.5 instead of 0? Because we do not want to include boundary region?
            if vri_s is not None:
                vri = vri_s;                del vri_s
        
                assert  normalize is not None    
                if normalize:
                    vri_t = N.zeros(vri.shape)
                    vri_f = N.isfinite(vri)
                    if vri_f.sum() > 0:     vri_t[vri_f] = (vri[vri_f] - vri[vri_f].mean()) / vri[vri_f].std()
                    vri = vri_t;                    del vri_f, vri_t
                else:
                    vri_t = N.zeros(vri.shape)
                    vri_f = N.isfinite(vri)
                    vri_t[vri_f] = vri[vri_f]
                    if vri_f.sum() > 0:     vri_t[N.logical_not(vri_f)] = vri[vri_f].mean()
                    vri = vri_t;                    del vri_f, vri_t

        if sum_local is None:
            sum_local = vri
        else:
            sum_local += vri

        nei_prod = DU.neighbor_product(vri)
        if neighbor_prod_sum is None:
            neighbor_prod_sum = nei_prod['p']
        else:
            neighbor_prod_sum += nei_prod['p']


    re = {'sum':sum_local, 'neighbor_prod_sum':neighbor_prod_sum, 'shift':nei_prod['shift']}

    re_key = self.cache.save_tmp_data(re, fn_id=self.task.task_id)
    assert re_key is not None

    return re_key




# modified from tomominer.dimension_reduction.util.pca__stack_dif__parallel()
def data_matrix_collect__parallel(self, data_json, segmentation_tg_op, normalize, n_chunk, voxel_mask_inds=None):
    start_time = time.time()

    data_json_copy = [_ for _ in data_json]
    inds = range(len(data_json_copy))

    tasks = []
    while data_json_copy:
        data_json_copy_t = data_json_copy[:n_chunk]
        inds_t = inds[:n_chunk]
        tasks.append(self.runner.task(module='tomominer.pursuit.multi.util', method='data_matrix_collect__local', kwargs={'data_json':data_json_copy_t, 'segmentation_tg_op':segmentation_tg_op, 'normalize':normalize, 'inds':inds_t, 'voxel_mask_inds':voxel_mask_inds}))
        data_json_copy = data_json_copy[n_chunk:]
        inds = inds[n_chunk:]

    # merge all results
    red = None
    for res in self.runner.run__except(tasks):

        with open(res.result) as f:     re = pickle.load(f)
        os.remove(res.result)

        if red is None:     red = N.zeros( [len(data_json), re['mat'].shape[1]] )

        red[re['inds'],:] = re['mat']



    print 'Calculated matrix of', len(data_json), 'subtomograms', '%2.6f sec'%(time.time() - start_time)

    return red


# modified from tomominer.dimension_reduction.util.pca__stack_dif()
def data_matrix_collect__local(self, data_json, inds, segmentation_tg_op, normalize, voxel_mask_inds=None):
    
    mat = None
    for i, rec in enumerate(data_json):
        if 'template' not in rec:            rec['template'] = None

        vi = impute_vol_keys(vk=rec, ang=rec['angle'], loc=rec['loc'], tk=rec['template'], align_to_template=True, normalize=normalize, cache=self.cache)['vi']

        if (segmentation_tg_op is not None) and (rec['template'] is not None) and ('segmentation' in rec['template']):
            #print 'apply template segmentation guide to the interpolated subtomgoram', rec['template']['segmentation']
            phi = IV.read_mrc(rec['template']['segmentation'])['value']

            vi_s = template_guided_segmentation(v=vi, m=(phi>0.5), op=segmentation_tg_op)      # why phi>0.5 instead of 0? Because we do not want to include boundary region?
            if vi_s is not None:
                vi = vi_s;                  del vi_s
                
                assert  normalize is not None
                if normalize:
                    vi_t = N.zeros(vi.shape)
                    vi_f = N.isfinite(vi)
                    vi_t[vi_f] = (vi[vi_f] - vi[vi_f].mean()) / vi[vi_f].std()
                    vi = vi_t;                    del vi_f, vi_t
                else:
                    vi_t = N.zeros(vi.shape)
                    vi_f = N.isfinite(vi)
                    vi_t[vi_f] = vi[vi_f]
                    vi_t[N.logical_not(vi_f)] = vi[vi_f].mean()
                    vi = vi_t;                    del vi_f, vi_t



        vi = vi.flatten()

        if voxel_mask_inds is not None:
             vi =  vi[voxel_mask_inds]

        if mat is None:
            mat = N.zeros([len(data_json), vi.size])
        
        mat[i, :] = vi    

    re = {'mat':mat, 'inds':inds}

    re_key = self.cache.save_tmp_data(re, fn_id=self.task.task_id)
    assert re_key is not None
    return re_key



# calculate average covariance between neighbor voxels, then gaussian smooth and segment to identify a small amount of voxels as features for PCA analysis
# paramters:    data_json_model: subtomograms for training PCA,     data_json_embed: subtomograms to project into low dimension with PCA
def covariance_filtered_pca(self, data_json_model=None, data_json_embed=None, normalize=None, segmentation_tg_op=None, n_chunk=100, max_feature_num=None, pca_op=None):
    print 'Dimension reduction'
   
    start_time = time.time()
    
    # If training set is None then use whole data_json_embed file to compute dimensions
    # But then make data_json_embed None to avoid repeated computation
    if data_json_model is None:
        assert      data_json_embed is not None
        data_json_model = data_json_embed
        data_json_embed = None      # just to avoid repeated computation


    cov_avg = None
    cov_avg__feature_num_cutoff = None
    if (max_feature_num==None) or (max_feature_num < 0):
        voxel_mask_inds = None
    else:
        # calculate average covariance, or load existing calculated result
        cov_avg = neighbor_covariance_avg__parallel(self=self, data_json=data_json_model, segmentation_tg_op=segmentation_tg_op, normalize=normalize, n_chunk=n_chunk)
        cov_avg_max = cov_avg.max()
        if (not N.isfinite(cov_avg_max)) or (cov_avg_max <= 0):
            raise Exception('cov_avg.max(): ' + repr(cov_avg_max))


        # restrict the max number of featrues to be processed to be less than max_feature_num, in order to avoid the dimension reduction to be over time consuming
        cov_avg_i = N.argsort((-cov_avg), axis=None)
        cov_avg__feature_num_cutoff = cov_avg.flatten()[cov_avg_i[min(max_feature_num, cov_avg_i.size-1)]]
        cov_avg__feature_num_cutoff = max(cov_avg__feature_num_cutoff, 0)       # such cutoff has to be at least 0
        voxel_mask_inds = N.flatnonzero(cov_avg > cov_avg__feature_num_cutoff)

        #print 'max covariance %f, cov_avg__feature_num_cutoff %f, %d features selected for dimension reduction'%(cov_avg.max(), cov_avg__feature_num_cutoff, len(voxel_mask_inds))


    # perform pca using only voxels / features indicated by vol_msk
    mat = data_matrix_collect__parallel(self=self, data_json=data_json_model, segmentation_tg_op=segmentation_tg_op, normalize=normalize, n_chunk=n_chunk, voxel_mask_inds=voxel_mask_inds)
   
    # centerize collected matrix, for PCA
    mat_mean = mat.mean(axis=0)
    for i in range(mat.shape[0]):   mat[i,:] -= mat_mean


    # weight missing value regions as 0
    empca_weight = N.zeros(mat.shape, dtype=float)
    empca_weight[N.isfinite(mat)] = 1.0

    if cov_avg is not None:
        # weight every feature according to its corresponding average correlation
        cov_avg_v = cov_avg.flatten()
        for i, ind_t in enumerate(voxel_mask_inds):
            empca_weight[:,i] *= cov_avg_v[ind_t]

    while True:
        try:
            import tomominer.dimension_reduction.empca as drempca
            pca = drempca.empca( data=mat, weights=empca_weight, nvec=pca_op['n_dims'], niter=pca_op['n_iter'] )     # note: need to watch out the R2 values to see how much variation can be explained by the estimated model, if the value is small, need to increase dims
        except N.linalg.linalg.LinAlgError as e:
            # in the native built version, sometimes such exception raises, strange, in such case just try again
            print 'Warning:', e
            print 'retry'
            pass
        else:
            break


    if (data_json_embed is None) or ():
        red = pca.coeff
    else:
        # if the embedding subtomograms is different to the subtomograms used for constructing PCA model, we collect these subtomograms and project them to low dimension using the trained PCA model
        ''' Because there might be junks inside. This is also useful to reduce computation time. Certainly the PCA may omit the important dimension of some classes and such classes may not form informative clusters, however we can always remove resulting clusters and use remaining subtomograms to repeat the process.'''
        mat_embed = data_matrix_collect__parallel(self=self, data_json=data_json_embed, segmentation_tg_op=segmentation_tg_op, normalize=normalize, n_chunk=n_chunk, voxel_mask_inds=voxel_mask_inds)
        red = N.dot(mat_embed, pca.eigvec.T)

    print "PCA with covariange thresholding  : %2.6f sec" % (time.time() - start_time)

    return {'red':red, 'cov_avg':cov_avg, 'cov_avg__feature_num_cutoff':cov_avg__feature_num_cutoff, 'voxel_mask_inds':voxel_mask_inds}



# missing wedge corrected pca, just for comparison purpose
def covariance_filtered_pca_with_wedge(self, data_json, n_chunk, op, pass_dir, max_feature_num=None):

    raise   Exception('todo: seperate data_json into data_json_model and data_json_embed, like done in covariance_filtered_pca()')

    data = CC.parse_data(data_json)

    gavg_re = vol_avg_parallel_global(self=self, data_json=data_json, n_chunk=n_chunk, out_dir=pass_dir)
    avg_key = gavg_re['template_keys']['subtomogram']
    avg_mask_key = gavg_re['template_keys']['mask']

    cfp_re = DU.covariance_filtered_pca(self=self, avg_key=(avg_key, avg_mask_key), data=data, n_chunk=n_chunk, dims=op['dim_reduction']['pca']['n_dims'], pca_op=op['dim_reduction']['pca'], out_dir=pass_dir, max_feature_num=max_feature_num)

    return cfp_re


# just convert labels and data_json into a cluster dictionary
# divide data_json amonf cluster labels
def labels_to_clusters(data_json, labels, cluster_mode=None, ignore_negative_labels=True):
    clusters = {}
    for l, d in zip(labels, data_json):
        if ignore_negative_labels and (l < 0):      continue

        if l not in clusters:       clusters[l] = {}
        if 'cluster_mode' not in clusters[l]:     clusters[l]['cluster_mode'] = cluster_mode
        if 'data_json' not in clusters[l]:     clusters[l]['data_json'] = []
        clusters[l]['data_json'].append(d)

    return clusters



def kmeans_clustering(x, k):


    warnings.filterwarnings('once')
    from sklearn.cluster import KMeans
    warnings.filterwarnings('error')
    
    if False:
        import multiprocessing
        km = KMeans(n_clusters=k, n_jobs=multiprocessing.cpu_count())
    else:
        km = KMeans(n_clusters=k, n_jobs = 10)       # by default, such kmeans method uses kmeans++ initialization, and try 10 times to reduce the chance of falling to local minima. In such case, hopfully we have a high chance to get clusters inside dense regions of data points, becase this makes the objective function have a higher value

    labels = km.fit_predict(x)

    # re-index labels, just to prevent the label number to be too large
    # labels them starting from 0 to (num_of_unique_labels -1)
    labels_t = N.zeros(labels.shape, dtype=N.int) + N.nan
    label_count = 0
    for l in N.unique(labels):
        labels_t[labels == l] = label_count
        label_count += 1

    labels = labels_t.astype(N.int)


    return labels



def model_based_clustering__spherical(x, k):
    import tomominer.cluster.model_based as CM

    best = {'bic':float('-inf')}

    for G in range(2, k+1):
        lbl_init = N.floor( N.random.rand(x.shape[0]) * G ) 
        re = CM.spherical_model({'x':x, 'lbl_init':lbl_init, 'G':G, 'max_iteration_num':1000, 'stop_tolorence':1e-5})

        print G, re['bic']

        if re['bic'] > best['bic']:
            best = re

    return best



# calculate spectial snr and infer fourier shell correlation for each cluster
def cluster_ssnr_fsc(self, clusters, n_chunk, op=None):
    start_time = time.time()

    sps = SS.ssnr_parallel(self=self, clusters=clusters, n_chunk=n_chunk, op=op)

    #print 'cluster_ssnr_fsc took time %.1f'%(time.time() - start_time)

    return sps


def vol_avg__local(self, data_json, op=None, return_key=True):

    vol_sum = None
    mask_sum = None

    # temporary collection of local volume, and mask.
    for rec in data_json:
        if self.work_queue.done_tasks_contains(self.task.task_id):            raise Exception('Duplicated task')         # terminate alignment process when the corresponding task is already completed by another worker

        if op['with_missing_wedge'] or op['use_fft']:
            # in this case, do not interpolate, only use impute_vol_keys() rotate subtomogram and corresponding masks
            # impute_vol_keys below impute the missing region with template, as tk=None in call below, no template is used and just rotated volume is retured
            in_re = impute_vol_keys(vk=rec, ang=rec['angle'], loc=rec['loc'], tk=None, align_to_template=True, normalize=False, cache=self.cache)
            # here v_r is rotated subtomogram according to angle info stored in its data_json
            vt = in_re['v_r']
        # else use templates (cluster average from last iteration) stored as part of data_json, to impute the subtomograms and interpolate and use those subtomograms for averaging
        else:
            raise Exception('following options need to be re-considered')
            # in this case, rotate, interpolate and normalize subtomogram
            if 'template' not in rec:       rec['template'] = None
            in_re = impute_vol_keys(vk=rec, ang=rec['angle'], loc=rec['loc'], tk=rec['template'], align_to_template=True, normalize=True, cache=self.cache)
            vt = in_re['vi']

        if vol_sum is None:     vol_sum = N.zeros(vt.shape, dtype=N.float64, order='F')
        vol_sum += vt

        if mask_sum is None:        mask_sum = N.zeros(in_re['vm_r'].shape, dtype=N.float64, order='F')
        mask_sum += in_re['vm_r']

        

    re = {'vol_sum':vol_sum, 'mask_sum':mask_sum, 'vol_count':len(data_json), 'op':op}

    if return_key:
        re_key = self.cache.save_tmp_data(re, fn_id=self.task.task_id)
        assert re_key is not None
        return {'key':re_key}

    else:

        return re


def cluster_averaging_vols(self, clusters, op={}):

    start_time = time.time()

    if op['centerize_loc']:
        clusters_cnt = copy.deepcopy(clusters)
        # here clusters is list of (list of subtomogram in cluster) OR complete data_json in each cluster is passes as element of list
        for c in clusters_cnt:
            # create 2d array loc [num of subotmograms in cluster c][3]
            loc = N.zeros( (len(clusters_cnt[c]), 3) )
            # load loc parameter from data_json to loc array just created
            for i, rec in enumerate(clusters_cnt[c]):      loc[i,:] = N.array(rec['loc'], dtype=N.float)
            # loc - mean_of_loc OR centerize loc
            loc -= N.tile( loc.mean(axis=0), (loc.shape[0], 1) )       # substract mean so that mean of loc is equal to zero vector. The purpose of centerize_loc is to reduce the chance of clipping brought by large displacements
            assert N.all(N.abs(loc.mean(axis=0)) <= 1e-10)
            
            # update new centerized loc to data_json in clusters
            for i, rec in enumerate(clusters_cnt[c]):      rec['loc'] = loc[i]

        # update clusters
        clusters = clusters_cnt


    # parallel averaging
    tasks = []
    for c in clusters:
        while clusters[c]:
            part = clusters[c][:op['n_chunk']]

            op_t = copy.deepcopy(op)
            op_t['cluster'] = c
            tasks.append(self.runner.task(module='tomominer.pursuit.multi.util', method='vol_avg__local', kwargs={'data_json':part, 'op':op_t, 'return_key':True}))

            clusters[c] = clusters[c][op['n_chunk']:]
    
    # calculate averages for each cluster
    # cluster_sizes = dict((c,len(clusters[c])) for c in clusters)
    cluster_sums = {}
    cluster_mask_sums = {}
    cluster_sizes = {}
    for res in self.runner.run__except(tasks):

        with open(res.result['key']) as f:     re = pickle.load(f)
        os.remove(res.result['key'])
        
        oc = re['op']['cluster']
        ms = re['mask_sum']
        s = re['vol_sum']
        vc = re['vol_count']

        if oc not in cluster_sums:      
            cluster_sums[oc] = s
        else:
            cluster_sums[oc] += s

        if oc not in cluster_mask_sums:
            cluster_mask_sums[oc] = ms
        else:
            cluster_mask_sums[oc] += ms

        if oc not in cluster_sizes:
            cluster_sizes[oc] = vc
        else:
            cluster_sizes[oc] += vc

        del oc, ms, s, vc       # just to prevent this temp variables to be reused outside the block



    cluster_avg_dict = {}

    for c in cluster_sums:
        assert cluster_sizes[c] > 0
        assert cluster_mask_sums[c].max() > 0

        if op['use_fft']:

            ind = cluster_mask_sums[c] >= op['mask_count_threshold']
            if ind.sum() == 0:      continue            # we ignore small clusters, because their resulting average is purly zeros


            cluster_sums_fft = NF.fftshift( NF.fftn(cluster_sums[c]) )
            
            cluster_avg = N.zeros(cluster_sums_fft.shape, dtype = N.complex)

            cluster_avg[ind] = cluster_sums_fft[ind] / cluster_mask_sums[c][ind]
            
            cluster_avg = N.real( NF.ifftn(NF.ifftshift(cluster_avg)) )

        else:
            cluster_avg = cluster_sums[c] / cluster_sizes[c]


        if op['mask_binarize']:
            # just make a binary mask
            cluster_mask_avg = cluster_mask_sums[c] >= op['mask_count_threshold']
        else:
            cluster_mask_avg = cluster_mask_sums[c] / cluster_sizes[c]


        cluster_avg_dict[c] = {'vol':cluster_avg, 'mask':cluster_mask_avg}

    
    # smooth averages, if needed
    if 'smooth' in op:
        print 'smoothing', op['smooth']['mode'], 
        for c in cluster_avg_dict:
            if c not in op['smooth']['fsc']:    continue

            cluster_avg_dict[c]['smooth'] = {'vol-original':cluster_avg_dict[c]['vol']}         # just make a backup, in case there is a need to use such data
            s = cluster_averaging_vols__smooth(v=cluster_avg_dict[c]['vol'], fsc=op['smooth']['fsc'][c], mode=op['smooth']['mode'])
            cluster_avg_dict[c]['vol'] = s['v']

            if 'fit' in s:
                cluster_avg_dict[c]['smooth']['fit'] = s['fit']
                print 'average', c, 'sigma ', cluster_avg_dict[c]['smooth']['fit']['c'], '    ',

        print

    #print "Averaging subtomogram sets %2.6f sec" % (time.time() - start_time)

    return {'cluster_avg_dict':cluster_avg_dict, 'cluster_sums':cluster_sums, 'cluster_mask_sums':cluster_mask_sums, 'cluster_sizes':cluster_sizes}


def cluster_averaging_vols__smooth(v, fsc, mode):

    re = {}

    assert  N.all(fsc >= 0)
    if fsc.max() == 0.0:   return {'v':v}     

    if mode == 'fsc_direct':
        band_pass_curve = fsc
    elif mode == 'fsc_gaussian':
        import tomominer.fitting.gaussian.one_dim as FGO
        bands = N.array(range(len(fsc)))
        fit = FGO.fit__zero_mean__multi_start(x=bands, y=fsc)

        if fit['c'] < bands.max():
            # in such case, the gaussian function fitting looks reasonable
            re['fit'] = fit
            band_pass_curve = FGO.fit__zero_mean__gaussian_function(x=bands, a=fit['a'], c=fit['c'])
        else:
            # if c is too wide, we do not perform band pass filtering
            re['v'] = v
            return re
        
    else:
        raise Exception('mode')

    import tomominer.filter.band_pass as IB
    re['v'] = IB.filter_given_curve(v=v, curve=band_pass_curve)
    return re



def cluster_averaging(self, clusters, op={}):

    #print 'cluster_averaging()', op
   
    cav =  cluster_averaging_vols(self, clusters=clusters, op=op)

    if not os.path.isdir(op['out_dir']):      os.makedirs(op['out_dir'])
    # This cluster.pickle is saved under directory clus_avg
    with open(os.path.join(op['out_dir'], 'cluster.pickle'), 'wb') as f:      pickle.dump(cav, f, protocol=-1)
############################################################################################
    # save averages
    cluster_avg_dict = cav['cluster_avg_dict']
    template_keys = {}
    for c in cluster_avg_dict:
        template_keys[c] = {}
        template_keys[c]['cluster'] = c

        vol_avg_out_key = os.path.join(op['out_dir'], 'clus_vol_avg_%03d.mrc'%(c))
        IV.put_mrc(N.array(cluster_avg_dict[c]['vol'], order='F'), vol_avg_out_key)
        template_keys[c]['subtomogram'] = vol_avg_out_key               # use key 'subtomogram' for the consistency with json record in data_config, for calling align_keys()

        if 'smooth' in cluster_avg_dict[c]:
            vol_avg_original_out_key = os.path.join(op['out_dir'], 'clus_vol_avg_original_%03d.mrc'%(c))
            IV.put_mrc(N.array(cluster_avg_dict[c]['smooth']['vol-original'], order='F'), vol_avg_original_out_key)
            template_keys[c]['subtomogram-original'] = vol_avg_original_out_key
        else:
            vol_avg__riginal_out_key = None

        mask_avg_out_key = os.path.join(op['out_dir'], 'clus_mask_avg_%03d.mrc'%(c))
        IV.put_mrc(N.array(cluster_avg_dict[c]['mask'], order='F'), mask_avg_out_key)
        template_keys[c]['mask'] = mask_avg_out_key

        if 'pass_i' in op:      template_keys[c]['pass_i'] = op['pass_i']


    return {'template_keys':template_keys}


def cluster_average_select(self, data_json, labels, template_keys, op):

    template_keys = copy.deepcopy(template_keys)      # since this dict is going to be modified, make a deep copy to keep original copy

    # construct inverse mapping
    template_keys_inv = {}
    for c in template_keys:     template_keys_inv[template_keys[c]['subtomogram']] = c



    clusters = defaultdict(list)
    for c, rec in zip(labels, data_json):
        if c < 0:   continue
        clusters[c].append(rec)


    # use vol_cluster_label for recording of cluster membership 
    vol_cluster_label = {}
    for c in clusters:
        for rec in clusters[c]:
            vol_cluster_label[rec['subtomogram']] = c



    # remove small clusters
    for c in clusters:
        if len(clusters[c]) < op['min_cluster_size']:
            template_keys.pop(c, None)      # delete from tomominer.template_keys, if any


    cluster_sizes = [(-len(clusters[c]), c) for c in template_keys];    cluster_sizes.sort();       print 'Cluster sizes ' + repr(cluster_sizes)
    ordered_clusters = [c for s,c in cluster_sizes]         ; print 'Ordered clusters ' + repr(ordered_clusters)


    if len(template_keys) <= 2:
        return {'selected_templates':template_keys}


    # calculate pairwise alignment
    tasks = []
    for c0 in template_keys:
        for c1 in template_keys:
            if c1 <= c0: continue
            tasks.append(self.runner.task(module='tomominer.pursuit.multi.util', method='align_keys', kwargs={'v1k':template_keys[c0], 'v2k':template_keys[c1], 'op':op['align_op']}))

    ress = [_ for _ in self.runner.run__except(tasks)]
    
    align_score = defaultdict(dict)
    for res in ress:
        c0 = template_keys_inv[res.result['v1_key']['subtomogram']]
        c1 = template_keys_inv[res.result['v2_key']['subtomogram']]
        align_score[c0][c1] = res.result['score']
        align_score[c1][c0] = res.result['score']

        if res.result['err'] is not None:     print 'cluster_average_select() alignment error ' + repr( res.result['err'] )


    print 'Template alignment correlation scores ' + repr(align_score)


    # perform hierarchical clustering, and get optimal set of clusters, then for each cluster, select one average from largest subtomogram cluster
    align_score_dist_mat = N.zeros( (len(ordered_clusters), len(ordered_clusters)) )
    for i0, c0 in enumerate(ordered_clusters):
        for i1, c1 in enumerate(ordered_clusters):
            if c0 == c1: continue

            score_t = align_score[c0][c1]
            if not N.isfinite(score_t):     score_t = -1.0
            
            assert      N.abs(score_t) <= 1.0
            align_score_dist_mat[i0, i1] = N.sqrt(2 - 2 * score_t)

    import tomominer.cluster.hierarchy as H
    align_score_dist_mat__hc = H.cluster_dist_mat(align_score_dist_mat)
    osb = H.optimal_silhouette(align_score_dist_mat, align_score_dist_mat__hc, dist_threhold_max=N.sqrt(2 - 2 * op['align_corr_threshold']))


    for i0, c0 in enumerate(ordered_clusters):
        for i1, c1 in enumerate(ordered_clusters):
            if i1 <= i0: continue

            if osb['labels'][i0] != osb['labels'][i1]:      continue

            template_keys.pop(c1, None)     # remove c1, since it is from a smaller subtomogram cluster compared to c0, and it has similar shape like c0 according to hierarchical clustering



    '''
    for i0, c0 in enumerate(ordered_clusters):
        if c0 not in template_keys:     continue
        
        # remove all templates that have identical shape with the template of c0
        for i1, c1 in enumerate(ordered_clusters):
            if i1 <= i0: continue
            if c1 not in template_keys: continue
            if align_score[c0][c1] >= op['align_corr_threshold']:       template_keys.pop(c1, None)     # delete from tomominer.template_keys, if any        
    '''

    assert      len(template_keys.keys()) > 0
    print 'Selected clusters for subtomogram alignment ' + repr(template_keys.keys())

    return {'selected_templates':template_keys, 'align_score':align_score, 'ordered_clusters':ordered_clusters, 'align_score_dist_mat':align_score_dist_mat, 'align_score_dist_mat__hc':align_score_dist_mat__hc, 'optimal_silhouette':osb}


# select the largest set of non-overlapping clusters with highest FSC scores
def cluster_average_select_fsc(self, cluster_info, cluster_info_stat, op=None, debug=False):

    ci = []
    for i in cluster_info:
        for c in cluster_info[i]:
            if 'fsc' not in cluster_info[i][c]: continue
            if len(cluster_info[i][c]['data_json']) < op['cluster']['size_min']:      continue
            if (not op['keep_non_specific_clusters']) and ('is_specific' in cluster_info_stat[i][c]) and (cluster_info_stat[i][c]['is_specific'] is not None):       continue      # ignore non-specific clusters, if keep_non_specific_clusters==False

            ci.append( cluster_info[i][c] )        

    ci = sorted(ci, key = lambda x:float(- x['fsc'].sum()))         # use float to make sure that the value used is a single scalar, not a numpy array, in order to avoid confusion in sorting

    # select those non-overlap clusters with highest FSC
    ci_cover = []
    covered_set = set()
    for ci_t in ci:
        if 'template_key' not in ci_t:      continue

        subtomograms_t = set(_['subtomogram'] for _ in ci_t['data_json'])

        overlap_ratio_t = (float(len( covered_set.intersection(subtomograms_t) )) / len(subtomograms_t))

        if overlap_ratio_t <= op['cluster']['overlap_ratio']:
            # note: op['cluster']['overlap_ratio'] has to be small, say less than 0.1, so that the clusters will not be made non-specific due to member overlaps!!!!!
            ci_cover.append(ci_t)
            covered_set.update(subtomograms_t)

        if debug:            print   ci_t['pass_i'], ci_t['cluster'], len(ci_t['data_json']), ci_t['fsc'].sum(), overlap_ratio_t


    del ci

    #print 'cluster_average_select_fsc() initially selected %d averages covering %d subtomograms '%(len(ci_cover), len(covered_set))
    print 'Set sizes', sorted([len(_['data_json']) for _ in ci_cover ])

    tk = {}
    for i, cc in enumerate(ci_cover):
        tk[i] = copy.deepcopy(cc['template_key'])
        tk[i]['id'] = i
        assert  tk[i]['pass_i'] == cc['pass_i']
        assert  tk[i]['cluster'] == cc['cluster']

    tk_selected = set(tk[_]['subtomogram'] for _ in tk)

    tk_info = {}
    for i in cluster_info:
        for c in cluster_info[i]:
            if 'template_key' not in cluster_info[i][c]:        continue

            tk_subtomogram = cluster_info[i][c]['template_key']['subtomogram']
            if tk_subtomogram not in tk_selected:   continue                # to reduce storage, we only keep information of templates selected in current pass
            tk_info[tk_subtomogram] = cluster_info[i][c]



    if op['keep_non_specific_clusters']:
        # RESET specificity of selected clusters!!!
        for i, tkt in tk.iteritems():   cluster_info_stat[tkt['pass_i']][tkt['cluster']]['is_specific'] = None

    assert  len(tk) > 0

    return {'selected_templates':tk, 'tk_info':tk_info}




# perform pairwise alignment and use hierarchical clustering to seperate averages into different average clusters, then for each average cluster select the best average in terms of FSC
def cluster_average_select_fsc_alignment_clustering(self, cluster_info, op=None):

    casf = cluster_average_select_fsc(self=self, cluster_info=cluster_info, op=op)
    tk = casf['selected_templates']


    # calculate pairwise alignment
    tasks = []
    for i0,tk0 in enumerate(tk):
        for i1,tk1 in enumerate(tk):
            if i1 <= i0: continue
            tasks.append(self.runner.task(module='tomominer.pursuit.multi.util', method='align_keys', kwargs={'v1k':tk0, 'v2k':tk1, 'op':op['align_op']}))


    pair_align = defaultdict(dict)
    for res_t in self.runner.run__except(tasks):
        res = res_t.result

        vk0 = res['v1_key']['subtomogram']
        vk1 = res['v2_key']['subtomogram']
        pair_align[vk0][vk1] = res

        if N.abs(res['score']) > 1.0:        
            print 'cluster_average_select_fsc()  warning: score %f'%(res['score'])
            res['score'] = N.sign(res['score'])

        if not N.isfinite(res['score']):     res['score'] = -1.0

        # store reversed alignment of c1 and c0
        pair_align[vk1][vk0] = copy.deepcopy(res)
        ang_rev, loc_rev = AAL.reverse_transform_ang_loc(ang=res['angle'], loc_r=res['loc'])
        pair_align[vk1][vk0].update( {'angle':ang_rev, 'loc':loc_rev} )

        if res['err'] is not None:     print 'cluster_average_select_fsc() alignment error ' + repr( res['err'] )


    aligned_keys = defaultdict(dict)

    # perform hierarchical clustering, and get optimal set of clusters, then for each cluster, select one best average
    align_score_mat = N.zeros( (len(tk), len(tk)) )
    for i0, tk0 in enumerate(tk):
        for i1, tk1 in enumerate(tk):
            if i0 == i1:
                align_score_mat[i0, i1] = 1.0
                continue

            vk0 = tk0['subtomogram']
            vk1 = tk1['subtomogram']

            score_t = pair_align[vk0][vk1]['score']
            
            align_score_mat[i0, i1] = score_t

    align_score_dist_mat = N.sqrt(2 - 2 * align_score_mat)

    import tomominer.cluster.hierarchy as H
    align_score_dist_mat__hc = H.cluster_dist_mat(align_score_dist_mat)
    osb = H.optimal_silhouette(align_score_dist_mat, align_score_dist_mat__hc)
   
    if 'labels' in osb: 
        template_labels = osb['labels']
    else:
        print 'warning: no optimal labeles obtained'
        template_labels = [i for i,_ in enumerate(tk)]


    # ---------------------------------------------------------------------------------
    # for each set of clusters, select only the first one (i.e. the one with largest FSC)
    tk_sel_flag = [True] * len(tk)

    for i0, tk0 in enumerate(tk):
        for i1, tk1 in enumerate(tk):
            if i1 <= i0: continue
            if template_labels[i0] != template_labels[i1]:      continue

            tk_sel_flag[i1] = False


    stk = {}
    for i, f in enumerate(tk_sel_flag):
        if f == True:   
            assert template_labels[i] not in stk
            stk[template_labels[i]] = tk[i]


    return {'selected_templates':stk, 'tk_fsc':casf['tk_fsc'], 'tk_pass':casf['tk_pass']}




#  if a cluster C has members that best match the average Da of another cluster D, use the MannWightny test to see if the match scores of C's members against Ca is significantly higher than against Da. If not significantly higher, and if D has a higher fsc than C does, then mark cluster C as nonspecific and ignore further processing. Also, re-assign the best alignments on only remaining specific clusters
# parameters:       ci: cluster_info     cis: cluster_info_stat     al: template alignment results      tk: selected templates
# results: if the is_specific field is None, it means the selected cluster is specific, otherwise, it is non-specific, and it will have a record of the respecting other cluster that makes this cluster non-specific

def cluster_removal_according_to_center_matching_specificity(ci, cis, al, tk, significance_level, test_type=0, test_sample_num_min=10):

    tk = copy.deepcopy(tk)
    tkd = {tk[_]['subtomogram']:_ for _ in tk}      # mapping from tomominer.template key to cluster id


    # collect FSC for each matched templates
    tk_fsc = {}
    for pass_i in ci:
        for c in ci[pass_i]:
            ci0 = ci[pass_i][c]
            if 'template_key' not in ci0:   continue
            tk0 = ci0['template_key']['subtomogram']

            if tk0 not in tkd:     continue

            tk_fsc[tk0] = ci0['fsc'].sum()


    # remove non-specific (based on Mann Whitney test) and low FSC clusters
    non_specific_clusters = []
    wilcoxion_stat = defaultdict(dict)      # this keeps a record of wilcoxion statistics

    for pass_i in ci:
        for ci_c0 in ci[pass_i]:
            ci0 = ci[pass_i][ci_c0]
            cis0 = cis[pass_i][ci_c0]

            if 'is_specific' not in cis0:     cis0['is_specific'] = None
            if cis0['is_specific'] is not None:      continue       # ignore non-specific cluster
            
            if 'template_key' not in ci0:       continue
            tk0 = ci0['template_key']['subtomogram']
            if tk0 not in tkd:  continue        # this means the cluster is not selected in current pass, ignore this cluster

            c0 = tkd[tk0]

            # collect the alignment scores of all the subtomograms of ci0
            ci0s = set( str(_['subtomogram']) for _ in ci0['data_json']  )      # set of subtomograms of ci0
            al0 = [_ for _ in al if _['vol_key']['subtomogram'] in ci0s]        # collection of alignment scores of subtomograms of ci0

            ss = {}     # scores collected for each template
            for c1 in tk:
                tk1 = tk[c1]['subtomogram']

                if tk1 not in tkd:  continue        # this means the template is already deleted because its corresponding cluster is not specific

                if tk_fsc[tk1] < tk_fsc[tk0]:      continue     # we require cluster c0 to be WORSE than c1, but we also need to include ss[c0], so instead of using <=, we use <
                  
                ss[c1] = N.array( [_['align'][c1]['score'] for _ in al0]  )

            # compare alignment scores between matched templates
            for c1 in ss:
                if c1 == c0:    continue

                tk1 = tk[c1]['subtomogram']

                assert tk1 is not tk0
                assert tk1 in tkd
                assert  tk_fsc[tk1] > tk_fsc[tk0]        # we require cluster c0 to be WORSE than c1

                ind_t = N.logical_and(N.isfinite(ss[c0]), N.isfinite(ss[c1]))
                if ind_t.sum() < test_sample_num_min:    continue
                if N.all(ss[c0][ind_t] > ss[c1][ind_t]):      continue        # this means the template matching scores are all higher in ss[c0], best matchings are consistant


                is_specific = None


                if test_type == 0:
                    t_, p_ = SCS.wilcoxon(ss[c1][ind_t], ss[c0][ind_t])
                    if p_ > significance_level:
                        # this means that alignment scores in ss[c1] are similiar to ss[c0], which means that the best matchings to tk0 is not significantly better than the best matchings to tk1. IMPORTANT: the larger the significance level value, the less chance that c0 is determined as non-specific respect to c1, the more chance that there are redundant clusters
                        is_specific = {'tk0':tk[c0], 'tk1':tk[c1], 'stat':{'t':t_, 'p':p_}}           # Need to record original pass and cluster id, and selection cluster id and wilcoxion test results
                    elif N.median(ss[c1][ind_t]) > N.median(ss[c0][ind_t]):
                        # this means that alignment scores in ss[c1] are significantly higher than ss[c0], which means that the best matchings to tk0 is not significantly better than the best matchings to tk1
                        is_specific = {'tk0':tk[c0], 'tk1':tk[c1], 'stat':{'t':t_, 'p':p_}}


                elif test_type == 1:
                    t_, p_ = SCS.mannwhitneyu(ss[c1][ind_t], ss[c0][ind_t])     # test if ss[c1] are significantly less than ss[c0]
                    if p_ >= significance_level:
                        # this means that alignment scores in ss[c1] are not significantly smaller than ss[c0], which means that the best matchings to tk0 is not significantly better than the best matchings to tk1. IMPORTANT: the larger the significance level value, the less chance that c0 is determined as non-specific respect to c1, the more chance that there are redundant clusters
                        is_specific = {'tk0':tk[c0], 'tk1':tk[c1], 'stat':{'t':t_, 'p':p_}}
                else:
                    raise AttributeError('test_type')

                wilcoxion_stat[c0][c1] = {'t':t_, 'p':p_, 'median_c0':N.median(ss[c0][ind_t]), 'median_c1':N.median(ss[c1][ind_t])}

                if is_specific is not None:
                    cis0['is_specific'] = is_specific
                    del tkd[tk0]        # remove the average tk0 of cluster c0 from the matched key list
                    non_specific_clusters.append(ci0)
                    break

    # remove non-specific clusters' averages
    none_specific_cluster_ids = []
    for c in tk.keys():
        if tk[c]['subtomogram'] in tkd:     continue        # this means the cluster is specific

        del tk[c]
        none_specific_cluster_ids.append(c)

    # re-calculate best alignment
    for al_ in al:

        # select the best
        best = {}
        best['score'] = (-N.inf)
        best['template_id']   = None
        best['angle'] = N.random.random(3) * (N.pi * 2)      # mxu: randomly assign an angle
        best['loc'] = N.zeros(3)

        for c in tk:
            al_c = al_['align'][c]

            if al_c['score'] > best['score']:
                best['score']  = al_c['score']
                best['angle']  = al_c['angle']
                best['loc']    = al_c['loc']
                best['template_id'] = c     # since there are templates in original form and aligned into a common frame, we use template id instead of template key

        al_['best'] = best


    print len(non_specific_clusters), 'redundant averages detected', none_specific_cluster_ids      ;       sys.stdout.flush()

    return {'non_specific_clusters':non_specific_clusters, 'wilcoxion_stat':wilcoxion_stat}



def cluster_average_align_common_frame__pairwise_alignment(self, template_keys, align_op):

    # construct inverse mapping
    template_keys_inv = {}
    for c in template_keys:     template_keys_inv[template_keys[c]['subtomogram']] = c

    # calculate pairwise alignment
    tasks = []
    for c0 in template_keys:
        for c1 in template_keys:
            if c1 <= c0: continue
            tasks.append(self.runner.task(module='tomominer.pursuit.multi.util', method='align_keys', kwargs={'v1k':template_keys[c0], 'v2k':template_keys[c1], 'op':align_op}))

    
    pair_align = defaultdict(dict)
    for res_t in self.runner.run__except(tasks):
        res = res_t.result

        c0 = template_keys_inv[res['v1_key']['subtomogram']]
        c1 = template_keys_inv[res['v2_key']['subtomogram']]
        pair_align[c0][c1] = res
        pair_align[c0][c1].update({'c0':c0, 'c1':c1})

        # store reversed alignment of c1 and c0
        pair_align[c1][c0] = copy.deepcopy(res)
        ang_rev, loc_rev = AAL.reverse_transform_ang_loc(ang=res['angle'], loc_r=res['loc'])
        pair_align[c1][c0].update( {'angle':ang_rev, 'loc':loc_rev, 'c0':c1, 'c1':c0} )

        if res['err'] is not None:     print 'cluster_average_align_common_frame__pairwise() alignment error ' + repr( res['err'] )


    return pair_align




'''
simply align averages against the best average, this is useful when the averages are all similiar

we assume a cluster with smaller id has better quality
'''
def cluster_average_align_common_frame__single_best(self, tk, align_op, pass_dir):
    print 'cluster_average_align_common_frame__single_best()'

    out_dir = os.path.join(pass_dir, 'common_frame')
    if not os.path.isdir(out_dir):      os.mkdir(out_dir)

    # first perform pairwise alignment
    pa = cluster_average_align_common_frame__pairwise_alignment(self=self, template_keys=tk, align_op=align_op)

    c_best = int(N.min([_ for _ in tk]))
    
    # produce alignments
    tka = {}
    for c in tk:
        tka[c] = {}
        tka[c]['id'] = c
        tka[c]['pass_i'] = tk[c]['pass_i']
        tka[c]['cluster'] = tk[c]['cluster']
        tka[c]['subtomogram'] = os.path.join(out_dir, 'clus_vol_avg_%03d.mrc'%(c,))
        tka[c]['mask'] = os.path.join(out_dir, 'clus_mask_avg_%03d.mrc'%(c,))

        if c > c_best:
            # align c with the selected template
            rotate_key(self=self, vk=tk[c], vk_out=tka[c], angle=pa[c_best][c]['angle'], loc=pa[c_best][c]['loc'])
            print c_best, '-', c, ':', '%0.3f'%(pa[c_best][c]['score'],), N.linalg.norm( pa[c_best][c]['loc'] ),  '\t',
        else:
            # do not rotate, just make a copy
            shutil.copyfile(tk[c]['subtomogram'], tka[c]['subtomogram'])
            shutil.copyfile(tk[c]['mask'], tka[c]['mask'])
            print 'copy(%d)'%(c,), '\t',

    print       ;       sys.stdout.flush()

    return {'tka':tka, 'pa':pa}






'''
Perform pariwise alignment of templates, then sort the alignment scores from high to low. For score of each pair, say c0 and c1, where the c0 is a better cluster than c1. If the translation is larger than a threshold, ignore this alignemt. Otherwise, if c1 has not been aligned to any subtomogram, then align c1 with c0 and record this alignment. 

The main advantage of such approach is it can deal with the case of more than two clusters
'''

def cluster_average_align_common_frame__multi_pair(self, tk, align_op, loc_r_max, pass_dir):
    print 'align averages to common frames'
    #print 'cluster_average_align_common_frame__multi_pair()'

    out_dir = os.path.join(pass_dir, 'common_frame')
    if not os.path.isdir(out_dir):      os.mkdir(out_dir)

    # first perform pairwise alignment
    pa = cluster_average_align_common_frame__pairwise_alignment(self=self, template_keys=tk, align_op=align_op)

    
    # prepare a list of alignments
    pal = []
    for c0 in pa:
        for c1 in pa[c0]:
            if c0 >= c1: continue           # we assume a cluster with smaller id is better, c0 is better than c1
            if N.linalg.norm( pa[c0][c1]['loc'] ) > loc_r_max:        continue        # ignore those alignments with large shift
            pal.append(pa[c0][c1])

    # order the alignments according to decrease of alignment scores
    pal = sorted(pal, key=lambda _:(- _['score']))

    # determine selected alignments
    align_to_clus = {}
    unrotated_clus = set()
    for i in range(len(pal)):
        palt = pal[i]
        c0 = palt['c0']
        c1 = palt['c1']
        assert  c0 < c1
        if c1 in unrotated_clus: continue   # this means some template alredy needs to align to c1, c1 should not be rotated
        if c0 in align_to_clus:  continue   # this means c1 is already aligned to some template
        if c1 in align_to_clus:  continue        # this means c1 is already aligned to some template, it has already been processed
        unrotated_clus.add(c0)
        align_to_clus[c1] = c0

    assert  len(unrotated_clus.intersection(align_to_clus.keys())) == 0


    # produce alignments
    tka = {}
    for c in tk:
        tka[c] = {}
        tka[c]['id'] = c
        tka[c]['pass_i'] = tk[c]['pass_i']
        tka[c]['cluster'] = tk[c]['cluster']
        tka[c]['subtomogram'] = os.path.join(out_dir, 'clus_vol_avg_%03d.mrc'%(c,))
        tka[c]['mask'] = os.path.join(out_dir, 'clus_mask_avg_%03d.mrc'%(c,))

        if c in align_to_clus:
            # align c with the selected template
            rotate_key(self=self, vk=tk[c], vk_out=tka[c], angle=pa[align_to_clus[c]][c]['angle'], loc=pa[align_to_clus[c]][c]['loc'])
            print align_to_clus[c], '-', c, ':', '%0.3f'%(pa[align_to_clus[c]][c]['score'],), N.linalg.norm( pa[align_to_clus[c]][c]['loc'] ),  '\t',
        else:
            # do not rotate, just make a copy
            shutil.copyfile(tk[c]['subtomogram'], tka[c]['subtomogram'])
            shutil.copyfile(tk[c]['mask'], tka[c]['mask'])
            print 'copy(%d)'%(c,), '\t',

    print       ;       sys.stdout.flush()


    return {'tka':tka, 'unrotated_clus':unrotated_clus, 'align_to_clus':align_to_clus, 'pa':pa, 'pal':pal}



'''
parameters: c: cluster id,      tk: template key,       op: option
'''
def template_segmentation__single(c, tk, op):
    out_path = os.path.join(os.path.dirname(tk['subtomogram']), 'clus_vol_seg_phi_%03d.mrc'%(c,))

    v = IV.read_mrc(tk['subtomogram'])['value']

    phi = template_segmentation__single_vol(v=v, op=op)

    if phi is not None:
        tk['segmentation'] = out_path
        IV.put_mrc(phi, out_path, overwrite=True)


    return {'c':c, 'tk':tk}



def template_segmentation__single_vol(v, op):

    if not op['density_positive']:      v = -v
    del op['density_positive']

    if ('normalize_and_take_abs' in op) and op['normalize_and_take_abs']:
        v -= v.mean()
        v = N.abs(v)

    if 'gaussian_smooth_sigma' in op:
        vg = FG.smooth(v=v, sigma=float(op['gaussian_smooth_sigma']))
        del op['gaussian_smooth_sigma']
    else:
        vg = v

    phi = SACS.segment_with_postprocessing(vg, op)

    if phi is None:
        sys.stderr.write( 'Warning: segmentation failed ' + out_path + '\n' )

    if phi is not None:
        '''
        important!! sometimes, especially in the experimental data, a cluster average can still be very noisy, in such case, the none structural region may occationaly have a higher mean intensity than strucural region. We assume the structural region should have less overlap with the boundary of the subvolume than the none structural region does. We use this as a constrain to discard those bad quality segmentations
        '''
        bm = MU.boundary_mask(phi.shape)
        if bm[phi<0].sum() < bm[phi>0].sum():   
            phi = None
            sys.stderr.write( 'Warning: segmentation of the following cluster average violates boundary condition ' + out_path + '\n')

    return phi


# perform segmentation and store segmented masks
def template_segmentation(self, tk, op, multiprocessing=False):
    #print 'template_segmentation()', op

    if multiprocessing:
        if self.pool is None:       self.pool = Pool()
        
        pool_results = [self.pool.apply_async(func=template_segmentation__single, kwds={'c':c, 'tk':tk[c], 'op':op}) for c in tk]
        for r in pool_results:
            r = r.get(999999)
            tk[r['c']] = r['tk']        # get updated template keys, which includes segmented files

        self.pool.close()
        self.pool = None
    else:
        for c in tk:
            r = template_segmentation__single(c=c, tk=copy.deepcopy(tk[c]), op=copy.deepcopy(op))       # just make a copy to make sure that modifications do not affect original copy
            tk[r['c']] = r['tk']        # get updated template keys, which includes segmented files




'''
parameters: a subtomogram v and a template segment mask, they are aligned
'''
def template_guided_segmentation(v, m, op):

    op = copy.deepcopy(op)

    v_org = N.copy(v)           # keep a copy of original v

    if not op['density_positive']:        v = -v_org
    del op['density_positive']

    if 'gaussian_smooth_sigma' in op:
        vg = FG.smooth(v=v, sigma=float(op['gaussian_smooth_sigma']))
        del op['gaussian_smooth_sigma']
    else:
        vg = v



    #calculate the mean intensity of vg inside mask and outside mask, then fix these two mean value and apply chan-vese to obtain structural regions
    if  (m > 0.5).sum() == 0:   return
    if  (m < 0.5).sum() == 0:   return

    op['mean_values'] = [vg[m<0.5].mean(), vg[m>0.5].mean()]
    phi = SACS.segment_with_postprocessing(vg, op)

    if phi is None: return

    if (phi > 0).sum() == 0:    return
    if (phi < 0).sum() == 0:    return

    # get connected components of structural regions, and seperate the components into two groups: those have overlap with m and those does not
    struc = N.array((phi > 0), dtype=N.int32, order='F')        # structural regions
    mcr = core.connected_regions(  struc  )

    struc_tem = N.zeros(vg.shape)
    for l in range(1, mcr['max_lbl'] + 1):
        if N.logical_and((mcr['lbl'] == l), m).sum() > 0:
            struc_tem[mcr['lbl'] == l] = 1
        else:
            struc_tem[mcr['lbl'] == l] = 2

    if (struc_tem == 1).sum() == 0:     return

    # perform watershed on phi using these two groups as seed.
    sws = SW.segment(vol_map=phi, vol_lbl=struc_tem)

    seg = sws['vol_seg_lbl'] == 1
    if seg.sum() == 0:      return

    seg = N.logical_and((phi > (phi[seg].max() * op['phi_propotion_cutoff'])), seg)       # restrict to a margin according to phi

    if seg.sum() == 0:  return

    vs = N.zeros(v.shape) + N.nan
    vs[seg] = v_org[seg]

    return vs

'''
alternatively, can use such template phi as a WEIGHT when calculating the difference between mean and subtomogram intensities
'''



def align_to_templates__pair_align(c, t_key, v, vm, align_op):

    if align_op['with_missing_wedge']:
        t = IV.get_mrc(t_key['subtomogram'])
        tm = IV.get_mrc(t_key['mask'])
        at_re = align_vols_with_wedge(v1=t, m1=tm, v2=v, m2=vm, op=align_op)
    else:
        t = IV.get_mrc(t_key['subtomogram'])
        at_re = align_vols_no_wedge(v1=t, v2=vi, op=align_op)

    at_re['c'] = c

    return at_re


'''
given a subtomogram imputed using its old template, align it with new templates and get best alignment score
rec is one record in json format


mxu140311: I have tried to increased template_wedge_cutoff from 0.1 to 0.5. I often find there are wedge biases of aligning subtomograms against averages of small sized particles. There are multiple averages of the same particle with different missing wedge regions. If template_wedge_cutoff is low, some of the missing wedge regions will be included into alignment, which biases the alignment, and becomes an image feature that can be discriminated by cluster_removal_according_to_center_matching_specificity(). The cluster_removal_according_to_center_matching_specificity() does not effectively eliminate redundancies due to aligment bias. Result: ~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_specificity_analysis_with_alignment.py shows that template_wedge_cutoff=0.5 and template_wedge_cutoff=0.1 make only very very small difference!!!

'''

def align_to_templates(self, rec=None, segmentation_tg_op=None, tem_keys=None, template_wedge_cutoff=0.1, align_op=None, multiprocessing=False):
    #print 'align_to_templates()', 'rec', rec, 'segmentation_tg_op,', segmentation_tg_op, 'tem_keys', tem_keys, 'align_op', align_op         ;      sys.stdout.flush()

    vi = None    
    if align_op['with_missing_wedge']:
        v =  self.cache.get_mrc(rec['subtomogram'])
        vm =  self.cache.get_mrc(rec['mask'])
    else:
        raise Exception('following options are need to be doube checked')
        if 'template' not in rec:  rec['template'] = None
        vi = impute_vol_keys(vk=rec, ang=rec['angle'], loc=rec['loc'], tk=rec['template'], align_to_template=False, normalize=True, cache=self.cache)['vi']     # important: set align_to_template=False, so that subtomogram is not rotated, but template is rotated

    if (segmentation_tg_op is not None) and ('template' in rec) and ('segmentation' in rec['template']):
        v = align_to_templates__segment(rec=rec, v=v, segmentation_tg_op=segmentation_tg_op)['v']

    # record pairwise alignment results
    if multiprocessing:
        # Thread / process pool for parallel computing.
        
        if self.pool is None:       self.pool = Pool()

        #print "Active threads:", threading.active_count()
        pool_results = [self.pool.apply_async(func=align_to_templates__pair_align, kwds={'c':c, 't_key':tem_keys[c], 'v':v, 'vm':vm, 'align_op':align_op}) for c in tem_keys]
        align_re = {}
        for r in pool_results:
            at_re = r.get(999999)
            c = at_re['c']
            align_re[c] = at_re

            if N.isnan(align_re[c]['score']):
                if self.logger is not None:      self.logger.warning("alignment failed: rec %s, template %s, error %s ", repr(rec), repr(tem_keys[c]), repr(align_re[c]['err']) )

        self.pool.close()
        self.pool = None

    else:

        align_re = {}
        for c in tem_keys:

            if self.work_queue.done_tasks_contains(self.task.task_id):            raise Exception('Duplicated task')         # terminate alignment process when the corresponding task is already completed by another worker

            #print tem_keys[c];      sys.stdout.flush()

            align_re[c] = align_to_templates__pair_align(c=c, t_key=tem_keys[c], v=v, vm=vm, align_op=align_op)

            if N.isnan(align_re[c]['score']):
                if self.logger is not None:      self.logger.warning("alignment failed: rec %s, template %s, error %s ", repr(rec), repr(tem_keys[c]), repr(align_re[c]['err']) )

    return {'vol_key':rec, 'align':align_re}


def align_to_templates__segment(rec, v, segmentation_tg_op):
    #print 'apply segmentation before alignment', rec['template']['segmentation_tg'], segmentation_tg_op

    phi = IV.read_mrc_vol(rec['template']['segmentation'])
    phi_m = phi > 0.5             # why phi>0.5 instead of 0? Because we do not want to include boundary region?

    ang_inv, loc_inv = AAL.reverse_transform_ang_loc(rec['angle'], rec['loc'])
    phi_mr = GR.rotate(phi_m, angle=ang_inv, loc_r=loc_inv, default_val=0)

    v_s = template_guided_segmentation(v=v, m=phi_mr, op=segmentation_tg_op)

    if (v_s is not None) and (v_s[N.isfinite(v_s)].std() > 0):
        v = v_s;            del v_s
   
        v_t = N.zeros(v.shape) 
        v_f = N.isfinite(v)
        v_t[v_f] = (v[v_f] - v[v_f].mean()) / v[v_f].std()      # I suppose the normalization should not harm alignment, since alignment calculates correlation
        v = v_t;            del v_f, v_t

    #IV.put_mrc(v, os.path.join('/tmp', 'vol-seg-'+str(uuid.uuid1())+'.mrc'))    # for debugging inspection
    return {'v':v, 'phi_m':phi_m, 'phi_mr':phi_mr}


def align_to_templates__batch(self, op, data_json, segmentation_tg_op, tmp_dir, tem_keys):

    if ('template' in op) and ('match' in op['template']) and ('priority' in op['template']['match']):
        task_priority = op['template']['match']['priority']
    else:
        task_priority = 2000 + N.random.randint(100)        # add a random number so that this batch of tasks stay together

    print 'align against templates', 'segmentation_tg_op', (segmentation_tg_op if op['template']['match']['use_segmentation_mask'] else None), 'task priority', task_priority         ;       sys.stdout.flush()

    # if there are matching results stored in the tmp_dir, load them
    at_ress = []

    for f in os.listdir(tmp_dir):
        if not f.endswith('.pickle'):        continue
        res_file = os.path.join(tmp_dir, f)
        if not os.path.isfile(res_file):    continue
        with open(res_file, 'rb') as f:     at_ress_t = pickle.load(f)
        at_ress.append(at_ress_t)

    if len(at_ress) > 0:        print 'loaded previous', len(at_ress), ' resutlts'         ;       sys.stdout.flush()

    completed_subtomogram_set = set([_.result['vol_key']['subtomogram'] for _ in at_ress])

    tasks = []
    for rec in data_json:
        if rec['subtomogram'] in completed_subtomogram_set:     continue            # ignore those completed alignments

        tasks.append(self.runner.task( priority=task_priority, module='tomominer.pursuit.multi.util', method='align_to_templates', kwargs={'rec':rec, 'segmentation_tg_op':(segmentation_tg_op if op['template']['match']['use_segmentation_mask'] else None), 'tem_keys':tem_keys, 'align_op':op['align'], 'multiprocessing':False} ))

    for at_ress_t in self.runner.run__except(tasks):
        at_ress.append(at_ress_t)
        res_file = os.path.join(tmp_dir, '%s.pickle'%(at_ress_t.task_id))
        with open(res_file, 'wb') as f:     pickle.dump(at_ress_t, f, protocol=0)           # use protocol 0 (ASCII) for the ease of human inspection

    return at_ress



# given a set of subtomogram keys in tk, find corresponding cluster in cluster_info ci
def cluster_info_filter_by_subtomogram_key(tk, ci):
    ci_s = {}
    for k in tk:
        ci_selected = None
        for i in ci:
            for c in ci[i]:
                if str(ci[i][c]['template_key']['subtomogram']) != tk:      continue
                assert      ci_selected is None
                ci_selected = ci[i][c]
                break

        assert ci_selected is not None
        ci_s[k] = ci_selected
    return ci_s



# given template alignment and best match information, sequencially include subtomograms of best alignment and sequentially calculate FSC scores. Then look at local maxima of the FSC sequence, and choose the one that has best correlation with the template
# see idea 1.1 of:                https://docs.google.com/document/d/1NmYhEczJ8A--KI5PStf8hyuPM5QPfRLz7RftOfkn8_k/edit
# parameters:       dj: data_json       ci: cluster_info used to collect the corresponding FSC score of original templates. If ci is not None, we require the FSC score of refined clusters to be larger than original cluster.
# return:           selected data_json, and corresponding cluster labels
def cluster_formation_alignment_fsc__by_local_maxima(self, dj, ci=None, op=None):
    print 'cluster_formation_alignment_fsc__by_local_maxima()'

    if 'debug' not in op:       op['debug'] = False

    dj = copy.deepcopy(dj)      # IMPORTANT: must make a copy of dj, otherwise the following generated volumeric data will be inserted into the original dj, and eat up storage in following processings

    # collect the records that maximally match each template
    djm = defaultdict(list)
    for r in dj:
        if 'template' not in r:     continue
        djm[ str(r['template']['subtomogram']) ].append(r)

    # collect templates information
    ts = {}     # collection of templates (and masks)
    for r in dj:
        if 'template' not in r:     continue

        tk = r['template']
        k = str(tk['subtomogram'])
        if k in ts:    continue

        ts[k] = tk
        ts[k]['v'] = IV.read_mrc_vol(str(tk['subtomogram']))
        ts[k]['m'] = IV.read_mrc_vol(str(tk['mask']))
        ts[k]['fv'] = AFGF.foster_transform_single_vol(v=ts[k]['v'], m=ts[k]['m'])

        if ci is not None:          ts[k]['fsc'] = ci[ts[k]['pass_i']][ts[k]['cluster']]['fsc'].sum()
 

    # order subtomograms according to decrease of alignment scores, and restrict the maximum length of each list
    for k in djm:
        djmt = djm[k]
        djmt = sorted( djmt, key=lambda _ : float(_['score']), reverse=True )
        if ('max_expansion_size' in op) and (len(djmt) > op['max_expansion_size']):
            djmt = djmt[:op['max_expansion_size']]

        djm[k] = djmt
    
    # parallel calculate sequential ssnr
    ssnr_sequential_op = copy.deepcopy(op['ssnr_sequential'])
    ssnr_sequential_op['n_chunk'] = op['n_chunk']
    ssnr_s = SS.ssnr_sequential_parallel(self=self, data_json_dict=djm, op=ssnr_sequential_op)



    # retrive FSC score sequence
    fsc_s = {}
    for k in ssnr_s:        fsc_s[k] = N.array([N.sum(_) for _ in ssnr_s[k]['fsc']])

    import scipy.ndimage.filters as SDF
    if 'gaussian_smooth_sigma' in op:
        for k in fsc_s:        fsc_s[k] = SDF.gaussian_filter1d(fsc_s[k], op['gaussian_smooth_sigma'])           # sometimes there are too many local maxima, can smooth a bit to reduce


    # find local maxima and collect corresponding data_json
    import tomominer.filter.local_extrema as FL
    dj_lm = {}
    c = 0       # temporary cluster label
    for k in fsc_s:

        inds = FL.local_maxima(fsc_s[k])
        inds = inds[0]
 
        if ('min_expansion_size' in op):        inds = inds[inds > op['min_expansion_size']]        # remove those local maxima that corresponds to small sets
        if 'fsc' in ts[k]:      inds = [_ for _ in inds if (fsc_s[k][_] > ts[k]['fsc'])]         # we require the refined set average to have a higher FSC score than the original average, in such case, 

        if len(inds) == 0:  continue

        if op['debug']:        print 'template', k, 'number of local maxima', len(inds)
        
        for i in inds:
            dj_lm[c] = {'c':c, 'k':k, 'i':i, 'data_json':copy.deepcopy(djm[k][:(i+1)]), 'fsc':fsc_s[k][i]}
            c += 1

    if op['debug']:
        with open('/tmp/cluster_formation_alignment_fsc__by_local_maxima.pickle', 'wb') as f:       pickle.dump({'ssnr_s':ssnr_s, 'fsc_s':fsc_s, 'dj_lm':dj_lm}, f, protocol=-1)

    if len(dj_lm) == 0:     return {'best':{}, 'dj_lm':dj_lm, 'djm':djm}

    if op['debug']:        print '# parallel averaging (without storing the averages)'
    averaging_op = copy.deepcopy(op['averaging'])
    averaging_op['n_chunk'] = op['n_chunk']
    cav = cluster_averaging_vols(self=self, clusters={_:dj_lm[_]['data_json'] for _ in dj_lm}, op=averaging_op)
    cav = cav['cluster_avg_dict']

    for c in dj_lm:        dj_lm[c]['average'] = cav[c]


    # calculate correlation between original templates and new averages
    for c in dj_lm:
        vf_maximum = AFGF.foster_transform_single_vol(v=dj_lm[c]['average']['vol'], m=dj_lm[c]['average']['mask'])
        dj_lm[c]['score'] = N.sum(ts[k]['fv'] * vf_maximum)        # correlation score

    # select those achieve maximum correlations, record results
    dj_lm_best = {}
    for k in djm:
        best = None
        for c in dj_lm:
            if dj_lm[c]['k'] != k:     continue
            if (best is None) or (best['score'] < dj_lm[c]['score']):                best = dj_lm[c]
        if best is None:        continue
        dj_lm_best[k] = best

        print k, 'best match size', len(djm[k]), 'new set size', len(dj_lm_best[k]['data_json'])


    return {'best':dj_lm_best, 'dj_lm':dj_lm, 'djm':djm}



    
# given template alignment and best match information, sequencially include subtomograms of best alignment and sequentially calculate FSC scores. Then look at global maximum of the FSC sequencev# parameters:       dj: data_json
# return:           selected data_json, and corresponding cluster labels
def cluster_formation_alignment_fsc__by_global_maximum(self, dj, op=None):
    #print 'cluster_formation_alignment_fsc__by_global_maximum()'

    if 'debug' not in op:       op['debug'] = False

    dj = copy.deepcopy(dj)      # IMPORTANT: must make a copy of dj, otherwise the following generated volumeric data will be inserted into the original dj, and eat up storage in following processings

    # collect the records that maximally match each template
    djm = defaultdict(list)
    for r in dj:
        if 'template' not in r:     continue
        djm[ str(r['template']['subtomogram']) ].append(r)

    djm_org = copy.deepcopy(djm)

    # order subtomograms according to decrease of alignment scores, and restrict the maximum length of each list
    for k in djm:
        djmt = djm[k]
        djmt = sorted( djmt, key=lambda _ : float(_['score']), reverse=True )
        if ('max_expansion_size' in op) and (len(djmt) > op['max_expansion_size']):
            djmt = djmt[:op['max_expansion_size']]

        djm[k] = djmt
    
    # parallel calculate sequential ssnr
    ssnr_sequential_op = copy.deepcopy(op['ssnr_sequential'])
    ssnr_sequential_op['n_chunk'] = op['n_chunk']
    ssnr_s = SS.ssnr_sequential_parallel(self=self, data_json_dict=djm, op=ssnr_sequential_op)



    # retrive FSC score sequence
    fsc_sum = {}
    for k in ssnr_s:        fsc_sum[k] = N.array([N.sum(_) for _ in ssnr_s[k]['fsc']])

    import scipy.ndimage.filters as SDF
    if 'gaussian_smooth_sigma' in op:
        for k in fsc_sum:        fsc_sum[k] = SDF.gaussian_filter1d(fsc_sum[k], op['gaussian_smooth_sigma'])           # sometimes there are too many local maxima, can smooth a bit to reduce

    if 'min_expansion_size' in op:
        for k in copy.deepcopy(fsc_sum.keys()):
            if len(fsc_sum[k]) < op['min_expansion_size']:
                del fsc_sum[k]
                continue
            fsc_sum[k][:op['min_expansion_size']] = N.min(fsc_sum[k])-1              # remove those maximum that corresponds to small sets


    # find global maximum and collect corresponding data_json
    dj_gm = {}
    for k in fsc_sum:

        i = N.argmax(fsc_sum[k])

        if op['debug']:        print 'template', k, 'original subtomogram num', len(djm_org[k]), 'global maximum', i
        
        dj_gm[k] = {'k':k, 'i':i, 'data_json':copy.deepcopy(djm[k][:(i+1)]), 'fsc':ssnr_s[k]['fsc'][i], 'fsc_sum':fsc_sum[k][i]}


    return {'dj_gm':dj_gm, 'djm':djm, 'ssnr_s': ssnr_s}


 
'''
first perform pursuit based on kmeans mode, then cut the hierarchical clusters so that the result grouping is most consistent with hierarchical clustering.

see
https://docs.google.com/document/d/1xU1AWj--Uol-2fx_xhgToneeRjmOEtgAOW04_PF_uCU/edit
~/ln/tomominer/tomominer/pursuit/multi/tests/clustering/hierarchy/cut_according_to_last_selected_clusters.py

parameters:     dj: data_json,      cas: cluster selection result,       cis: cluster_info_stat,     hc_re: hierarchical clustering result
'''
def hierarchical_clustering__optimal_consistent_dist_cutoff(dj, cas, cis, hc_re, verbose=False):
    print 'hierarchical_clustering__optimal_consistent_dist_cutoff()'

    dj_sel = {}
    for i in cas['selected_templates']:
        k = cas['selected_templates'][i]
        if cis[k['pass_i']][k['cluster']]['is_specific'] is not None:       continue
        dj_sel[i] = cas['tk_info'][k['subtomogram']]['data_json']

    if verbose:        print 'size distribution of selected clusters', sorted([len(dj_sel[_]) for _ in dj_sel], reverse=True)      # size distribution


    dj_sel_dict = {}            # the dictionary of subtomogram's set id, notice that it may have a little ambiguity when we allow a little bit of overlap between sets in the selection process
    for i in dj_sel:
        for _ in dj_sel[i]:
            dj_sel_dict[_['subtomogram']] = i

    label_sel = [N.nan] * len(dj)
    for i, _ in enumerate(dj):
        if _['subtomogram'] not in dj_sel_dict:     continue
        label_sel[i] = dj_sel_dict[_['subtomogram']]


    if verbose:        print 'select cutoff for most consistent clustering'

    import tomominer.cluster.hierarchy.hierarchy as CHH
    best_cut = CHH.optimal_consistent_cut(z=hc_re['clusters'], l0=label_sel, find_min_dist=True, consistency_matric='nmi', verbose=verbose)

    if verbose:
        # inspect cluster size distribution
        chh_cc = CHH.cluster_cut(r=hc_re['root_id'], hi=hc_re['info'], dist_cutoff=best_cut['t'])
        print 'cluster size distribution after hierarchy cut', sorted([hc_re['info'][_]['size'] for _ in chh_cc], reverse=True)


    return best_cut




def rigid_transform_difference(data_json, data_json_new):

    ang = {}
    loc = {}
    for r in data_json:
        ang[r['subtomogram']] = N.array(r['angle'])
        loc[r['subtomogram']] = N.array(r['loc'])


    dif = {}
    for r in data_json_new:     dif[r['subtomogram']] = AAL.angle_zyz_translation_difference(N.array(r['angle']), N.array(r['loc']), ang[r['subtomogram']], loc[r['subtomogram']])


    dif_mean = N.array([dif[i] for i in dif]).mean()

    return dif_mean




