
# functions for dimension reduction
import os
import pickle
import time
import random

import numpy as np
from numpy.fft import fftn, ifftn, fftshift, ifftshift

import tomominer.io.file as iv

import tomominer.average.util as avgu
import tomominer.image.vol.util as uv



# mxu: original parallel PCA dimension reduction and k-means clustering written by Zach, moved here from classify.py
def standard_pca_clustering_parallel(data, self, worker, pass_dir, n_chunk=100, dim_reduction_dims=-1, cluster_kmeans_k=-1, avg_key=None):


    # parallel dimension reduction and k-means clustering

    # TODO: we are double computing everything.  Try this strategy:
    #
    # send coordinate to worker, and get back sub matrix directly.  Then we
    # can assign opposite coordinate the transpose of the data we recieved.
    # Only send out half the work.
    # 
    # [ (0,0), (0,1), (0,2) ]
    # [      , (1,1), (1,2) ]
    # [      ,      , (2,2) ]
    # 
    # We send out only the blocks pictured.
    # The we have a rule: when we get back (i,j), if i != j, then also assign to (j,i)
    #


    start_time = time.time()
    # Calculate matrix elements.
    tasks = []
    cnt = 0
    mat_structure = []

    vmal1 = [_ for _ in data]
    i = 0
    while vmal1:
        part1 = vmal1[:n_chunk]
        vmal1 = vmal1[n_chunk:]

        mat_row = []
        j = 0
        vmal2 = [_ for _ in data]
        while vmal2:
            part2 = vmal2[:n_chunk]
            vmal2 = vmal2[n_chunk:]

            out_path = pass_dir + '/' + 'dim_reduce_submatrix_%d_%d.npy' % (i, j)
            mat_row.append(out_path)
            tasks.append(self.runner.task('dim_reduce_submatrix', part1, part2, avg_key, out_path))
            j+=1

        i+=1
        mat_structure.append(mat_row)

    for res in self.runner.run(tasks):
        #print "got res:", res
        if res.error:
            print "Computation Failed!"
            print res
            print "task_id      :", res.task_id
            print "method       :", res.method
            print "args         :", res.args
            print "kwargs       :", res.kwargs
            print "error_msg    :", res.error_msg
            raise Exception
    print "Build Correlation sub-matrix: %2.6f sec" % (time.time() - start_time)

    # combine matrix from pieces and do dimension reduction.
    # also do kmeans.

    if False:
        start_time = time.time()
        tasks = [ self.runner.task('dimension_reduction', mat_structure, opt.dim_reduction_dims, opt.cluster_kmeans_k)]
        #wait = [_ for _ in self.runner.run(tasks)]
        for res in self.runner.run(tasks):
            if res.error:
                print "Computation Failed!"
                print res
                print "task_id      :", res.task_id
                print "method       :", res.method
                print "args         :", res.args
                print "kwargs       :", res.kwargs
                print "error_msg    :", res.error_msg
                raise Exception
            labels = res.result
    else:
        labels = dimension_reduction(mat_structure, dim_reduction_dims, cluster_kmeans_k)       # the dimension reduction is very memory consuming, run locally 

    print "Building matrix, Dimension Reduction, k-means: %2.6f sec" % (time.time() - start_time)

    return {'labels':labels}
   
# mxu: calculate missing wedge masked difference (according to the dimension reduction paper), and use NIPALS to perform dimension reduction
def nipals_clustering(op):

    data = op['data']
    runner = op['runner']

    vmal = [_ for _ in data]
    tasks = []
    while vmal:
        tasks.append(self.runner.task('masked_vol_difference', op['avg_key'], vmal[:n_chunk]))
        vmal = vmal[n_chunk:]

    res = [_ for _ in self.runner.run__except(tasks)]

    # merge all results
    vol_dif_d = {}
    for res_t in res:
        vol_dif_d = dict(vol_dif_d.items() + res_t.result.items())

    labels = worker.dimension_reduction_nipals([vol_dif_d[_[0]] for _ in data], opt.dim_reduction_dims, opt.cluster_kmeans_k)

# mxu: calculate missing wedge masked difference (according to the dimension reduction paper), and use randomized PCA to perform dimension reduction
def rand_pca_clustering(op):

    data = op['data']
    runner = op['runner']


    vmal = [_ for _ in data]
    tasks = []
    while vmal:
        tasks.append(self.runner.task('masked_vol_difference', op['avg_key'], vmal[:['n_chunk']]))
        vmal = vmal[op['n_chunk']:]

    res = [_ for _ in self.runner.run__except(tasks)]

    # merge all results
    vol_dif_d = {}
    for res_t in res:
        vol_dif_d = dict(vol_dif_d.items() + res_t.result.items())

    labels = worker.dimension_reduction_randomized_pca([vol_dif_d[_[0]] for _ in data], opt.dim_reduction_dims, opt.cluster_kmeans_k)



# mxu: parallel calculate missing wedge masked difference (according to the dimension reduction paper), and use randomized PCA to perform dimension reduction
# WARNING: op['dims'] cannot be too large, say <= 20 is safer. Otherwise the system will break due to heavy data transfer in queue
def pca_parallel(self, data, dims, avg_key, n_chunk, sample_num=None, pca_op=None):   

    if (sample_num==None) or (len(data) <= 2*sample_num):
        mat = pca__stack_dif__parallel(self=self, avg_key=avg_key, data=data, n_chunk=n_chunk)
        trn_res = pca__train(dims=dims, mat=mat, transform=True, pca_op=pca_op)      # we can also choose   transform=True,  the transform takes similiar amount of time as the training step, therefore we still make the transform parallel
        red = trn_res['red']
    else:
        data_spl = random.sample(data, sample_num)

        mat = pca__stack_dif__parallel(self=self, avg_key=avg_key, data=data_spl, n_chunk=n_chunk)
        trn_res = pca__train(dims=dims, mat=mat, pca_op=pca_op)        

        # parallel dimension reduction.     When (len(data) <= 2*sample_num), we can also choose pca__train
        start_time = time.time()

        pca = trn_res['pca']
        vmal = [_ for _ in data]
        inds = range(len(data))

        tasks = []
        while vmal:
            vmal_t = vmal[:n_chunk]
            inds_t = inds[:n_chunk]
            tasks.append(self.runner.task('pca__transform', pca, avg_key, vmal_t, inds_t))
            vmal = vmal[n_chunk:]
            inds = inds[n_chunk:]

        res = [_ for _ in self.runner.run__except(tasks)]

        # merge all results
        red = np.zeros([len(data), dims])
        for res_t in res:
            red[res_t.result['inds'],:] = res_t.result['red']


        print "PCA transform : %2.6f sec" % (time.time() - start_time)


    return red





# mxu: compute masked difference between subtomgrams and the global average
def masked_vol_difference(avg_key, vmal_in, tmp_dir=None):
    avg = avgu.load_avg(avg_key)

    v_dif = {}
    for vk, mk, ang, loc in vmal_in:
        # load vol/mask.
        vol  = iv.get_mrc_cache_fs(vk, tmp_dir)
        mask = iv.get_mrc_cache_fs(mk, tmp_dir)

        v_r_msk_dif, m_r = masked_difference_given_vol_avg_fft(v=vol, m=mask, ang=ang, loc=loc, vol_avg_fft=avg['vol_avg_fft'], vol_mask_avg=avg['vol_mask_avg'])

        v_dif[vk] = v_r_msk_dif

    return v_dif    



# mxu: IMPORTANT: correction of dimension reduction procedure according to the paper  Heumann12
def masked_difference_given_vol_avg_fft(v, m, ang, loc, vol_avg_fft, vol_mask_avg, smoothing_gauss_sigma=None, mean_standardize=False, std_standardize=False):
    import tomominer.core as tomo

    v_r = tomo.rotate_vol_pad_mean_py(v,ang,loc)
    m_r = tomo.rotate_mask_py(m, ang)


    v_r_fft = fftshift(fftn(v_r))
    v_r_msk_dif = np.real(ifftn(ifftshift( (v_r_fft - vol_avg_fft) * m_r * vol_mask_avg )))

    if smoothing_gauss_sigma > 0:
        # when the subtomograms are very noisy, PCA alone may not work well. Because PCA only process vectors but does consider spetial structures. In such case we can make use spatial information as an additional layer of constrain (or approsimation model) by apply certain amount of gaussian smoothing
        import tomominer.filter.gaussian as fg
        v_r_msk_dif = fg.smooth(v=v_r_msk_dif, sigma=smoothing_gauss_sigma)
 
 
    if mean_standardize:
        # mxu: actually I do not know why following standardizations are needed, according to the paper
        v_r_msk_dif -= np.mean(v_r_msk_dif)
        if std_standardize:     v_r_msk_dif /= np.std(v_r_msk_dif)      # mxu: standardization added according to the paper. See quote: "More recently, we have stopped adjusting variance", therefore std_standardize=False by default
    
    return (v_r_msk_dif, m_r)




def dim_reduce_submatrix(vmal_1_in, vmal_2_in, avg_key, mat_out_key, tmp_dir=None, vol_mask_avg_threshold=0.1):

    mat = np.zeros( (len(vmal_1_in), len(vmal_2_in)), order='F')

    avg = avgu.load_avg(avg_key)

    # load volumes.
    vols1 = []
    masks1 = []

    # compute normalized volume for all entries in first volume list.
    for vk,mk,ang,loc in vmal_1_in:

        v = iv.get_mrc_cache_fs(vk, tmp_dir)
        m = iv.get_mrc_cache_fs(mk, tmp_dir)

        vrd, mr = masked_difference_given_vol_avg_fft(v=vol, m=mask, ang=ang, loc=loc, vol_avg_fft=avg['vol_avg_fft'], vol_mask_avg=avg['vol_mask_avg'])

        vols1.append(vrd)
        masks1.append(mr)


    vols2 = []
    masks2 = []
    # compute normalized volume for all entries in second volume list.
    for vk,mk,ang,loc in vmal_2_in:

        v = iv.get_mrc_cache_fs(vk, tmp_dir)
        m = iv.get_mrc_cache_fs(mk, tmp_dir)

        vrd, mr = masked_difference_given_vol_avg_fft(v=vol, m=mask, ang=ang, loc=loc, vol_avg_fft=vol_avg_fft, vol_mask_avg=vol_mask_avg)

        vols2.append(vrd)
        masks2.append(mr)

    # calculate pairwise scoring. which is simple dot product for normalized values.
    for i, v1 in enumerate(vols1):
        for j, v2 in enumerate(vols2):

            # mxu: assertion added
            cor = np.sum(v1*v2)
            if not np.isfinite(cor):
                self.logger.error("Exception: %s", traceback.format_exc())
                return False
            mat[i,j] = cor

    # save output back to disk.
    iv.np_save(mat_out_key, mat)

    return True

def dimension_reduction(mat_structure, dims, k):

    from scipy.cluster.vq import kmeans, vq

    mat = np.vstack( np.hstack(iv.np_load(_) for _ in row) for row in mat_structure )

    # TODO: figure out a way to do this smarter.  matrix is symmetric.
    mat = np.asfortranarray(mat)

    if mat.shape[1] < dims:     # mxu: bug correction, was mat.shape[0] < dims
        red = mat
    else:
        #red = tomo.dimension_reduction(mat, dims)
        red = dimension_reduction_numpy(mat, dims)


    centroids, _ = kmeans(red, k)
    labels,    _ = vq(red, centroids)

    return labels


# mxu: converted the dimension reduction from the cpp version to numpy version
def dimension_reduction_numpy(mat, dims):
    if dims > mat.shape[0]:
        self.logger.error("dimension reduction called with frewer dimensions then rows.")
        raise

    eig_val, eig_vec = np.linalg.eigh(mat)
    #print eig_val
    
    # sort descending
    S = abs(eig_val)
    idx = np.argsort( - S )     # use minus sign for descending order sort
    S_sqrt = np.sqrt(S)

    V = np.zeros([eig_vec.shape[0], dims])
    for dim_i in range(dims):
        V[:, dim_i] = eig_vec[:, idx[dim_i]] * S_sqrt[idx[dim_i]]
        #print eig_val[idx[dim_i]]
    
    return V    


   

def dimension_reduction_test():
    data_dir = '/home/rcf-47/mxu/tmp-shared/classification/pass_001'
    block_num = 20

    mat_structure = []
    for i in range(block_num):
        mat_row = []
        for j in range(block_num):
            mat_row.append('%s/dim_reduce_submatrix_%d_%d.npy'%(data_dir, i, j))
        mat_structure.append(mat_row)
    

    print np.vstack( np.hstack(iv.np_load(_) for _ in row) for row in mat_structure )

    dimension_reduction(mat_structure, 30, 5)
            
            
# use non-linear iterative partial least squares to perform dimension reduction
# then use k-means to perform clustering
# input: vs is a list of subvolums (masked difference between a subtomogram and global average, according to the dimension reduction paper)
def dimension_reduction_nipals(vs, dims, k):

    mat = np.zeros([len(vs), vs[0].size])
    for i in range(len(vs)):    mat[i,:] = vs[i].flatten()


    import mdp
    pca = mdp.nodes.NIPALSNode(output_dim=dims, conv=1E-2, max_it=1000)
    #pca.train(mat)
    red = pca.execute(mat, n=dims)


    from scipy.cluster.vq import kmeans, vq

    centroids, _ = kmeans(red, k)
    labels,    _ = vq(red, centroids)

    return labels

# calculate masked differences between subtomograms and gloval average, and form a matrix stacking all differences
def pca__stack_dif(avg_key, vmal_in, inds, tmp_dir=None, smoothing_gauss_sigma=None, voxel_mask_inds=None):
    
    avg = avgu.load_avg(avg_key)

    mat = None
    for i, (vk, mk, ang, loc) in enumerate(vmal_in):
        # load vol/mask.
        vol  = iv.get_mrc_cache_fs(vk, tmp_dir)
        mask = iv.get_mrc_cache_fs(mk, tmp_dir)

        v_dif, m_r = masked_difference_given_vol_avg_fft(vol, mask, ang, loc, avg['vol_avg_fft'], avg['vol_mask_avg'], smoothing_gauss_sigma=smoothing_gauss_sigma)
        v_dif = v_dif.flatten()

        if voxel_mask_inds is not None:
             v_dif =  v_dif[voxel_mask_inds]

        if mat is None:
            mat = np.zeros([len(vmal_in), v_dif.size])
        
        mat[i, :] = v_dif    

    return {'mat':mat, 'inds':inds}

def pca__stack_dif__parallel(self, avg_key, data, n_chunk, smoothing_gauss_sigma=None, voxel_mask_inds=None):

    start_time = time.time()

    vmal = [_ for _ in data]
    inds = range(len(data))

    tasks = []
    while vmal:
        vmal_t = vmal[:n_chunk]
        inds_t = inds[:n_chunk]
        tasks.append(self.runner.task('pca__stack_dif', avg_key=avg_key, vmal_in=vmal_t, smoothing_gauss_sigma=smoothing_gauss_sigma, inds=inds_t, voxel_mask_inds=voxel_mask_inds))
        vmal = vmal[n_chunk:]
        inds = inds[n_chunk:]

    # merge all results
    red = None
    for res_t in self.runner.run__except(tasks):
        if red is None:     red = np.zeros( [len(data), res_t.result['mat'].shape[1]] )

        red[res_t.result['inds'],:] = res_t.result['mat']


    print "Calculated matrix of masked difference to global : %2.6f sec" % (time.time() - start_time)

    return red



def pca__train(dims, mat, transform=False, pca_op=None):      # just set a large pca_job_num to make use all cpu cores


    start_time = time.time()

    if pca_op['mode'] == 0:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=dims)
    elif pca_op['mode'] == 1:
        from sklearn.decomposition import RandomizedPCA
        pca = RandomizedPCA(n_components=dims)
    elif pca_op['mode'] == 2:
        import multiprocessing
        from sklearn.decomposition import SparsePCA
        pca = SparsePCA(n_components=dims, alpha=pca_op['alpha'], n_jobs=min(multiprocessing.cpu_count(), pca_op['job_num']), verbose=True)    

    pca.fit(mat)
    print "PCA training : %2.6f sec" % (time.time() - start_time)

    if transform:
        start_time = time.time()
        red = pca.transform(mat)        # mxu: should use transform(), not fit_transform()!!
        print "RandomizedPCA transform : %2.6f sec" % (time.time() - start_time)
        return {'pca':pca, 'red':red}
    else:
        return {'pca':pca}    

# parallel transform all subtomograms (indexed by inds according order of data) given trained pca
def pca__transform(pca, avg_key, vmal_in, inds):
    assert len(vmal_in) == len(inds)

    pst_re = pca__stack_dif(avg_key=avg_key, vmal_in=vmal_in, inds=inds)
    mat = pst_re['mat']
    assert mat.shape[0] == len(inds)

    red = pca.transform(mat)        # mxu: should use transform(), not fit_transform()!!
    assert red.shape[0] == mat.shape[0]
    assert red.shape[1] == pca.get_params()['n_components']

    return {'inds':inds, 'red':red}



# use randomized PCA for dimension reduction
# then use k-means to perform clustering
# input: vs is a list of subvolums (masked difference between a subtomogram and global average, according to the dimension reduction paper)
def randomized_pca(vs, dims, k):

    mat = np.zeros([len(vs), vs[0].size])
    for i in range(len(vs)):    mat[i,:] = vs[i].flatten()


    from sklearn.decomposition import RandomizedPCA
    pca = RandomizedPCA(n_components=dims)
    pca.fit(mat)
    red = pca.transform(mat)


    from scipy.cluster.vq import kmeans, vq

    centroids, _ = kmeans(red, k)
    labels,    _ = vq(red, centroids)

    return labels



# calculate the covariance between neighbor voxels, then take average
def neighbor_covariance_avg__parallel(self, avg_key, data, n_chunk, smoothing_gauss_sigma=None):

    start_time = time.time()

    vmal = [_ for _ in data]
    inds = range(len(data))

    tasks = []
    while vmal:
        vmal_part = vmal[:n_chunk]
        inds_t = inds[:n_chunk]
        tasks.append(self.runner.task('neighbor_covariance__collect_info', avg_key, vmal_part, smoothing_gauss_sigma=smoothing_gauss_sigma))
        vmal = vmal[n_chunk:]
        inds = inds[n_chunk:]


    sum_global = None
    neighbor_prod_sum = None
    for res_t in self.runner.run__except(tasks):
        if sum_global is None:
            sum_global = res_t.result['sum']
        else:
            sum_global += res_t.result['sum']
        
        if neighbor_prod_sum is None:
            neighbor_prod_sum = res_t.result['neighbor_prod_sum']    
        else:
            neighbor_prod_sum += res_t.result['neighbor_prod_sum']    
    
    avg_global = sum_global / len(data)
    neighbor_prod_avg = neighbor_prod_sum / len(data)


    shift = res_t.result['shift']
    cov = np.zeros(neighbor_prod_avg.shape)
    for i in range(shift.shape[0]):
        cov[:,:,:,i] = neighbor_prod_avg[:,:,:,i] - avg_global * uv.roll(avg_global, shift[i,0], shift[i,1], shift[i,2])

    cov_avg = np.mean(cov, axis=3)

    print "Calculated neighbor covariance : %2.6f sec" % (time.time() - start_time)

    return cov_avg


# collecting information for calculating the neighbor covariance, calculated at worker side
def neighbor_covariance__collect_info(avg_key, vmal, smoothing_gauss_sigma=None, tmp_dir=None):

    avg = avgu.load_avg(avg_key)    

    sum_local = np.zeros(avg['vol_avg'].shape)
    neighbor_prod_sum = None

    for vk, mk, ang, loc in vmal:
        v  = iv.get_mrc_cache_fs(vk, tmp_dir)
        m = iv.get_mrc_cache_fs(mk, tmp_dir)
               
        v_r_msk_dif, m_r = masked_difference_given_vol_avg_fft(v, m, ang, loc, avg['vol_avg_fft'], avg['vol_mask_avg'], smoothing_gauss_sigma=smoothing_gauss_sigma)

       
        sum_local += v_r_msk_dif

        nei_prod = neighbor_product(v_r_msk_dif)
        if neighbor_prod_sum is None:
            neighbor_prod_sum = nei_prod['p']
        else:
            neighbor_prod_sum += nei_prod['p']

    return {'sum':sum_local, 'neighbor_prod_sum':neighbor_prod_sum, 'shift':nei_prod['shift']}

# calculate product of one voxel and all its neighbors
def neighbor_product(v):

    siz = list(v.shape)
    siz.append(26)

    p = np.zeros(siz, dtype=np.float32)    # use float32 to reduce data storage,
    shift = np.zeros( (26, 3), dtype=np.int8 )

    i = 0
    for s0 in range(-1,2):
        for s1 in range(-1, 2):
            for s2 in range(-1, 2):
                if (s0 == 0) and (s1 == 0) and (s2 == 0):     continue
                
                p[:,:,:, i] = v * uv.roll(v, s0, s1, s2)
                shift[i,:] = np.array([s0, s1, s2])

                i += 1 


    assert i == 26

    return {'p':p, 'shift':shift}


# calculate average covariance between neighbor voxels, then gaussian smooth and segment to identify a small amount of voxels as features for PCA analysis
# related paper: Sparse PCA via Covariance Thresholding
def covariance_filtered_pca(self, avg_key, data, n_chunk, diff_smoothing_gauss_sigma=None, gauss_sigma=None, cov_avg_min_cutoff_ratio=0.0, cutoff_num=10, max_feature_num=None, dims=None, do_segmentation=False, pca_op=None, out_dir=None):

    start_time = time.time()


    cov_avg_for_masking__cutoff = None

    # calculate average covariance, or load existing calculated result
    cov_avg = None
    cov_avg_g = None
    cov_avg_for_masking = None
    cov_avg_for_masking__cutoff = None

    scores = None
    score_best = None
    score_best_cutoff = None

    voxel_mask_inds = None

    if ((max_feature_num is not None) and (max_feature_num > 0)) or do_segmentation:        # in this case, need to calculated neighbor voxel covariance

        cov_avg_file = None
        if out_dir is not None:
            cov_avg_file = os.path.join(out_dir, 'cov_avg.pickle')
            if os.path.exists(cov_avg_file):
                with open(cov_avg_file, 'rb') as f:  cov_avg = pickle.load(f)

        if cov_avg is None:
            cov_avg = neighbor_covariance_avg__parallel(self=self, avg_key=avg_key, data=data, n_chunk=n_chunk, smoothing_gauss_sigma=diff_smoothing_gauss_sigma)
            if cov_avg_file is not None:
                with open(cov_avg_file, 'wb') as f: pickle.dump(cov_avg, f)

        
        # perform certain smoothing
        if gauss_sigma is None:
            cov_avg_for_masking = cov_avg
        else:
            import tomominer.filter.gaussian as fg
            cov_avg_g = fg.smooth(v=cov_avg, sigma=gauss_sigma)
            cov_avg_for_masking = cov_avg_g

        assert      cov_avg_for_masking.max() > 0

        # restrict the max number of featrues to be processed to be less than max_feature_num, in order to avoid the dimension reduction to be over time consuming
        cov_avg_for_masking_i = np.argsort((-cov_avg_for_masking), axis=None)
        cov_avg_for_masking__feature_num_cutoff = cov_avg_for_masking.flatten()[cov_avg_for_masking_i[min(max_feature_num, cov_avg_for_masking_i.size-1)]]
        cov_avg_for_masking__cutoff = max(cov_avg_for_masking__cutoff, cov_avg_for_masking__feature_num_cutoff)


        if do_segmentation:
            # segment out a region of high covariance
            import tomominer.segmentation.thresholding as st
            cutoffs = np.linspace(start=max(cov_avg_min_cutoff_ratio*cov_avg_for_masking.max(), cov_avg_for_masking.min()), stop=cov_avg_for_masking.max(), num=cutoff_num)
            cutoffs = cutoffs[1:];      cutoffs = cutoffs[:(len(cutoffs)-1)]
            scores = st.ms_cutoffs(v=cov_avg_for_masking, cutoffs=cutoffs)
            score_best = st.ms_mo_best(scores)
            score_best_cutoff = score_best['best']['cutoff']
            assert  score_best_cutoff >= 0

            cov_avg_for_masking__cutoff = max(cov_avg_for_masking__cutoff, score_best_cutoff)

        assert      cov_avg_for_masking__cutoff is not None
        voxel_mask_inds = np.flatnonzero(cov_avg_for_masking > cov_avg_for_masking__cutoff)

        print 'max covariance %f, score_best cutoff %f, cov_avg_for_masking__feature_num_cutoff %f, cov_avg_for_masking__cutoff %f, %d features selected for dimension reduction'%(cov_avg_for_masking.max(), score_best_cutoff, cov_avg_for_masking__feature_num_cutoff, cov_avg_for_masking__cutoff, len(voxel_mask_inds))

    # perform pca using only voxels / features indicated by vol_msk
    mat = pca__stack_dif__parallel(self=self, avg_key=avg_key, data=data, smoothing_gauss_sigma=diff_smoothing_gauss_sigma, n_chunk=n_chunk, voxel_mask_inds=voxel_mask_inds)
   
   
    if pca_op['mode'] == 0:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=dims)
        pca.fit(mat)
        red = pca.transform(mat)

    elif pca_op['mode'] == 1:
        from sklearn.decomposition import RandomizedPCA
        pca = RandomizedPCA(n_components=dims) 
        pca.fit(mat)
        red = pca.transform(mat)

    elif pca_op['mode'] == 2:
        import multiprocessing
        from sklearn.decomposition import SparsePCA
        pca = SparsePCA(n_components=dims, alpha=pca_op['alpha'], n_jobs=min(multiprocessing.cpu_count(), pca_op['job_num']), verbose=True)
        pca.fit(mat)
        red = pca.transform(mat)

    elif pca_op['mode'] == 3:
        mat_t = np.copy(mat)
        mat_t[np.isnan(mat)] = 0.0

        # weight missing value regions as 0
        empca_weight = np.isfinite(mat)

        if cov_avg_for_masking is not None:
            # weight every feature according to its corresponding average correlation
            cov_avg_for_masking__v = cov_avg_for_masking.flatten()
            for i, ind_t in enumerate(voxel_mask_inds):
                empca_weight[:,i] *= cov_avg_for_masking__v[ind_t]

        import tomominer.dimension_reduction.empca as drempca
        pca = drempca.empca( data=mat_t, weights=empca_weight, nvec=dims, niter=pca_op['n_iter'] )     # note: need to watch out the R2 values to see how much variation can be explained by the estimated model, if the value is small, need to increase dims

        red = pca.coeff


    print "PCA with covariange thresholding  : %2.6f sec" % (time.time() - start_time)

    return {'red':red, 'cov_avg':cov_avg, 'cov_avg_g':cov_avg_g, 'cov_avg_for_masking':cov_avg_for_masking, 'cov_avg_for_masking__cutoff':cov_avg_for_masking__cutoff, 'scores':scores, 'score_best':score_best, 'score_best_cutoff':score_best_cutoff, 'voxel_mask_inds':voxel_mask_inds}






