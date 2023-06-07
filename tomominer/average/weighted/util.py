
'''
~/ln/tomominer/tomominer/average/weighted/util.py


utility functions for averaging

'''


import os
import sys
import copy
import time
import traceback
import warnings
import cPickle as pickle


import numpy as N
import numpy.fft as NF


import tomominer.io.file as IV
import tomominer.geometry.rotate as GR
import tomominer.align.util as AU




def cluster_averaging(self, clusters, op={}):

    print 'average.weighted.util.cluster_averaging()', op

    cav =  cluster_averaging_vols(self, clusters=clusters, op=op)

    if not os.path.isdir(op['out_dir']):      os.makedirs(op['out_dir'])
    with open(os.path.join(op['out_dir'], 'cluster.pickle'), 'wb') as f:      pickle.dump(cav, f, protocol=-1)



    # save averages
    cluster_avg_dict = cav['cluster_avg_dict']
    template_keys = {}
    for c in cluster_avg_dict:
        vol_avg_out_key = os.path.join(op['out_dir'], 'clus_vol_avg_%03d.mrc'%(c))
        IV.put_mrc(N.array(cluster_avg_dict[c]['vol'], order='F'), vol_avg_out_key)

        mask_avg_out_key = os.path.join(op['out_dir'], 'clus_mask_avg_%03d.mrc'%(c))
        IV.put_mrc(N.array(cluster_avg_dict[c]['mask'], order='F'), mask_avg_out_key)

        template_keys[c] = {'cluster':c, 'subtomogram':vol_avg_out_key, 'mask':mask_avg_out_key}     # use key 'subtomogram' for the consistency with json record in data_config, for calling align_keys()
        if 'pass_i' in op:      template_keys[c]['pass_i'] = op['pass_i']


    return {'template_keys':template_keys}



def cluster_averaging_vols(self, clusters, op={}):

    start_time = time.time()

    if op['centerize_loc']:
        clusters_cnt = copy.deepcopy(clusters)

        for c in clusters_cnt:
            loc = N.zeros( (len(clusters_cnt[c]), 3) )
            for i, rec in enumerate(clusters_cnt[c]):
                if 'loc' in rec:
                    assert 'align' not in rec
                    loc[i,:] = N.array(rec['loc'], dtype=N.float)
                else:
                    assert 'loc' not in rec
                    loc[i,:] = N.array(rec['align'][0]['loc'], dtype=N.float)

            loc_mean= loc.mean(axis=0)       # substract mean so that mean of loc is equal to zero vector. The purpose of centerize_loc is to reduce the chance of clipping brought by large displacements
            #assert N.all(N.abs(loc_mean.mean(axis=0)) <= 1e-10)
            
            for i, rec in enumerate(clusters_cnt[c]):
                if 'align' in rec:
                    assert 'loc' not in rec
                    for a in rec['align']:
                        a['loc'] -= loc_mean
                else:
                    assert 'loc' in rec
                    rec['loc'] -= loc_mean

        clusters = clusters_cnt


    # parallel averaging
    tasks = []
    for c in clusters:
        while clusters[c]:
            part        = clusters[c][:op['n_chunk']]

            op_t = copy.deepcopy(op)
            op_t['cluster'] = c
            tasks.append(self.runner.task(module='tomominer.average.weighted.util', method='vol_avg__local', kwargs={'data_json':part, 'op':op_t, 'return_key':True}))

            clusters[c] = clusters[c][op['n_chunk']:]
    

    # calculate averages for each cluster
    #cluster_sizes = dict((c,len(clusters[c])) for c in clusters)
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




        ind = cluster_mask_sums[c] >= op['mask_count_threshold']
        if ind.sum() == 0:      continue            # we ignore small clusters, because their resulting average is purly zeros


        cluster_sums_fft = NF.fftshift( NF.fftn(cluster_sums[c]) )
        
        cluster_avg = N.zeros(cluster_sums_fft.shape, dtype = N.complex)

        cluster_avg[ind] = cluster_sums_fft[ind] / cluster_mask_sums[c][ind]
        
        cluster_avg = N.real( NF.ifftn(NF.ifftshift(cluster_avg)) )

        cluster_mask_avg = cluster_mask_sums[c] / cluster_sizes[c]


        cluster_avg_dict[c] = {'vol':cluster_avg, 'mask':cluster_mask_avg}

    
    print "Averaging cluster templates %2.6f sec" % (time.time() - start_time)

    return {'cluster_avg_dict':cluster_avg_dict, 'cluster_sums':cluster_sums, 'cluster_mask_sums':cluster_mask_sums, 'cluster_sizes':cluster_sizes}





def vol_avg__local(self, data_json, op=None, return_key=True):

    vol_sum = None
    mask_sum = None

    # temporary collection of local volume, and mask.
    for rec in data_json:

        if self.work_queue.done_tasks_contains(self.task.task_id):            raise Exception('Duplicated task')         # terminate alignment process when the corresponding task is already completed another worker

        v  = self.cache.get_mrc(rec['subtomogram'])
        if not N.all(N.isfinite(v)):        raise Exception('error loading', rec['subtomogram'])

        vm = self.cache.get_mrc(rec['mask'])
        if not N.all(N.isfinite(vm)):        raise Exception('error loading', rec['mask'])


        if 'align' in rec:
            assert 'angle' not in rec
            assert 'loc' not in rec

            score_sum = N.sum([_['score'] for _ in rec['align']])
            assert N.isfinite(score_sum)


            for a in rec['align']:
                v_r = GR.rotate_pad_mean(v, angle=a['angle'], loc_r=a['loc'])       ;           assert N.all(N.isfinite(v_r))
                vm_r = GR.rotate_mask(vm, angle=a['angle'])                    ;           assert N.all(N.isfinite(vm_r))

                v_r -= v_r.mean()               # make sure that the subtomogram has zero mean!! this makes weighting more reasonable!!!

                weight = a['score'] / score_sum

                if vol_sum is None:     vol_sum = N.zeros(v_r.shape, dtype=N.float64, order='F')
                vol_sum += v_r * weight

                if mask_sum is None:        mask_sum = N.zeros(vm_r.shape, dtype=N.float64, order='F')
                mask_sum += vm_r * weight
        else:
            assert 'align' not in rec

            v_r = GR.rotate_pad_mean(v, angle=rec['angle'], loc_r=rec['loc'])       ;           assert N.all(N.isfinite(v_r))
            vm_r = GR.rotate_mask(vm, angle=rec['angle'])                    ;           assert N.all(N.isfinite(vm_r))

            if vol_sum is None:     vol_sum = N.zeros(v_r.shape, dtype=N.float64, order='F')
            vol_sum += v_r

            if mask_sum is None:        mask_sum = N.zeros(vm_r.shape, dtype=N.float64, order='F')
            mask_sum += vm_r


    re = {'vol_sum':vol_sum, 'mask_sum':mask_sum, 'vol_count':len(data_json), 'op':op}

    if return_key:
        re_key = self.cache.save_tmp_data(re, fn_id=self.task.task_id)
        assert re_key is not None
        return {'key':re_key}

    else:

        return re




'''
adapted from tomominer.pursuit.multi.util.align_to_templates()
'''

def align_to_templates(self, rec=None, segmentation_op=None, tem_keys=None, template_wedge_cutoff=0.1, align_op=None, multiprocessing=False):
    #print 'align_to_templates()', rec;      sys.stdout.flush()

    v =  self.cache.get_mrc(rec['subtomogram'])
    vm =  self.cache.get_mrc(rec['mask'])

    align_re = {}
    for c in tem_keys:

        if self.work_queue.done_tasks_contains(self.task.task_id):            raise Exception('Duplicated task')         # terminate alignment process when the corresponding task is already completed by another worker

        #print tem_keys[c];      sys.stdout.flush()

        align_re[c] = align_to_templates__pair_align(t_key=tem_keys[c], v=v, vm=vm, align_op=align_op)

        if align_re[c]['err'] is not None:
            if self.logger is not None:      self.logger.warning("alignment failed: rec %s, template %s, error %s ", repr(rec), repr(tem_keys[c]), repr(align_re[c]['err']) )


    return {'vol_key':rec, 'align':align_re}


'''
adapted from tomominer.pursuit.multi.util.align_to_templates__pair_align()
'''
def align_to_templates__pair_align(t_key, v, vm, align_op):

    t = IV.get_mrc(t_key['subtomogram'])
    tm = IV.get_mrc(t_key['mask'])
    at_re = align_vols_with_wedge(v1=t, m1=tm, v2=v, m2=vm, op=align_op)
    
    return at_re



'''
adapted from tomominer.pursuit.multi.util.align_vols_with_wedge()
'''
def align_vols_with_wedge(v1, m1, v2, m2, op=None, logger=None):

    err = None
    try:
        al = AU.align_vols__multiple_rotations(v1=v1, m1=m1, v2=v2, m2=m2, L=op['L'])
    except Exception as e:      # mxu: just bypass alignment errors
        err = traceback.format_exc()


    if err is not None:        al = [{'score': float('nan'), 'loc': N.zeros(3), 'angle': N.random.random(3) * (N.pi * 2)}]      # mxu: randomly assign an angle

    return {'align':al, 'err':err}


