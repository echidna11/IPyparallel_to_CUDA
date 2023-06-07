
# calculate Spectial Signal to Noise Ratio

'''
~/ln/tomominer/tomominer/statistics/ssnr.py


for derivations, see google doc: parallel SSNR
https://docs.google.com/document/d/1vPVAPQStkjrzjiwRqHwfUrJEzrST4mr1sKGJdnfTrVU/edit


'''




import os
import sys
import time
import cPickle as pickle
import copy

import numpy as N
import numpy.fft as NF


from multiprocessing.pool import Pool as Pool


import tomominer.image.vol.util as GV
import tomominer.io.file as IV
import tomominer.geometry.rotate as GR
import tomominer.geometry.ang_loc as AAL


def ssnr_to_fsc(ssnr):
    fsc = ssnr / (2.0 + ssnr)
    return fsc


def var__local(self, data_json, labels=None, mask_cutoff=0.5, return_key=True, segmentation_tg_op=None):
    #print 'var__local()', 'segmentation_tg_op', segmentation_tg_op

    if labels is None:        labels = [0] * len(data_json)

    sum_v = {}
    prod_sum_v = {}
    mask_sum = {}

    for i, r in enumerate(data_json):
        if (self is not None) and (self.work_queue is not None) and self.work_queue.done_tasks_contains(self.task.task_id):            raise Exception('Duplicated task')         # terminate alignment process when the corresponding task is already completed by another worker

        v = IV.read_mrc_vol(r['subtomogram']);       v = GR.rotate_pad_mean(v, angle=N.array(r['angle'], dtype=N.float), loc_r=N.array(r['loc'], dtype=N.float) )
        m = IV.read_mrc_vol(r['mask']);       m = GR.rotate_mask(m, N.array(r['angle'], dtype=N.float) )

        if (segmentation_tg_op is not None) and ('template' in r) and ('segmentation' in r['template']):
            phi = IV.read_mrc(r['template']['segmentation'])['value']

            import tomominer.pursuit.multi.util as PMU
            v_s = PMU.template_guided_segmentation(v=v, m=(phi>0.5), op=segmentation_tg_op)

            if v_s is not None:
                v = v_s;                  del v_s
                
                v_t = N.zeros(v.shape)
                v_f = N.isfinite(v)
                v_t[v_f] = v[v_f]
                v_t[N.logical_not(v_f)] = v[v_f].mean()
                v = v_t;                    del v_f, v_t


        v = NF.fftshift( NF.fftn(v) )
        v[m < mask_cutoff] = 0.0        # in such case F_{is} <-- M_{is}F_{is} , see the derivation at the end of this file

        if labels[i] not in sum_v:
            sum_v[labels[i]] = v
        else:
            sum_v[labels[i]] += v

        if labels[i] not in prod_sum_v:
            prod_sum_v[labels[i]] = v * N.conj(v)
        else:
            prod_sum_v[labels[i]] += v * N.conj(v)

        if labels[i] not in mask_sum:
            mask_sum[labels[i]] = N.zeros(m.shape, dtype=N.int)

        mask_sum[labels[i]][m >= mask_cutoff] += 1


    re = {'sum':sum_v, 'prod_sum':prod_sum_v, 'mask_sum':mask_sum}

    if return_key:
        re_key = self.cache.save_tmp_data(re, fn_id=self.task.task_id)
        assert re_key is not None
        return {'key':re_key}
    else:
        return re


def var__global(self, clusters, n_chunk, segmentation_tg_op=None):
    #print 'var__global()', segmentation_tg_op


    data_json_copy = []
    labels_copy = []
    for l in clusters:
        for d in clusters[l]:
            data_json_copy.append(d)
            labels_copy.append(l)

    tasks = []
    while data_json_copy:
        data_json_copy_part = data_json_copy[:n_chunk]
        labels_copy_part = labels_copy[:n_chunk]

        tasks.append(self.runner.task(module='tomominer.statistics.ssnr', method='var__local', kwargs={'data_json':data_json_copy_part, 'labels':labels_copy_part, 'return_key':True, 'segmentation_tg_op':segmentation_tg_op}))

        data_json_copy = data_json_copy[n_chunk:]
        labels_copy = labels_copy[n_chunk:]


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

# get index within certain frequency band
def ssnr_rad_ind(rad, r, band_width_radius):
    return ( abs(rad - r) <= band_width_radius )

# get a volume that contains SSNR at different frequency bands
def ssnr_vol(ssnr, siz, band_width_radius=1.0):

    v = N.zeros(siz)

    rad = ssnr__get_rad(siz)
    for r in range(len(ssnr)):
        ind = ssnr_rad_ind(rad, r, band_width_radius=band_width_radius)
        v[ind] = ssnr[r]

    return v





def ssnr__given_stat(sum_v, prod_sum, mask_sum, rad=None, op=None):

    op = copy.deepcopy(op)

    if 'band_width_radius' not in op:
        op['band_width_radius'] = 1.0

    if 'mask_sum_threshold' not in op:
        op['mask_sum_threshold'] = 2.0
    else:
        op['mask_sum_threshold'] = N.max((op['mask_sum_threshold'], 2.0))         # should be at least 2

    if 'debug' not in op:
        op['debug'] = False

    if op['debug']:     print 'ssnr__given_stat()', op

    siz = N.array(sum_v.shape)

    subtomogram_num = mask_sum.max()

    avg = N.zeros(sum_v.shape, dtype=N.complex) + N.nan
    ind = mask_sum > 0
    avg[ind] = sum_v[ind]  /  mask_sum[ind];         avg_abs_sq = N.real(  avg * N.conj( avg ) )
    del ind

    var = N.zeros(sum_v.shape, dtype=N.complex) + N.nan
    ind = mask_sum >= op['mask_sum_threshold']
    var[ind] = ( prod_sum[ind] - mask_sum[ind] * ( avg[ind] * N.conj( avg[ind] ) ) ) / ( mask_sum[ind] - 1 );          var = N.real(var)        # mxu 150223: bug correction, should use (mask_sum - 1) instead of (mask_sum + 1)
    del ind

    if rad is None:     rad = ssnr__get_rad(siz)

    vol_rad = int( N.floor( N.min(siz) / 2.0 ) + 1)
    ssnr = N.zeros(vol_rad) + N.nan     # this is the SSNR of the AVERAGE image

    # the interpolation can also be performed using scipy.ndimage.interpolation.map_coordinates()
    for r in range(vol_rad):

        ind = ssnr_rad_ind(rad=rad, r=r, band_width_radius=op['band_width_radius'])         # in order to use it as an index or mask, must convert to a bool array, not integer array!!!!
        ind[mask_sum < op['mask_sum_threshold']] = False
        ind[N.logical_not(N.isfinite(avg))] = False
        ind[N.logical_not(N.isfinite(var))] = False

        if op['method'] == 0:
            if var[ind].sum() > 0:
                ssnr[r] = avg_abs_sq[ind].sum() / (var[ind] / mask_sum[ind]).sum()      # mxu 150227: this is modified version reverting to hg version of 72. Details see the google document referred in the head of this file
            else:
                ssnr[r] = 0.0

        elif op['method'] == 1:
            if var[ind].sum() > 0:
                ssnr[r] = (mask_sum[ind] * avg_abs_sq[ind]).sum() / var[ind].sum()
            else:
                ssnr[r] = 0.0

        elif op['method'] == 2:
            if var[ind].sum() > 0:
                ssnr[r] = subtomogram_num * avg_abs_sq[ind].sum() / var[ind].sum()
            else:
                ssnr[r] = 0.0


        elif op['method'] == 3:
            # a kind of total SSNR, such way of calculation is still NOT throughly analysed

            mask_total = mask_sum[ind].sum()
            avg_total = sum_v[ind].sum() / mask_total           ;           avg_total_abs_sq = N.real(avg_total * N.conj(avg_total))
            var_total = (prod_sum[ind].sum() - mask_total * avg_total * N.conj(avg_total)) / (mask_total - 1);      var_total = N.real(var_total)
            
            if var_total > 0:
                ssnr[r] = subtomogram_num * avg_total_abs_sq / var_total           # we multiply subtomogram_num because SSNR of average is subtomogram_num times of SSNR of individual subtomogram
            else:
                ssnr[r] = 0.0

        else:
            raise Exception('method')
            
        del ind

    assert N.all(N.isfinite(ssnr))

    # calculate corresponding Fourier Shell Correlation
    fsc = ssnr_to_fsc(ssnr)

    return {'ssnr':ssnr, 'fsc':fsc}





def ssnr_parallel(self, clusters, n_chunk=None, band_width_radius=1.0, op=None):

    c = var__global(self=self, clusters=clusters, n_chunk=n_chunk, segmentation_tg_op=(op['segmentation_tg'] if 'segmentation_tg_op' in op else None))

    ssnr = {}
    fsc = {}
    avg = {}
    var = {}

    for l in c['sum']:

        ssnr_re = ssnr__given_stat(sum_v=c['sum'][l], prod_sum=c['prod_sum'][l], mask_sum=c['mask_sum'][l], op=op['ssnr'])

        ssnr[l] = ssnr_re['ssnr']
        fsc[l] = ssnr_re['fsc']
    
    return {'ssnr':ssnr, 'fsc':fsc, 'avg':avg, 'var':var}




# this function is used for parallel calculation for sequential ssnr
def ssnr_sequential___batch_data_collect(self, data_json, op):

    if 'save_tmp_data' not in op:       op['save_tmp_data'] = False

    d = [None] * len(data_json)
    for i, r in enumerate(data_json):
        if self.work_queue.done_tasks_contains(self.task.task_id):            raise Exception('Duplicated task')         # terminate alignment process when the corresponding task is already completed by another worker

        d[i] = ssnr_sequential___individual_data_collect(self, r, op)

    if op['save_tmp_data']:
        re_key = self.cache.save_tmp_data(d={'data':d, 'op':op}, fn_id=self.task.task_id)
        assert re_key is not None       # if re_key is None, the some other work already saved the result of same task, we just terminate this function so that it will be ignored. IMPORTANT: do not return None, otherwise it may suppress the real successful copy!!

        # in such case, the result is successfully saved
        op['re_key'] = re_key
        return {'op':op}

    else:
        return {'data':d, 'op':op}



# this function is used to collect information of individual subtomograms, for sequential ssnr calculation purpose
def ssnr_sequential___individual_data_collect(self, r, op):

    if self is None:
        get_mrc_func = IV.get_mrc
    else:
        get_mrc_func = self.cache.get_mrc


    # ------------------------
    # collect statistics
    v = get_mrc_func(r['subtomogram'])
    if 'angle' in r:        v = GR.rotate_pad_mean(v, angle=N.array(r['angle'], dtype=N.float), loc_r=N.array(r['loc'], dtype=N.float) )

    if (op is not None) and ('segmentation_tg' in op) and ('template' in r) and ('segmentation' in r['template']):
        # segment v if needed
        
        phi = IV.read_mrc_vol(r['template']['segmentation'])
        phi_m = phi > 0.5             # why phi>0.5 instead of 0? Because we do not want to include boundary region?
        del phi

        ang_inv, loc_inv = AAL.reverse_transform_ang_loc(r['angle'], r['loc'])
        phi_mr = GR.rotate(phi_m, angle=ang_inv, loc_r=loc_inv, default_val=0)
        del phi_m
        del ang_inv, loc_inv

        import tomominer.pursuit.multi.util as PMU
        v_s = PMU.template_guided_segmentation(v=v, m=phi_mr, op=op['segmentation_tg'])         # why phi>0.5 instead of 0? Because we do not want to include boundary region?
        del phi_mr

        if v_s is not None:
            v_f = N.isfinite(v_s)
            if v_f.sum() > 0:
                v_s[N.logical_not(v_f)] = v_s[v_f].mean()
                v = v_s
            del v_s
       

    v = NF.fftshift( NF.fftn(v) )

    m = get_mrc_func(r['mask'])
    if 'angle' in r:    m = GR.rotate_mask(m, angle=N.array(r['angle'], dtype=N.float) )

    v[m < op['mask_cutoff']] = 0.0      # in such case F_{is} <-- M_{is}F_{is} , see the derivation at the end of this file

    return {'v':v, 'm':m}



# given a list of subtomograms and alignment information, for every n, calculate the ssnr for the first n subtomograms
def ssnr_sequential(self=None, data_json=None, op=None):

    assert data_json is not None

    if 'mask_cutoff' not in op:     op['mask_cutoff'] = 0.5
    if 'band_width_radius' not in op:   op['band_width_radius'] = 1.0


    sum_v = None
    prod_sum = None
    mask_sum = None

    rad = None

    ssnr_s = [None] * len(data_json)
    fsc_s = [None] * len(data_json)

    for i, r in enumerate(data_json):
        tmp = ssnr_sequential___individual_data_collect(self=self, r=r, op=op)
        v = tmp['v']
        m = tmp['m']
        del tmp


        if sum_v is None:
            sum_v = v
        else:
            sum_v += v

        if prod_sum is None:
            prod_sum = v * N.conj(v)
        else:
            prod_sum += v * N.conj(v)

        if mask_sum is None:
            mask_sum = N.zeros(m.shape, dtype=N.int)

        mask_sum[m >= op['mask_cutoff']] += 1


        # -----------------------
        # calculate SSNR and FSC

        if rad is None:          rad = ssnr__get_rad(N.array(sum_v.shape))
        
        ssnr_re = ssnr__given_stat(rad=rad, sum_v=sum_v, prod_sum=prod_sum, mask_sum=mask_sum, method=op['ssnr'])

        ssnr_s[i] = ssnr_re['ssnr']
        fsc_s[i] = ssnr_re['fsc']

        print '\r', '%d   %0.3f   '%(i, float(i) / len(data_json)), '     ',  fsc_s[i].sum(), '       ',
        sys.stdout.flush()

    return {'ssnr':ssnr_s, 'fsc':fsc_s, 'op':op}


# use multiprocessing to parallelly calculate sequential SSNR
def ssnr_sequential_multiprocessing(self=None, data_json_dict=None, band_width_radius=1.0, mask_cutoff=0.5, ssnr_sequential_op=None):

    # parallel calculationg
    pool = Pool()
    pool_results = []

    for c in data_json_dict:
        ssnr_sequential_op_t = copy.deepcopy(ssnr_sequential_op)
        ssnr_sequential_op_t['set_id'] = c
        pool_results.append( pool.apply_async(func=ssnr_sequential, kwds={'data_json':data_json_dict[c], 'op':ssnr_sequential_op_t}) )
        
    ssnr_s = {}
    for r in pool_results:
        r = r.get(999999)
        ssnr_s[r['op']['set_id']] = r

    pool.close()

    return ssnr_s



# collect data in parallel that is used for ssnr sequential calculation
# then load the data in correct order and calculate sequential SSNR
# parameters:       data_json_dict: sets of subtomograms
def ssnr_sequential_parallel(self, data_json_dict, op):

    if 'verbose' not in op:     op['verbose'] = False
    if op['verbose']:   print 'ssnr_sequential_parallel()'

    if 'mask_cutoff' not in op:     op['mask_cutoff'] = 0.5
    if 'band_width_radius' not in op:   op['band_width_radius'] = 1.0


    tasks = []
    for c in data_json_dict:
        dj = copy.deepcopy(data_json_dict[c])
        inds = list(range(len(dj)))

        sequential_id = 0

        while dj:
            dj_part = dj[:op['n_chunk']]
            inds_t = inds[:op['n_chunk']]

            opt = copy.deepcopy(op)
            opt['cluster'] = c
            opt['sequential_id'] = sequential_id
            opt['inds'] = inds_t
            opt['save_tmp_data'] = True

            tasks.append(self.runner.task(module='tomominer.statistics.ssnr', method='ssnr_sequential___batch_data_collect', kwargs={'data_json':dj_part, 'op':opt}))

            dj = dj[op['n_chunk']:]
            inds = inds[op['n_chunk']:]

            sequential_id += 1


    # collect and arrange the returned information
    data_op = {}
    for re in self.runner.run__except(tasks):
        re = re.result

        if re is None: continue
        re = re['op']
        if re['cluster'] not in data_op:      data_op[re['cluster']] = {}
        data_op[ re['cluster'] ][ re['sequential_id'] ] = re


    # for all clusters / sets, load result in sequence, and calculate sequential SSNR
    ssnr_s = {}
    for c in data_op:
        sum_v = None
        prod_sum = None
        mask_sum = None

        rad = None

        ssnr_s[c] = {}
        ssnr_s[c]['ssnr'] = [None] * len(data_json_dict[c])
        ssnr_s[c]['fsc'] = [None] * len(data_json_dict[c])

        current_ind = -1

        for data_op_i in range(len(data_op[c])):
            opt = data_op[c][data_op_i]

            with open(opt['re_key'], 'rb') as f:     data_t = pickle.load(f)
            os.remove(opt['re_key'])

            data_t = data_t['data']

            assert  len(data_t) == len(opt['inds'])

            for ind_i, ind in enumerate(opt['inds']):

                assert  ind == (current_ind + 1)     # make sure that the data is in correct sequence
                current_ind += 1

                v = data_t[ind_i]['v']
                m = data_t[ind_i]['m']

                if sum_v is None:
                    sum_v = v
                else:
                    sum_v += v

                if prod_sum is None:
                    prod_sum = v * N.conj(v)
                else:
                    prod_sum += v * N.conj(v)

                if mask_sum is None:
                    mask_sum = N.zeros(m.shape, dtype=N.int)

                mask_sum[m >= op['mask_cutoff']] += 1


                # -----------------------
                # calculate SSNR and FSC

                if rad is None:          rad = ssnr__get_rad(N.array(sum_v.shape))
                
                ssnr_re = ssnr__given_stat(rad=rad, sum_v=sum_v, prod_sum=prod_sum, mask_sum=mask_sum, op=op['ssnr'])

                assert  ssnr_s[c]['ssnr'][ind] is None              ;               ssnr_s[c]['ssnr'][ind] = ssnr_re['ssnr']
                assert  ssnr_s[c]['fsc'][ind] is None               ;               ssnr_s[c]['fsc'][ind] = ssnr_re['fsc']
                del ssnr_re

                print '\r', '%d   %0.3f   '%(ind, float(ind) / len(data_json_dict[c])), '     ',  ssnr_s[c]['fsc'][ind].sum(), '       ',
                sys.stdout.flush()
        
    for c in ssnr_s:
        for t in ssnr_s[c]['ssnr']:         assert t is not None
        for t in ssnr_s[c]['fsc']:          assert t is not None

    return ssnr_s




'''

# test code

# first run DataProcessingSSNR.sequential_subtomogram_generate__commands() to generate test data
import json

with open('/home/rcf-47/mxu/tmp-staging/ssnr-test/info.json') as f:     dj = json.load(f)



# just get a generic object that can contain a cache
class Object:
    pass
self = Object()
 
from tomominer.io.cache import Cache
self.cache = Cache()

ssnr_sequential(self=self, data_json=dj)

'''




