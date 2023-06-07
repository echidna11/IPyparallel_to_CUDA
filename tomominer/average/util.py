
# utility functions to support averaging

import os
import time
import uuid
import copy

import numpy as np
from numpy.fft import fftn, ifftn, fftshift, ifftshift

import shutil

import tomominer.core as tomo

import tomominer.io.file as iv






def load_avg(avg_key, vol_mask_avg_threshold=0.0):
    vol_avg_key, vol_mask_avg_key = avg_key

    vol_avg = iv.get_mrc(vol_avg_key)       # we do not cache average
    vol_avg_fft = fftshift(fftn(vol_avg))

    vol_mask_avg = iv.get_mrc(vol_mask_avg_key)
    if vol_mask_avg_threshold is not None:
        vol_mask_avg = np.array((vol_mask_avg > vol_mask_avg_threshold), dtype=np.float)    # just discritize to denote Fourier space regions that with or without information 

    return {'vol_avg':vol_avg, 'vol_avg_fft':vol_avg_fft, 'vol_mask_avg':vol_mask_avg}



# given a set of aligned subtomograms arbitrarily separate to two equal size sets, then calculate averages of each set
def average_halves(op):
    data = op['data']

    reminders = [0, 1]

    data_h = [None] * len(reminders)
    for reminder in reminders:
        data_h[reminder] = [data[i] for i in range(len(data)) if ((i%2) == reminder)]

    avgs = [None] * len(reminders)
    for reminder in reminders:
        avgs[reminder] = average_rotate_parallel({'runner':op['runner'], 'data':data_h[reminder], 'vol_shape':op['vol_shape'], 'tmp_dir':op['tmp_dir']})


    return {'data_halves':data_h, 'averages':avgs}
    

# rotate subtomograms according to given rotation, then average them in parallel
def average_rotate_parallel(op):

    data = op['data']

    n_chunk = 50

    vmal = [_ for _ in data]
    n_vol = len(vmal)
    tasks = []
    cnt = 0
    while vmal:
        #tasks.append(op['runner'].task('vol_avg_local', vol_shape, vmal[:n_chunk], v_out_key, m_out_key))
        tasks.append(op['runner'].task('vol_avg_fft_parallel_local', op['vol_shape'], vmal[:n_chunk]))

        vmal = vmal[n_chunk:]
        cnt += 1

    [_ for _ in op['runner'].run__except(tasks)]
    assert False        # need to correct




    vol_avg_out_key = os.path.join(op['tmp_dir'], 'vol-avg-global--' + str(uuid.uuid1()) + '.mrc')
    mask_avg_out_key = os.path.join(op['tmp_dir'], 'maks-avg-global--' + str(uuid.uuid1()) + '.mrc')

    #tasks = [op['runner'].task('vol_avg_global', vol_shape, n_vol, local_vm, vol_avg_out_key)]
    vol_avg_fft_parallel_global(vol_shape=op['vol_shape'], n_vol=n_vol, vm_in=local_vm, vol_avg_key=vol_avg_out_key, mask_avg_key=mask_avg_out_key)


    vol_avg = iv.get_mrc(vol_avg_out_key, tmp_dir)
    mask_avg = iv.get_mrc(mask_avg_out_key, tmp_dir)

    # clean up temporary data
    for loval_vm_t in local_vm:
        for loval_vm_tt in loval_vm_t:
            try: 
                os.remove(loval_vm_tt)
            except:
                pass    

    try:
        os.remove(vol_avg_out_key)
        os.remove(mask_avg_out_key)
    except:
        pass

    return {'vol':vol_avg, 'mask':mask_avg}



# mxu: calculating global averages written by Zach, moved from classify.py
def global_average_parallel(data, vol_shape, n_chunk, pass_dir, runner, use_fft_avg, centerize_loc=False):

    vol_avg_out_key = os.path.join(pass_dir, 'vol_avg.mrc')
    mask_avg_out_key = os.path.join(pass_dir, 'mask_avg.mrc')

    if os.path.exists(vol_avg_out_key) and os.path.exists(mask_avg_out_key):    return {'vol_avg_out_key':vol_avg_out_key, 'mask_avg_out_key':mask_avg_out_key}     # just save a little bit of computation

    start_time = time.time()


    if centerize_loc:
        data_t = copy.deepcopy(data)

        loc = np.zeros( (len(data_t), 3) )
        for i, (v,m,a,l) in enumerate(data_t):      loc[i,:] = l
        loc -= np.tile( loc.mean(axis=0), (loc.shape[0], 1) )       # substract mean so that mean of loc is equal to zero vector. The purpose of centerize_loc is to reduce the chance of clipping brought by large displacements
        assert np.all(np.abs(loc.mean(axis=0)) <= 1e-10)
        
        for i, (v,m,a,l) in enumerate(data_t):      data_t[i] = (v,m,a,loc[i])

    else:
        data_t = data;


    
    if use_fft_avg:
        vol_avg_fft_parallel_global(vol_shape=vol_shape, data=data_t, runner=runner, n_chunk=n_chunk, vol_avg_key=vol_avg_out_key, mask_avg_key=mask_avg_out_key)
    else:
        vol_avg_parallel_global(vol_shape=vol_shape, data=data_t, runner=runner, n_chunk=n_chunk, vol_avg_key=vol_avg_out_key, mask_avg_key=mask_avg_out_key)


    print "Global average volume computed: %2.6f sec" % (time.time() - start_time)

    return {'vol_avg_out_key':vol_avg_out_key, 'mask_avg_out_key':mask_avg_out_key}



def vol_avg_fft_parallel_local(vol_shape, vmal_in, tmp_dir=None):
    """
    Local work for computing average of all volumes.  This calculates a sum
    over the given subset of volumes.

    Several of these are done and combined with vol_avg_global.  This is
    the map() stage of the global averaging step.
    """

    # TODO: decide if we will ever use non-standard weights in this step.
    weight = 1.0

    # temporary collection of local volume, and mask.
    vol_sum  = np.zeros(vol_shape, dtype=np.complex, order='F')
    mask_sum = np.zeros(vol_shape, order='F')

    # iterate over all volumes/masks and incorporate data into averages.
    # volume_key, mask_key, angle offset, location offset.
    for vk, mk, ang, loc in vmal_in:
        # load vol/mask.
        vol  = iv.get_mrc_cache_fs(vk, tmp_dir)
        mask = iv.get_mrc_cache_fs(mk, tmp_dir)

        # rotate vol and mask according to angle/loc.
        vol  = tomo.rotate_vol_pad_mean_py(vol, ang, loc)
        mask = tomo.rotate_mask_py(mask, ang)

        # TODO: apparently this casts to real?
        vol_sum  += weight * (fftshift(fftn(vol)) * mask)
        mask_sum += weight * mask


    return {'vol_sum':vol_sum, 'mask_sum':mask_sum}


def vol_avg_fft_parallel_global(vol_shape, data, runner, n_chunk, vol_avg_key, mask_avg_key):
    """
    Reduce step of global average.
    """
    vmal = [_ for _ in data]
    n_vol = len(vmal)
    tasks = []
    cnt = 0
    while vmal:
        tasks.append(runner.task('vol_avg_fft_parallel_local', vol_shape=vol_shape, vmal_in=vmal[:n_chunk]))

        vmal = vmal[n_chunk:]
        cnt += 1


    # TODO: what is the correct value?
    mask_threshold = 1.0

    fft_sum  = np.zeros(vol_shape, dtype=np.complex, order='F')
    mask_sum = np.zeros(vol_shape, order='F')

    for res in runner.run__except(tasks):
        fft_sum += res.result['vol_sum']
        mask_sum += res.result['mask_sum']


    # mxu: corrected RuntimeWarning: invalid value encountered in divide
    flag_t = (mask_sum >= mask_threshold)
    tem_fft = np.zeros(vol_shape, dtype=np.complex, order='F')
    tem_fft[flag_t] = fft_sum[flag_t] / mask_sum[flag_t]
    #tem_fft[mask_sum < mask_threshold] = 0

    vol_avg = np.asfortranarray(np.real(ifftn(ifftshift(tem_fft))))
    mask_avg = mask_sum / n_vol

    iv.put_mrc(vol_avg, vol_avg_key)
    iv.put_mrc(mask_avg, mask_avg_key)

#        # stash the avg volume into the cache since we will be using it next.
#        if add_to_cache:
#            self.mrc_cache[vol_avg_key] = vol_avg

    return True


def template_avg_fft_parallel_global(vol_shape, runner, n_vol, vol_avg_key, mask_avg_key, add_to_cache, template_mask_ratio_cutoff):
    """
    Reduce step of global average when processing a template.
    """

    # TODO: what is the correct value?
    mask_threshold = 1.0

    fft_sum  = np.zeros(vol_shape, dtype=np.complex, order='F')
    mask_sum = np.zeros(vol_shape)

    for res in runner.run__except(tasks):
        fft_sum += res.result['vol_sum']
        mask_sum += res.result['mask_sum']

    # mxu: solved RuntimeWarning: invalid value encountered in divide
    tem_fft = np.zeros(vol_shape, dtype=np.complex, order='F')
    mask_threshold_mask = (mask_sum >= mask_threshold)
    tem_fft[mask_threshold_mask] = fft_sum[mask_threshold_mask] / mask_sum[mask_threshold_mask]

    vol_avg = np.asfortranarray(np.real(ifftn(ifftshift(tem_fft))))

    if False:
        mask_tem = np.zeros(mask_sum.shape, order='F')
        mask_tem[mask_sum >= template_mask_ratio_cutoff] = 1
    else:
        mask_tem = mask_sum / n_vol;
        mask_tem = np.asfortranarray(mask_tem)    


    # stash the avg volume into the cache since we will be using it next.
    if add_to_cache:
        self.mrc_cache[vol_avg_key] = vol_avg

    if mask_avg_key and add_to_cache:
        self.mrc_cache[mask_avg_key] = mask_tem

    iv.put_mrc(vol_avg, vol_avg_key)

    if mask_avg_key:
        iv.put_mrc(mask_tem, mask_avg_key)

    return True

def vol_avg_parallel_local(vol_shape, vmal_in, tmp_dir=None):
    """
    Local work for computing average of all volumes.  This calculates a sum
    over the given subset of volumes.

    Several of these are done and combined with vol_avg_global.  This is
    the map() stage of the global averaging step.
    """

    # TODO: decide if we will ever use non-standard weights in this step.
    weight = 1.0

    # temporary collection of local volume, and mask.
    vol_sum  = np.zeros(vol_shape, dtype=np.float64, order='F')
    mask_sum = np.zeros(vol_shape, order='F')

    # iterate over all volumes/masks and incorporate data into averages.
    # volume_key, mask_key, angle offset, location offset.
    for vk, mk, ang, loc in vmal_in:
        # load vol/mask.
        vol  = iv.get_mrc_cache_fs(vk, tmp_dir)
        mask = iv.get_mrc_cache_fs(mk, tmp_dir)

        # rotate vol and mask according to angle/loc.
        vol  = tomo.rotate_vol_pad_mean_py(vol, ang, loc)
        mask = tomo.rotate_mask_py(mask, ang)

        # TODO: apparently this casts to real?
        vol_sum  += weight * vol
        mask_sum += weight * mask


    return {'vol_sum':vol_sum, 'mask_sum':mask_sum}


def vol_avg_parallel_global(vol_shape, data, runner, n_chunk, vol_avg_key, mask_avg_key):
    """
    Reduce step of global average.  The real space average

    :param vol_shape The shape of the volumes
    :type vol_shape array of 3 elements x,y,z coords of the volumes.
    :param n_vol The total number of volumes that are being average across the vol_leys.
    """
    # First calculate global average volume.
    vmal = [_ for _ in data]
    n_vol = len(vmal)
    tasks = []
    cnt = 0
    while vmal:

        tasks.append(runner.task('vol_avg_parallel_local', vol_shape=vol_shape, vmal_in=vmal[:n_chunk]))

        vmal = vmal[n_chunk:]
        cnt += 1


    vol_sum  = np.zeros(vol_shape, order='F')

    for res in runner.run__except(tasks):
        vol_sum += res.result['vol_sum']

    vol_avg = vol_sum / n_vol

    iv.put_mrc(vol_avg, vol_avg_key)
    iv.put_mrc(np.ones(vol_shape, order='F', dtype=np.float64), mask_avg_key)
    return True


def template_avg_parallel_global(vol_shape, runner, n_vol, vol_avg_key, mask_avg_key, add_to_cache, template_mask_ratio_cutoff):
    """
    Reduce step of global average when processing a template.
    """

    # TODO: what is the correct value?
    mask_threshold = 1.0

    vol_sum  = np.zeros(vol_shape, order='F')

    for res in runner.run__except(tasks):
        vol_sum += res.result['vol_sum']

    vol_avg = vol_sum / n_vol

    if add_to_cache:
        self.mrc_cache[vol_avg_key] = vol_avg

    iv.put_mrc(vol_avg, vol_avg_key)

    if mask_avg_key:
        self.put_mrc(np.ones(vol_shape, order='F', dtype=np.float64), mask_avg_key)

    return True


#    def vol_avg_global(self, vol_shape, n_vol, vol_keys, vol_avg_key):
#        """
#        Reduce step of global average.  The real space average
#
#        :param vol_shape The shape of the volumes
#        :type vol_shape array of 3 elements x,y,z coords of the volumes.
#        :param n_vol The total number of volumes that are being average across the vol_leys.
#        """
#
#        vol_sum  = np.zeros(vol_shape, order='F')
#
#        for vk,mk in vol_keys:
#            vol = self.np_load(vk)
#            vol_sum += vol
#
#        vol_avg = vol_sum / n_vol
#
#        self.put_mrc(vol_avg, vol_avg_key)
#        return True


