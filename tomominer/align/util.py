
import sys, copy
import numpy as N
import traceback

import tomominer.core as tomo

import tomominer.model.util as MU
import tomominer.io.file as iv

def align_vols(v1, m1, v2, m2, L):

    fail = False

    try:
        al = align_vols__multiple_rotations(v1=v1, m1=m1, v2=v2, m2=m2, L=L)

        # extract the score/displacement/angle from the first entry ret[0]
        score = al[0]['score']
        loc = al[0]['loc']
        angle = al[0]['angle']

    except Exception, err:
        print >>sys.stderr, traceback.format_exc()

        score = N.nan
        loc = N.zeros(3) + N.nan
        angle = N.zeros(3) + N.nan

        fail = True

    if not N.isfinite(score):  fail = True
    if len(loc) != 3:    fail = True
    if len(angle) != 3:    fail = True

    if not fail:
        return {'score':score, 'loc':loc, 'angle':angle}
    else:
        return {'score':float('nan'), 'loc':N.zeros(3), 'angle':N.random.random(3) * (N.pi * 2)}           # mxu: randomly assign an angle

def align_vols__multiple_rotations(v1, m1, v2, m2, L):

    if m1 is None:      m1 = MU.sphere_mask(v1.shape)
    if m2 is None:      m2 = MU.sphere_mask(v2.shape)

    assert(v1.shape == m1.shape)
    assert(v2.shape == m2.shape)
    if v1.shape != v2.shape:
        print v1.shape
        print v2.shape
        assert(v1.shape == v2.shape)

    cs = tomo.combined_search(v1.astype(N.float64), m1.astype(N.float64), v2.astype(N.float64), m2.astype(N.float64), L)

    al = [None] * len(cs)
    for i in range(len(cs)):        al[i] = {'score':cs[i][0], 'loc':N.copy(cs[i][1]), 'angle':N.copy(cs[i][2])}

    al = sorted(al, key=lambda _ : _['score'], reverse=True)            # make sure the alignment is in decreasing order

    return al



def align_keys(v1_key, m1_key, v2_key, m2_key, L, v2r_key=None, m2r_key=None, tmp_dir=None, logger=None):
    ea = N.array([0.0,0.0,0.0])
    dx = N.array([0.0,0.0,0.0])

    v1 = iv.get_mrc_cache_fs(v1_key, tmp_dir)
    m1 = iv.get_mrc_cache_fs(m1_key, tmp_dir)
    v2 = iv.get_mrc_cache_fs(v2_key, tmp_dir)
    m2 = iv.get_mrc_cache_fs(m2_key, tmp_dir)
    
    failure = False
    try:
        a_re = align_vols(v1, m1, v2, m2, L)
        score = a_re['score']
        loc = a_re['loc']
        ang = a_re['angle']

    except:      # mxu: just bypass alignment errors
        failure = True

    if failure:
        if logger is not None:      logger.warning("combined_search failed: v1_key %s, m1_key %s, v2_key %s, m2_key %s", v1_key, m1_key, v2_key, m2_key)

        score = float('nan')
        loc = N.zeros(3)
        ang = N.random.random(3) * (N.pi * 2)      # mxu: randomly assign an angle

    if v2r_key is not None:
        v2r = tomo.rotate_vol_pad_mean(v2, ang, loc)
        iv.put_mrc(v2r, v2r_key)

    if m2r_key is not None:
        m2r = tomo.rotate_mask(m2, ang)
        iv.put_mrc(m2r, m2r_key)

    
    return {'score':score, 'loc':loc, 'ang':ang, 'v1_key':v1_key, 'm1_key':m1_key, 'v2_key':v2_key, 'm2_key':m2_key, 'v2r_key':v2r_key, 'm2r_key':m2r_key}


def combined_search(v1_key, m1_key, v2_key, m2_key, L, v_out_key, m_out_key, tmp_dir=None, logger=None):

    ea = N.array([0.0,0.0,0.0])
    dx = N.array([0.0,0.0,0.0])

    v1 = iv.get_mrc_cache_fs(v1_key, tmp_dir)
    m1 = iv.get_mrc_cache_fs(m1_key, tmp_dir)
    v2 = iv.get_mrc_cache_fs(v2_key, tmp_dir)
    m2 = iv.get_mrc_cache_fs(m2_key, tmp_dir)

    try:
        ret = tomo.combined_search(v1, m1, v2, m2, L)

        # extract the score/displacement/angle from the first entry ret[0]
        score = ret[0][0]
        dx = N.copy(ret[0][1:4])
        ea = N.copy(ret[0][4:])

    except:      # mxu: just bypass alignment errors
        if logger is not None:      logger.warning("combined_search failed: v1_key %s, m1_key %s, v2_key %s, m2_key %s, v_out_key %s, m_out_key %s", v1_key, m1_key, v2_key, m2_key, v_out_key, m_out_key)

        score = -1
        dx = N.zeros(3)
        ea = N.random.random(3) * (N.pi * 2)      # mxu: randomly assign an angle
            
    v2 = tomo.rotate_vol_pad_mean(v2, ea, dx)
    m2 = tomo.rotate_mask(m2, ea)

    iv.put_mrc(v2, v_out_key)
    iv.put_mrc(m2, m_out_key)

    return True


def align_to_template(v, vm, t, tm, align_op=None, logger=None):

    e = None
    failure = False

    try:
        if align_op['mode'] == 0:
            #if logger is not None:      logger.debug("combined_search: v1_key %s, \n m1_key %s,\n v2_key %s,\n m2_key %s", v1_key, m1_key, v2_key, m2_key)

            re = align_vols(v1=t, m1=tm, v2=v, m2=vm, L=align_op['L'])
            ang = re['angle']
            loc = re['loc']
            score = re['score']

        elif align_op['mode'] == 1:
            raise
            #re = segment_and_align_one_vol_against_one_template(t=t, tg, t_seg, tm=tm, v=v, vg, vm=vm, L=align_op['L'])
            ang = re['angle']
            loc = re['loc']

        elif align_op['mode'] == 2:
            # align with interpolation
            raise RuntimeError('to be implemented')

        else:
            raise AttributeError( 'align_op[\'mode\'] : %d'%(align_op['mode']) )
  
    except Exception as e:      # mxu: just bypass alignment errors
        failure = True

    if failure:
        score = float('nan')
        loc = N.zeros(3)
        ang = N.random.random(3) * (N.pi * 2)      # mxu: randomly assign an angle


    return {'angle':ang, 'loc':loc, 'score':score, 'err':e}



# it seems tem_cutoff need to be low (say 0.01) to reduce the effect of missing wedge
def align_to_templates(vmal=None, vm_tem=None, align_op=None, tem_cutoff=None, transform=False, tmp_dir=None, logger=None):

    v_key, vm_key, ang, loc = vmal

    v = iv.get_mrc_cache_fs(v_key, tmp_dir)
    vm = iv.get_mrc_cache_fs(vm_key, tmp_dir)

    best_score = -2     # mxu: value changed to -2 to allow at least one matching
    best_template_key   = None
    best_template_mask_key = None
    best_ang = N.zeros(3)
    best_loc = N.zeros(3)

    for i,(t_key, tm_key) in enumerate(vm_tem):
        t = iv.get_mrc(t_key)      # we do not cache templates, to avoid confusion among multiple workers
        tm = iv.get_mrc(tm_key)    # we do not cache templates, to avoid confusion among multiple workers

        if tem_cutoff!=None:   tm = N.array((tm > tem_cutoff), dtype=N.float, order='F')

        at_re = align_to_template(v=v, vm=vm, t=t, tm=tm, align_op=align_op, logger=logger)
        if N.isnan(at_re['score']):
            if logger is not None:      logger.warning("alignment failed: v_key %s, vm_key %s, t_key %s, tm_key %s, error %s ", v_key, vm_key, t_key, tm_key, repr(at_re['err']) )



        if at_re['score'] > best_score:
            best_score  = at_re['score']
            best_ang    = at_re['angle']
            best_loc    = at_re['loc']
            best_template_key = t_key
            best_template_mask_key = tm_key


    re = {'vol_key':v_key, 'vol_mask_key':vm_key, 'angle':best_ang, 'loc':best_loc, 'score':best_score, 'template_key':best_template_key, 'template_mask_key':best_template_mask_key}

    if transform:
        vr = tomo.rotate_vol_pad_mean(v, best_ang, best_loc)
        vmr = tomo.rotate_mask(vm, best_ang)

        re['vol_r_key'] = v_key + '.aligned.mrc'
        re['vol_mask_r_key'] = vm_key + '.aligned.mrc'

        iv.put_mrc(vr, re['vol_r_key'])
        iv.put_mrc(vmr, re['vol_mask_r_key'])

                
    return re

'''
    # compute normalized volume for all entries in first volume list.
    for vk,mk,ang,loc in vmal_1_in:

        v = iv.get_mrc_cache_fs(vk)
        m = iv.get_mrc_cache_fs(mk)

        v = tomo.rotate_vol_pad_mean(v,ang,loc)
        m = tomo.rotate_mask(m, ang)

        v -= vol_avg
        v = N.real(ifftn(ifftshift(fftshift(fftn(v)) * m)))
        v -= N.mean(v)

        vols1.append(v)
        masks1.append(m)

    # compute normalized volume for all entries in second volume list.
    for vk,mk,ang,loc in vmal_2_in:

        v = iv.get_mrc_cache_fs(vk)
        m = iv.get_mrc_cache_fs(mk)

        v = tomo.rotate_vol_pad_mean(v,ang,loc)
        m = tomo.rotate_mask(m, ang)

        v -= vol_avg
        v = N.real(ifftn(ifftshift(fftshift(fftn(v)) * m)))
        v -= N.mean(v)

        vols2.append(v)
        masks2.append(m)

    # calculate pairwise scoring. which is simple dot product for normalized values.
    for i, v1 in enumerate(vols1):
        for j, v2 in enumerate(vols2):
            mat[i,j] = N.sum(v1*v2)

    # save output back to disk.
    iv.np_save(mat_out_key, mat)

    return True
'''


def align_to_templates__parallel(vmal_org=None, vm_tem=None, align_op=None, runner=None, n_chunk=None, tem_cutoff=None, transform=False):

    vmal = [_ for _ in vmal_org]
    inds = range(len(vmal_org))

    tasks = []
    while vmal:
        vmal_t = vmal[:n_chunk]
        inds_t = inds[:n_chunk]
        tasks.append(runner.task('align_to_templates__local', vmal_org=vmal_t, inds=inds_t, vm_tem=vm_tem, align_op=align_op, tem_cutoff=tem_cutoff, transform=transform))
        vmal = vmal[n_chunk:]
        inds = inds[n_chunk:]

    ress = [_ for _ in runner.run__except(tasks)]

    # record alignment result in the order of vmal_org
    align = [None] * len(vmal_org)
    for res in ress:
        rec = res.result
        for i, ind_t in enumerate(rec['inds']):
            assert (vmal_org[ind_t][0] == rec['align'][i]['vol_key'])
            align[ind_t] = rec['align'][i]

    return align

def align_to_templates__local(vmal_org=None, inds=None, vm_tem=None, align_op=None, tem_cutoff=None, transform=False, logger=None, tmp_dir=None):

    assert  len(vmal_org) == len(inds)

    align = [None] * len(vmal_org)
    for i, vmal_t in enumerate(vmal_org):
        at_re = align_to_templates(vmal=vmal_t, vm_tem=vm_tem, align_op=align_op, tem_cutoff=tem_cutoff, transform=transform, logger=logger, tmp_dir=tmp_dir)
        align[i] = at_re

    return {'align':align, 'inds':inds}


# align one volume against a set of templates, and transform according to the best template
# parameters: vm volume and its mask, vm_tem templates and their masks
def align_and_transform_vols_against_templates(vm, vm_tem, L):

    v_key, m_key = vm
    vmal = (v_key, m_key, N.zeros(3), N.zeros(3))

    re = align_to_templates(vmal, vm_tem, L, transform=True)
    
    
    return re

#-----------------------------------------------------------------------
# align with segmentation


def segment_and_align_vols_against_one_template__parallel(vmal, vm_tem, gauss_sigma, L):
    vmal = [_ for _ in vmal_org]
    inds = range(len(vmal_org))

    tasks = []
    while vmal:
        vmal_t = vmal[:n_chunk]
        inds_t = inds[:n_chunk]
        tasks.append(runner.task('segment_and_align_vols_against_one_template', vmal=vmal_t, inds=inds_t, vm_tem=vm_tem, gauss_sigma=gauss_sigma, L=L))
        vmal = vmal[n_chunk:]
        inds = inds[n_chunk:]

    ress = [_ for _ in runner.run__except(tasks)]

    # record alignment result in the order of vmal_org
    align = [None] * len(vmal_org)
    for res in ress:
        rec = res.result
        for i, ind_t in enumerate(rec['inds']):
            assert (vmal_org[ind_t][0] == rec['align'][i]['vol_key'])
            align[ind_t] = rec['align'][i]

    return align



# align with segmentation:
# apply gaussian smoothing to both template and subtomogram, then use template to guide segmentation of a subtomogram
# then align the segmented subtomogram to the original template 
def segment_and_align_vols_against_one_template(vmal=None, inds=None, vm_tem=None, gauss_sigma=None, L=None, cutoff_num=100, tmp_dir=None):
    import tomominer.filter.gaussian as fg
    from tomominer.segmentation.target_complex import Supervised as segts
    import tomominer.segmentation.thresholding as segt

    tk, tmk = vm_tem

    t = iv.get_mrc(tk)      # template, we do not cache
    tm = iv.get_mrc(tmk)    # mask, we do not cache

    tg = fg.smooth(v=t, sigma=gauss_sigma)
    mmb = segt.ms_mo_best( segt.ms_cutoffs(v=tg, cutoffs=N.linspace(0, tg.max(), cutoff_num)) )
    t_seg =     N.array( (tg > mmb['cutoff']), dtype=N.float, order='F' ) 


    align = [None] * len(vmal)
    for i, (vk, mk, ang, loc) in enumerate(vmal):
        v = iv.get_mrc_cache_fs(vk, tmp_dir)      # subtomogram
        vm = iv.get_mrc_cache_fs(vmk, tmp_dir)    # mask
        vg = fg.smooth(v=v, sigma=gauss_sigma)

        a_re = segment_and_align_one_vol_against_one_template(t=t, tg=tg, t_seg=t_seg, tm=tm, v=v, vg=vg, vm=vm, L=L)
        align[i] = {'vol_key':vk, 'vol_mask_key':mk, 'angle':a_re['ea'], 'loc':a_re['dx']}

    return {'align':vmal_new, 'inds':inds}



# parameters: t:template,   tg:smoothed template,   t_seg:segmented template region,    tm:template mask,   v:subtomogram,  vg:smoothed subtomogram,    vm:
def segment_and_align_one_vol_against_one_template(t, tg, t_seg, tm, v, vg, vm, L):

    a_vt_re = align_vols(vg, vm, tg, tm, L)       # align template against subtomogram in order to segment the subtomogram
    if not N.isfinite(a_vt_re['score']) is None:     return

    t_seg_r = tomo.rotate_vol_pad_zero(t_seg, a_vt_re['angle'], a_vt_re['loc'])
    t_seg_r = (t_seg_r>0.5)

    # find the target complex region and mask out the rest
    sc_re = segts.simple_cut(msk=t_seg_r, vol=vg)
    if sc_re is None:       return

    v_seg = N.copy(v, order='F')
    v_seg[sc_re['mask']==0] = v_seg[sc_re['mask']==1].mean()

    a_tv_re = align_vols(t, tm, v_seg, vm, L)
    if not N.isfinite(a_vt_re['score']) is None:     return

    re = a_tv_re
    re['t_seg_r'] = t_seg_r
    re['sc_re'] = sc_re
    re['v_seg'] = v_seg

    return re
    



