import os

import numpy as np

import tomominer.io.file as iv

import tomominer.filter.gaussian as filterg
from tomominer.segmentation.target_complex import Supervised as segts
import tomominer.align.util as align


import tomo




def save_figure(fname, op):

    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    f, sp = plt.subplots(2,2)

    import general_util.vol as gv

    sp[0,0].imshow( gv.cub_img(v=op['v'])['im'], cmap = cm.Greys_r )
    sp[0,1].imshow( gv.cub_img(v=op['seg'])['im'], cmap = cm.Greys_r )
    sp[1,0].imshow( gv.cub_img(v=op['vg'])['im'], cmap = cm.Greys_r )
    sp[1,1].imshow( gv.cub_img(v=op['vg_seg'])['im'], cmap = cm.Greys_r )

    for i in range(sp.shape[0]):
        for j in range(sp.shape[1]):
            sp[i,j].axis('off') # clear x- and y-axes

    plt.draw()

    import pylab
    pylab.savefig(fname)
    
    plt.close("all")






# parameters: t: template, tm: template mask, t_seg: template segmentation
def segment(vol_key, t, tm, t_seg, op, cache_dir):



    vk, vmk = vol_key
    v = iv.get_mrc_cache_fs(vk, cache_dir)
    vm = iv.get_mrc_cache_fs(vmk, cache_dir)

    vg = filterg.smooth(v=v, sigma=op['smooth']['gauss_sigma']);        vg = np.array(vg, order='F')

    a_vt_re = align.align_vols(vg, vm, t, tm, op['align']['L'])       # align template against subtomogram in order to segment the subtomogram
    t_seg_r = tomo.rotate_vol_pad_zero_py(t_seg, a_vt_re['ang'], a_vt_re['loc'])
    t_seg_r = (t_seg_r>0.5)

    # find the target complex region and mask out the rest
    sc_re = segts.simple_cut(msk=t_seg_r, vol=vg)
    v_seg = np.copy(v, order='F')
    v_seg[sc_re['mask']==0] = v_seg[sc_re['mask']==1].mean()        # replacing with mean value


    out_key = op['out_dir'] + vk;       assert(out_key != vk);      # use assertion to prevent overwritten of original file
    out_key_dir = os.path.dirname(out_key)

    if not os.path.isdir(out_key_dir):
        try:
            os.makedirs(out_key_dir)
        except:
            pass


    if not os.path.isfile(out_key): 
        iv.put_mrc(v_seg, out_key)

    out_png_key = out_key+'.png'
    if not os.path.isdir(out_png_key):
        vg_seg = np.copy(vg);       vg_seg[sc_re['mask']==0] = vg_seg[sc_re['mask']==1].mean()
        save_figure(out_png_key, {'v':v, 'vg':vg, 'vg_seg':vg_seg, 'seg':sc_re['mask']})

    return {'vol_key':out_key, 'mask_key':vmk}


    
