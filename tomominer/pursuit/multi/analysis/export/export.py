
# given a pursuit dump, export information that are needed for downstream analysis or inspection
import os
import shutil

import numpy as np

import tomominer.io.file as io
import tomo

def export_dump(d, out_dir):
    print out_dir

    if not os.path.exists(out_dir):     os.makedirs(out_dir)
   
    f_t = d['gavg_re']['vol_avg_out_key']
    fo_t = os.path.join(out_dir, os.path.basename(f_t) )      
    if not os.path.exists(fo_t):        shutil.copyfile( f_t, os.path.join(out_dir, os.path.basename(f_t) ) )

    f_t = d['gavg_re']['mask_avg_out_key']
    fo_t = os.path.join(out_dir, os.path.basename(f_t) )
    if not os.path.exists(fo_t):          shutil.copyfile( f_t, os.path.join(out_dir, os.path.basename(f_t) ) )

    fo_t = os.path.join(out_dir,'cov_avg_for_masking.mrc')
    if not os.path.exists(fo_t):
        cov_avg = np.array(d['cfp_re']['cov_avg_for_masking'], order='F')
        io.put_mrc(cov_avg, fo_t) 

    fo_t = os.path.join(out_dir,'cov_avg_for_masking_mask.mrc')
    if not os.path.exists(fo_t):
        cov_avg_m = np.zeros(cov_avg.size)
        cov_avg_m[d['cfp_re']['voxel_mask_inds']] = 1
        cov_avg_m = np.reshape(cov_avg_m, cov_avg.shape, order='F')
        io.put_mrc(cov_avg_m, fo_t) 

    #import pdb;     pdb.set_trace()


