
# functions to segment out target complexes

import numpy as np

import tomominer.core as tomo
import tomominer.segmentation.watershed as segw
import tomominer.segmentation.thresholding as segt

import scipy.ndimage.morphology as snm

# functions for supervised segmentation, guided by a template
class Supervised:

    # suppose the (low pass filtered) template is aligned with a (low pass filtered) subtomogram vol
    # and a mask (msk) is generated from tomominer.template to indicate target complex region
    # define a cutoff to define structural region, 
    # then find those strucural regions that do not interset with the masked region, as neighbor regions
    # then use watershed segmentation to remove the structural regions
    @staticmethod
    def simple_cut(msk, vol, strucural_region_cutoff_ratios=(0.1, 0.9), cutoff_num = 20):
        
        assert min(strucural_region_cutoff_ratios) > 0
        assert max(strucural_region_cutoff_ratios) < 1 

        max_intensity = vol[msk > 0].max();         assert  (max_intensity > vol.mean())
        cutoffs = np.linspace(max(vol.mean(), max_intensity * min(strucural_region_cutoff_ratios)), max_intensity * max(strucural_region_cutoff_ratios), cutoff_num)
        mmb = segt.ms_mo_best( segt.ms_cutoffs(v=vol, cutoffs=cutoffs) )
        cutoff = mmb['best']['cutoff']


        stru_rgn_msk = (vol > cutoff)
        stru_rgn_msk_e = snm.binary_erosion(stru_rgn_msk)     # a little bit of erosion to elimate small pieces

        if ( (stru_rgn_msk_e & msk).sum() > 0 ):    stru_rgn_msk = stru_rgn_msk_e       # just make sure that after erosion there is still a valid region for the target complex, otherwise, discard erosion

        stru_rgn_msk = np.array(stru_rgn_msk, dtype=np.int32, order='F')

        
        conn_rgn_lbl = np.zeros(msk.shape, dtype=np.int32, order='F')
        tomo.connected_regions(stru_rgn_msk, conn_rgn_lbl)

        re = {'stru_rgn_msk':stru_rgn_msk, 'conn_rgn_lbl':conn_rgn_lbl, 'max_intensity':max_intensity, 'cutoff':cutoff, 'mmb':mmb }

        if conn_rgn_lbl.max() > 1:
            lbl_target = 1
            lbl_nei = 2
            stru_lbl = np.zeros(msk.shape, dtype=np.int32, order='F')
            for lbl in range(1, conn_rgn_lbl.max() + 1):
                msk_intersect =    (conn_rgn_lbl==lbl) & (msk > 0)

                if msk_intersect.sum() > 0:
                    stru_lbl[conn_rgn_lbl==lbl] = lbl_target
                else:
                    stru_lbl[conn_rgn_lbl==lbl] = lbl_nei


            ws_re = segw.segmentation_cpp(vol_map=vol, vol_lbl=stru_lbl)

            re['stru_lbl'] = stru_lbl
            re['mask'] = (ws_re['vol_seg_lbl']==lbl_target)

        else:
            re['stru_lbl'] = conn_rgn_lbl
            re['mask'] = np.ones(msk.shape)

        return re

