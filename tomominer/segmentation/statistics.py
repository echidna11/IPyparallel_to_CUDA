# -*- coding: utf-8 -*-

# different types of statistics for segmentation

import numpy as np
import tomominer.core as tomo

# given a volume (v), and a segmentation (m), calculate simplified Mumford–Shah scores, 
# the scores is defined according to paper chan01
def mumford_shah_scores(v, m):
    m = np.array(m, dtype=np.int32, order='F')
    mf = np.isfinite(v)
    

    dif_mean_sq_sum = 0
    for s in np.unique(m):
        vt = v[mf & (m==s)]
        dif_mean_sq_sum += sum( (vt - vt.mean())**2 )


    bdry_mask = np.zeros(m.shape, dtype=np.int32, order='F')
    tomo.segment_boundary(m, bdry_mask)

    return {'dif_mean_sq_sum':dif_mean_sq_sum, 'bdry_area':bdry_mask.sum()}

