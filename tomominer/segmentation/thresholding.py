# -*- coding: utf-8 -*-

# functions for automatic thresholding

import numpy as np

import tomominer.statistics as sta


# given a volume (v), which is often smoothed,
# and a set of thredholds (cutoffs), calculate Mumford–Shah scores, 

def ms_cutoffs(v, cutoffs):

    s = {}
    for i, c in enumerate(cutoffs):
        # we need to seperate voxels to two classes
        assert np.any(v>c)
        assert np.any(v <= c)

        s_t = sta.mumford_shah_scores(v, (v>c))
        s[i] = s_t
        s[i]['cutoff'] = c


    return s


import numpy as np
# given results generated from ms_cutoffs(),  use multiobjective optimization to determine optimal threshold
def ms_mo_best(scores):
    
    x = np.zeros( [len(scores), 2] )
    for i in scores:
        s = scores[i]
        x[i, 0] = s['dif_mean_sq_sum']      # we want to minimize sum of difference
        x[i, 1] = s['bdry_area']            # we want to minimize boundary area

    import tomominer.optimization.multiobjective as omo
    pi = list(omo.pareto_min__set_sel(x))
    pmr_re = omo.pareto_min__regress(x[pi, :])

    pi_min = pi[np.argmin(pmr_re['score'])]

    return {'best':scores[ pi_min ], 'pi_min':pi_min, 'x':x, 'pi':pi, 'pmr_re':pmr_re}

