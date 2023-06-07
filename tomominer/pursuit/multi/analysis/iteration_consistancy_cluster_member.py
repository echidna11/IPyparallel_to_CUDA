#!/usr/bin/env python


# given classification iterations, inspect the cluster member consistancy between consective iterations

def label_consistancy(d0, d1):

    ld0 = { _['subtomogram']:int(_['cluster_label']) for _ in d0 }
    ld1 = { _['subtomogram']:int(_['cluster_label']) for _ in d1 }
    
    ss = ld0.keys()
    l0 = [ld0[s] for s in ss]
    l1 = [ld1[s] for s in ss]

    from sklearn import metrics
    
    return metrics.adjusted_rand_score(l0, l1), metrics.adjusted_mutual_info_score(l0, l1), metrics.homogeneity_score(l0, l1), metrics.completeness_score(l0, l1)


if __name__ == '__main__':

    import os
    import json

    pass_i = 0
    while True:
        dj0f  = os.path.join(os.getcwd(), 'pass_%03d' % (pass_i), 'data_config.json')
        dj1f  = os.path.join(os.getcwd(), 'pass_%03d' % (pass_i + 1), 'data_config.json')

        if not os.path.exists(dj0f) :       break
        if not os.path.exists(dj1f) :       break

        with open(dj0f) as f:   dj0 = json.load(f)
        with open(dj1f) as f:   dj1 = json.load(f)

        print pass_i, label_consistancy(dj0, dj1)

        pass_i += 1

