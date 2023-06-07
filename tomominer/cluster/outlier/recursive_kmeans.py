

'''
~/ln/tomominer/tomominer/cluster/outlier/recursive_kmeans.py
'''


import numpy as N
import sklearn.cluster as SC



'''
given a min cluster size, recursively perform kmeans clustering and remove the clusters that are smaller than that size

parameters:         x: data     k: cluster number for kmeans        n_jobs: number of parallel jobs        min_size: smallest cluster size to select       max_repeat: max number of repeats
'''


def cluster_with_small_set_removal(x, k=10, n_jobs=-1, n_init=10, verbose=False, large_cluster_number_proportion=0.99, min_size=10, max_repeat_num=N.inf):
    

    inds = N.array(range(x.shape[0]))
    rep_i = 0
    while rep_i < max_repeat_num:
        # perform kmeans clustering
        km = SC.KMeans(n_clusters=k, n_jobs=n_jobs, n_init=n_init, verbose=verbose)
        lbl_t = km.fit_predict(x[inds,:])

        # remove small clusters
        cluster_sizes = {}
        found_small_cluster = False
        for l in range(N.max(lbl_t) + 1):
            size_t = N.sum(lbl_t == l)
            cluster_sizes[l] = size_t
            if size_t >= min_size:       continue
            lbl_t[lbl_t == l] = -1
            found_small_cluster = True

        cluster_sizes_t = [cluster_sizes[_] for _ in cluster_sizes]
        cluster_sizes_t = N.array(      sorted(cluster_sizes_t, reverse=True)       )

        inds = inds[lbl_t >= 0]
        lbl_tf = lbl_t[lbl_t >= 0]

        print 'repeat', rep_i, 'selected sample num', len(inds), 'cluster sizes', cluster_sizes_t[cluster_sizes_t >= min_size]

        if not found_small_cluster:     break
        if N.sum(cluster_sizes_t >= min_size) >= (len(cluster_sizes_t) * large_cluster_number_proportion):      break
        if len(inds) < k:       break

        rep_i += 1


    # generate final labels according to inds and lbl_tf

    lbl = N.zeros(len(x), dtype=N.int) - 1
    if len(inds) == 0:      return None

    lbl_count = 0

    for l in range(N.max(lbl_tf) + 1):
        l_i = (lbl_tf == l)
        if sum(l_i) == 0:     continue
        for i in inds[l_i]:            lbl[i] = lbl_count

        lbl_count += 1


    return lbl



'''
# test code

import numpy as N

x_data = N.random.normal(size=(1000, 5))
x_outlier = N.random.normal(scale=10, size=(100, 5))
x = N.vstack(   (x_data, x_outlier)     )

import tomominer.cluster.outlier.recursive_kmeans as C
l = C.cluster_with_small_set_removal(x, k=10)


'''

