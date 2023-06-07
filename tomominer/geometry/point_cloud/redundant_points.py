# functions to identify redundant points

import numpy as N
import scipy.spatial.distance as SSD

# identify redundant points (when two points are closer than a threshold r, the point with larger id will be identified as redundant)
def identify(x, r):
    d = SSD.pdist(x)
    d = SSD.squareform(d)

    n = x.shape[0]
    f = [True] * n
    for i0 in range(n):
        w = N.where(d[i0,:].flatten() < r)[0].tolist()
        for i1 in w:
            if i1 <= i0:    continue
            f[i1] == False

    return f




# use special hasing to partition points, then identify redundant points (when two points are closer than a threshold r, the point with larger id will be identified as redundant)
def identify__special_hashing(x, r, partition_num=10):

    ids = N.array(      list(range(x.shape[0]))     )






