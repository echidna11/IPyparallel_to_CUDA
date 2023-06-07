# functions to support colocalization analysis

import numpy as N
import scipy.spatial.distance as SSD
from munkres import Munkres

# given two sets of points x and y, calculate the distance matrix between pairs of points, then calculate use Hungarian algorithm to find point matches
def point_pair_match_exclusive(x, y):

    d = SSD.cdist(x, y)

    m = Munkres()
    i = m.compute(d.tolist())
    d_f = [d[_[0], _[1]] for _ in i]

    return {'i':i, 'd':d_f}

# given two sets of points x and y, calculate the distance matrix between pairs of points, then for each point in x, find its closest point in y regardless whether y has been matched to another x or not
def point_pair_match_nonexclusive(x, y):

    d = SSD.cdist(x, y)

    mi = N.argmin(d, axis=1).flatten().tolist()

    i = [(_,mi[_]) for _ in range(len(mi))]
    d_f = [d[_[0], _[1]] for _ in i]

    return {'i':i, 'd':d_f}





