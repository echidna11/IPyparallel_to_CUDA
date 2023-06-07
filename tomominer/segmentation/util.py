# utility functions for segmentation


import scipy.ndimage as SN

def connected_components(m):
    l, n = SN.label(m)

    return {'label':l, 'num':n}



'''

# test code

import numpy as N
v = N.zeros((5,5,5))
v[0:2, 0:2, 0:2] = 1
v[3:5, 3:5, 3:5] = 1

import tomominer.segmentation.util as SU
r = SU.connected_components(v)

print v
print r


'''
