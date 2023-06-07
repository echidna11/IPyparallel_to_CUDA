


import numpy as np
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology

def local_maxima(arr):
    # http://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    """
    Takes an array and detects the troughs using the local maximum filter.
    Returns a boolean mask of the troughs (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    # define an connected neighborhood
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#generate_binary_structure
    neighborhood = morphology.generate_binary_structure(len(arr.shape),2)
    # apply the local maximum filter; all locations of maximal value 
    # in their neighborhood are set to 1
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.filters.html#maximum_filter
    local_max = (filters.maximum_filter(arr, footprint=neighborhood)==arr)
    # local_max is a mask that contains the peaks we are 
    # looking for, but also the background.
    # In order to isolate the peaks we must remove the background from the mask.
    # 
    # we create the mask of the background
    background = (arr==arr.min())           # mxu: in the original version, was         background = (arr==0)
    # 
    # a little technicality: we must erode the background in order to 
    # successfully subtract it from local_max, otherwise a line will 
    # appear along the background border (artifact of the local maximum filter)
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#binary_erosion
    eroded_background = morphology.binary_erosion(
        background, structure=neighborhood, border_value=1)
    # 
    # we obtain the final mask, containing only peaks, 
    # by removing the background from the local_max mask
    #detected_maxima = local_max - eroded_backround             # mxu: this is the old version, but the boolean minus operator is deprecated
    detected_maxima = np.bitwise_and(local_max, np.bitwise_not(eroded_background))          # Material nonimplication, see http://en.wikipedia.org/wiki/Material_nonimplication
    return np.where(detected_maxima)





def local_minima(arr):
    return local_maxima(-arr)




'''

# testing code

# randomly place several ones inside a 3D map, then see if detected local maxima are at corresponding places




import os,sys
sys.path.append(os.path.join( os.getenv('HOME'), 'ln/tomominer/tomominer' ))



siz = [64,64, 64]

import numpy as N
v = N.zeros(siz) + 0.00001

find_maxima=False

x = []
for i in range(5):
    xt = [  N.random.randint(siz[0]),  N.random.randint(siz[1]), N.random.randint(siz[2])  ]
    if find_maxima:
        v[xt[0], xt[1], xt[2]] = 1.0
    else:
        v[xt[0], xt[1], xt[2]] = -1.0
    x.append(xt)

import tomominer.filter.local_extrema as FLE
if find_maxima:
    x_p = FLE.local_maxima(v)
else:
    x_p = FLE.local_minima(v)

print N.array(x)
print N.array(x_p).T






'''


