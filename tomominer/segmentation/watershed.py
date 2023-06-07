
import numpy as np

import tomominer.core as core






#-----------------------------------------------------------------------
# c++ version of segmentation, using cython interface

def segmentation_cpp(vol_map, vol_lbl, max_overall_voxel_num=None, max_segment_voxel_num=None, queue_label=-1, conflict_lbl=-2):
    vol_map = np.array(vol_map, dtype=np.double, order='F')
    vol_lbl = np.array(vol_lbl, dtype=np.int32, order='F')

    if max_overall_voxel_num==None :   max_overall_voxel_num = vol_lbl.size + 1
    if max_segment_voxel_num==None :   max_segment_voxel_num = vol_lbl.size + 1
        

    return core.watershed_segmentation(vol_map, vol_lbl, max_overall_voxel_num, max_segment_voxel_num, queue_label, conflict_lbl)


'''

# test code



import os
import sys
sys.path.append(os.path.join( os.getenv('HOME'), 'ln/tomominer/tomominer' ))
sys.path.append(os.path.join( os.getenv('HOME'), 'ln/tomominer/build/lib.linux-x86_64-2.7/tomominer' ))


import numpy as N

v = N.random.random( (50,50,50) )
import tomominer.filter.gaussian as FG
v = FG.smooth(v, sigma=5)
v = N.array(v, order='F')

import tomominer.filter.local_extrema as FL
w = FL.local_maxima(v)
s = N.zeros(v.shape, dtype=int, order='F')

for i in range(len(w[0])):
    s[w[0][i], w[1][i], w[2][i]] = i + 1

import tomominer.segmentation.watershed as SW
seg = SW.segmentation_cpp(v, s)

import matplotlib
matplotlib.use('Qt4Agg')

import tomominer.image.vol.util as CV
CV.dspcub(v)
CV.dspcub(N.array(seg['vol_seg_lbl'], dtype=N.float))


import matplotlib.pyplot as plt
plt.show()



'''


#-------------------------------
# proxy

def segment(*args,**kwargs):
    return segmentation_cpp(*args,**kwargs)

    


