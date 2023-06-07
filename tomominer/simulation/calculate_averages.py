import os, math, gc
import numpy as N
import tomominer.io.file as IF
import scipy.ndimage as SN
from scipy.ndimage.filters import gaussian_filter as SNFG
import warnings
warnings.filterwarnings("ignore")

def eulerAnglesToRotationMatrix(theta):
    R_x = N.array([[1,         0,                  0                   ],
                    [0,         math.cos(theta[0]), -math.sin(theta[0]) ],
                    [0,         math.sin(theta[0]), math.cos(theta[0])  ]])

    R_y = N.array([[math.cos(theta[1]),    0,      math.sin(theta[1])  ],
                    [0,                     1,      0                   ],
                    [-math.sin(theta[1]),   0,      math.cos(theta[1])  ]])

    R_z = N.array([[math.cos(theta[2]),    -math.sin(theta[2]),    0],
                    [math.sin(theta[2]),    math.cos(theta[2]),     0],
                    [0,                     0,                      1]])

    R = N.dot(R_z, N.dot( R_y, R_x ))
    del R_x, R_y, R_z
    return R

def cluster_averaging_local(self, data_json, cluster, gaussian_filter, return_key=True):
    vol_sum = None
    mask_sum = None
    m = IF.get_mrc(data_json[0]['mask'])
    size_m = N.array(m.shape, dtype = N.float)
    c_m = (size_m -1 ) / 2.0
    size_v = None
    c_v = None
    for i, rec in enumerate(data_json):
        if not os.path.isfile(rec['subtomogram']):  continue
        v = IF.get_mrc(rec['subtomogram'])
        if i==0:
            size_v = N.array(v.shape, dtype = N.float)
            c_v = (size_v -1 ) / 2.0
        angle = N.array(rec['model']['angle'])
        misalign_angle =  N.array(rec['model']['misalign_angle'])
        rm1 = eulerAnglesToRotationMatrix(angle)
        rm1 = rm1.transpose() # reverse the rotation angles
        rm2 = eulerAnglesToRotationMatrix(misalign_angle)
        # We need to rotate back the subtomogram to original position because angle was used to rotate map for simulation. Also we need to rotate further for misalignment. Similar with mask.
        # Rotation matrix for both map and mask is same
        rm = N.matmul(rm1, rm2)
        # Rotate subtomogram
        c = -rm.dot(c_v) + c_v
        vr = SN.interpolation.affine_transform(input=v, matrix=rm, offset=c, output_shape=size_v.astype(N.int), cval=0.0)
        # Rotate mask
        c = -rm.dot(c_m) + c_m
        vm = SN.interpolation.affine_transform(input=m, matrix=rm, offset=c, output_shape=size_m.astype(N.int), cval=0.0)
        
        vr = SNFG(vr, gaussian_filter)
        if vol_sum is None:
            vol_sum = N.zeros(vr.shape, dtype=N.float64, order='F')
        vol_sum += vr
        if mask_sum is None:
            mask_sum = N.zeros(vm.shape, dtype=N.float64, order='F')
        mask_sum += vm
        del v, rm1, rm2, rm, vr, vm, c, angle, misalign_angle
    
    del m, size_m, size_v, c_v, c_m
    gc.collect()
    re = {'vol_sum':vol_sum, 'mask_sum':mask_sum, 'vol_count':len(data_json), 'cluster':cluster}
    if return_key:
        re_key = self.cache.save_tmp_data(re, fn_id=self.task.task_id)
        assert re_key is not None
        return {'key':re_key}
    else:
        return re