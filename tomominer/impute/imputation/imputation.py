
# functions for imputation
import numpy as np
from numpy.fft import fftn, ifftn, fftshift, ifftshift

# assume subtomogram(vol) is aligned, normalize then fill up missing value region of vol
def impute_subtomogram(tem, vol, mask):
    vol_v = vol.flatten()
    vol_v_mean = np.mean(vol_v)
    vol_v_std = np.std(vol_v)
    
    vol_fft = fftshift(fftn(vol))
    
    tem_v = tem.flatten()
    tem_v_mean = np.mean(tem_v)
    tem_v_std = np.std(tem_v)

    tem = (tem - tem_v_mean) / tem_v_std
    tem = tem * vol_v_std + vol_v_mean

    
    tem_fft = fftshift(fftn(tem))
    
    ind =  (mask == 0)
    vol_fft[ind] = tem_fft[ind]
    
    vol = real(ifftn(ifftn(vol_fft))) 
     
    return vol
    

def rotate_subtomogram_then_impute(tem, vol, mask, ang, loc):
    vol_r = rotate(vol, ang, loc)
    mask_r = rotate(mask, ang)
    
    return impute_subtomogram(tem, vol_r, mask_r)
         
def rotate_template_then_impute(tem, vol, mask, ang, loc):
    ang_inv, loc_inv = inverse_rigid_transform(ang, loc)    # suppose ang and loc are rigid transforms applied to vol to match tem, need to get inverse transform in order to rotate tem to match loc

    tem_r = rotate(tem, ang_inv, loc_inv)
    return impute_subtomogram(tem_r, vol, mask)
