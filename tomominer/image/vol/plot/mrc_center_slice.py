#!/usr/bin/env python


"""
Plot center slice of a number of mrc files, 

align them (against template / average) if necessary

"""

'''
~/ln/tomominer/tomominer/image/vol/plot/mrc_center_slice.py
'''


import os, sys
import numpy as N

import tomominer.io.file as IF
from tomominer.image.io import save_image_matplotlib

def main():

    v_threshold = float(sys.argv[1])
   
    for fn in os.listdir('./'):
        if not os.path.isfile(fn):      continue
        if not fn.endswith('.mrc'):     continue

        v = IF.read_mrc_vol(fn)
        v -= v.mean()
        v /= v.std()

        fn_root, fn_ext = os.path.splitext(fn)
        save_image_matplotlib(m=N.squeeze(v[v.shape[0]/2,:,:]), out_file='%s--0.png'%(fn_root,), vmin=-v_threshold, vmax=v_threshold)
        save_image_matplotlib(m=N.squeeze(v[:,v.shape[1]/2,:]), out_file='%s--1.png'%(fn_root,), vmin=-v_threshold, vmax=v_threshold)
        save_image_matplotlib(m=N.squeeze(v[:,:,v.shape[2]/2]), out_file='%s--2.png'%(fn_root,), vmin=-v_threshold, vmax=v_threshold)





if __name__ == '__main__':
    main()
