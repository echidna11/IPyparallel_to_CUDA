#!/usr/bin/env python


'''
Plot center slice of a number of mrc files, 

align them (against template / average) if necessary

~/ln/tomominer/tomominer/image/vol/plot/tomogram_center_slice.py
'''


import os, json
import numpy as N

import tomominer.io.file as IF
from tomominer.image.io import save_png


def main():
    
    with open('tomogram_center_slice__op.json') as f:     op = json.load(f)


    if op['tomogram'].endswith('.mrc'):   
        v = IF.read_mrc_vol(op['tomogram']).astype(N.float)
    elif op['tomogram'].endswith('.npy'):
        v = N.load(op['tomogram'])
    else:
        raise Exception('Unknown file type')

    v -= v.mean()
    v /= v.std()

    if 'threshold' in op:
        t = op['threshold']
        v[v < -t] = -t
        v[v > t] = t

    if 'c' in op:
        c = op['c']
    else:
        c = (N.array(v.shape) / 2.0).astype(N.int)

    if not os.path.isdir(op['out dir']):        os.makedirs(op['out dir'])

    save_png(N.squeeze(v[c[0],:,:]), os.path.join(op['out dir'], '0.png'))
    save_png(N.squeeze(v[:,c[1],:]), os.path.join(op['out dir'], '1.png'))
    save_png(N.squeeze(v[:,:,c[2]]), os.path.join(op['out dir'], '2.png'))


if __name__ == '__main__':
    main()






'''

related code

~/ln/tomominer/tomominer/image/vol/plot/mrc_center_slice.py

'''

