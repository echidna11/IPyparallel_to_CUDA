#!/usr/bin/env python



'''

change voxel size and generate new mrc volume

~/ln/tomominer/tomominer/image/vol/conversion/mrc_voxel_size_change.py

'''


import json

import numpy as N

from tomominer.io.file import read_mrc_vol, put_mrc
import tomominer.geometry.rotate as GR


def main():
    with open('image_vol_conversion_mrc_voxel_size_change__op.json') as f:      op = json.load(f)

    v = read_mrc_vol(op['input file'])

    siz_org = N.array(v.shape)
    siz = siz_org * (float(op['orginal voxel spacing']) / op['out voxel spacing'])
    siz = N.round(siz).astype(N.int)

    rm = N.eye(3) * float(op['out voxel spacing']) / op['orginal voxel spacing']
    
    if op['padding mode'] == 'zero':
        v = GR.rotate_pad_zero(v, rm=rm, siz2=siz)
    elif op['padding mode'] == 'mean':
        v = GR.rotate_pad_mean(v, rm=rm, siz2=siz)
    else:
        raise Exception('padding mode')

    
    put_mrc(v, op['output file'])


if __name__ == '__main__':
    main()

