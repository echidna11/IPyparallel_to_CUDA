#!/usr/bin/env python

'''
given a density map, and specification of CTF, convolute map with CTF, then generate corresponding template and a wedge mask

~/ln/tomominer/tomominer/template/generate/map_to_template.py
'''


import os, json

import numpy as N
from numpy.fft import fftn, ifftn, fftshift, ifftshift

import tomominer.io.file as IF
import tomominer.image.optics.ctf as IOC
import tomominer.image.vol.wedge.util as IVWU
import tomominer.model.util as MU
import tomominer.image.vol.util as IVU
import tomominer.geometry.rotate as GR



def convert(op):
    
    v = IF.read_mrc_vol(op['map_file'])

    ctf_op = op['ctf']

    if 'pix_size' not in ctf_op:
        with open(op['map_op_file']) as f:       map_op = json.load(f)
        ctf_op['pix_size'] = map_op['spacing'] / 10.0

    print 'ctf_op', ctf_op

    ctf = IOC.create(Dz=ctf_op['Dz'], size=v.shape, pix_size=ctf_op['pix_size'], voltage=ctf_op['voltage'], Cs=ctf_op['Cs'], sigma=(None if 'sigma' not in ctf_op else ctf_op['sigma']))['ctf']

    v = N.real(     ifftn(  ifftshift(  ctf * fftshift(fftn(v)) ) )     )           # convolute v with ctf

    if 'output_map_size' in op:
        print 'resizing'
        v = IVU.resize_center(v=v, s=op['output_map_size'], cval=v.mean())

    print 'output map size', v.shape

    IF.put_mrc(((-v) if op['reverse_intensity'] else v), path=op['template_file_out'], overwrite=True)


    if 'mask_file_out' in op:

        mask = MU.sphere_mask(shape=v.shape, center=IVU.fft_mid_co(siz=v.shape))
        IF.put_mrc(mask, path=op['mask_file_out'], overwrite=True)



if __name__ == '__main__':
    with open('map_to_template__op.json') as f:     op = json.load(f)
    
    convert(op)



