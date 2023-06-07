#!/usr/bin/env python

'''
given a density map, and specification of CTF, convolute map with CTF, then generate corresponding template and a wedge mask
'''

import os, json, copy
import numpy as N
from numpy.fft import fftn, ifftn, fftshift, ifftshift
import tomominer.io.file as IF
import tomominer.image.optics.ctf as IOC
import tomominer.image.vol.wedge.util as IVWU
import tomominer.model.util as MU
import tomominer.image.vol.util as IVU
import tomominer.geometry.rotate as GR

def convert(op):
    extension = 'map.mrc'
    pdb_path = {}
    for root, sub_folders, files in os.walk(op['pdb_dir']):
        for file_t in files:
            if not file_t.endswith(extension):  continue
            pdb_id = file_t[: len(file_t) - len(extension)]
            assert (pdb_id + extension) == file_t
            assert pdb_id not in pdb_path
            pdb_path[pdb_id] = os.path.join(root, file_t)

    ctf_op = op['ctf']
    if 'pix_size' not in ctf_op:
        with open(op['map_op_file']) as f:       map_op = json.load(f)
        ctf_op['pix_size'] = map_op['spacing'] / 10.0
    
    print 'ctf_op', ctf_op
    #just read on map file for shape/size
    v = IF.read_mrc_vol(pdb_path[pdb_id]) # this pdb_id comes from last iteration of above for loop
    if 'output_map_size' in op:
        print "output map size", op['output_map_size']
    else:
        print 'output map size', v.shape

    ctf_t = IOC.create(Dz=ctf_op['Dz'], size=v.shape, pix_size=ctf_op['pix_size'], voltage=ctf_op['voltage'], Cs=ctf_op['Cs'], sigma=(None if 'sigma' not in ctf_op else ctf_op['sigma']))['ctf']
    
    for pdb_id in pdb_path:
        ctf = copy.deepcopy(ctf_t)
        print pdb_id[0:4]
        v = IF.read_mrc_vol(pdb_path[pdb_id])
        v = N.real(     ifftn(  ifftshift(  ctf * fftshift(fftn(v)) ) )     )           # convolute v with ctf
        del ctf
        if 'output_map_size' in op:
            print 'resizing'
            v = IVU.resize_center(v=v, s=op['output_map_size'], cval=v.mean())
            
        template_file_out = "./" + pdb_id + "template.mrc"        
        if op['reverse_intensity']: v = -v
        IF.put_mrc(v, path=template_file_out, overwrite=True)
        
        if op['create_masks']:
            mask_file = "./" + pdb_id + "mask.mrc"
            mask = MU.sphere_mask(shape=v.shape, center=IVU.fft_mid_co(siz=v.shape))
            IF.put_mrc(mask, path=mask_file, overwrite=True)


if __name__ == '__main__':
    with open('map_to_template__batch_op.json') as f:     op = json.load(f)  
    convert(op)
