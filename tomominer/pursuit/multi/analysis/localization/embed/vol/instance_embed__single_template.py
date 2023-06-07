#!/usr/bin/env python




'''
embed instances in form of isosurface, generated using template_isosurface.py
assuming that there is only one template, and all subtomograms are aligned to that template. This is for template matching or averaging


~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/vol/instance_embed__single_template.py

'''


import sys, json

import numpy as N

import tomominer.io.file as IF
import tomominer.geometry.ang_loc as GA
import tomominer.geometry.rotate as GR
import tomominer.image.vol.util as IVU


'''
convert a model to a density map

parameters:   size: container size / map size     dj: instances that aligned to template    vt: template
modified from tomominer.simulation.tomogram.model_to_density_map.get_map()

'''

def get_map(size, dj, vt):

    v = N.zeros(shape=size, dtype=N.float) + N.nan

    for i, d in enumerate(dj):
        rev_rm, rev_loc_r = GA.reverse_transform(GA.rotation_matrix_zyz(d['angle']), d['loc'])
        vr = GR.rotate( v=vt, rm=rev_rm, loc_r=rev_loc_r, default_val=N.nan )

        IVU.add_to_whole_map(whole_map=v, vol=vr, c=N.array(d['mid_co']))

        print '\r %06d     %3.3f'%(i, float(i)/len(dj)),          ;       sys.stdout.flush()

    background = N.abs(v[N.isfinite(v)]).max()*1.1      # we use extreme value as background, for the purpose of intensity based segmentation
    v[N.isnan(v)] = background

    print 'backgroud', background, 'v[0,0,0]', v[0,0,0], 'v_max', v.max(), 'v_min', v.min(), 'v_mean', v.mean()
    print 'vt_mean', vt.mean(), 'vt_max', vt.max(), 'vt_abs_max', N.abs(vt).max()

    return v



def main():
    with open('instance_embed__single_template__op.json') as f:     op = json.load(f)

    vt = IF.read_mrc_vol(op['template file'])

    with open(op['data file']) as f:     dj = json.load(f)       # Note: we assume the dj is already filtered with respect to a particular tomogram_id

    for d in dj:        d['mid_co'] = d['peak']['loc']

    v = get_map(size=op['map size'], dj=dj, vt=vt)

    IF.write_mrc_by_chunk(v, op['out volume file'])

if __name__ == '__main__':
    main()



