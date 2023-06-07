#!/usr/bin/env python

'''
from particle picking result, generate an imod model file so that the imod can be used to display peaks

IMPORTANT: Because the particles are normally ELONGATED along z axis, the proper subtomogram size can be measured by inspect x-z and y-z planes. To do so, set view_direction=0, and choose proper radius_factor to see if the markers can cover the particles in x-z and y-z plane!!

'''


import tomominer.pursuit.multi.analysis.localization.particle_location_display_imod as PMALP


import numpy as N


def main():

    import os
    import subprocess
    import uuid

    import json
    with open('particle_picking_dog__display__imod__config.json') as f:       op = json.load(f)
    if 'view_direction' not in op:      op['view_direction'] = 0

    with open(op['particle_picking_dog__config_file']) as f:       pp_op = json.load(f)

    with open(op['peak_file']) as f:       pp = json.load(f)

    with open(pp_op['tomogram_info_file']) as f:    tom_info = json.load(f)
    tom_info = {  _['id']:_ for _ in tom_info }



    # select only top number of peaks with highest values
    pp = [_ for _ in pp if (_['tomogram_id'] == op['tomogram_id'])]
    if pp_op['find_maxima']:
        pp = sorted(pp, key=lambda _:(-_['val']))      # sort pp according to decrease of peak value
    else:
        pp = sorted(pp, key=lambda _:(_['val']))      # sort pp according to increase of peak value

    if 'top_num' in op:     pp = pp[:op['top_num']]
    assert      len(pp) > 0

    # Modification for displaying top or bottom few particles
    if op['top_bottom']==0:
        pp.reverse()

    # End modification

    sigma1 = pp_op['sigma1'] / tom_info[op['tomogram_id']]['voxel_spacing']

    rad = sigma1
    if 'radius_factor' in op:   rad *= op['radius_factor']

    x = N.array([_['x'] for _ in pp])


    if 'projection_slices' in op:
        p_start = op['projection_slices']['start']
        p_end = op['projection_slices']['end']
        p_slice = op['projection_slices']['slice']

        if op['view_direction'] == 0:
            p_dim = 2
        elif op['view_direction'] == 1:
            p_dim = 1
        elif op['view_direction'] == 2:
            p_dim = 0

        if p_end < 0:   p_end = x[:,p_dim].max()

        print 'project all peaks within a range of slices', p_start, p_end, 'onto one slice', p_slice, ', just for demostration purpose'
        x = x[N.logical_and((x[:,p_dim]>=p_start), (x[:,p_dim]<=p_end)), :]
        x[:,p_dim] = p_slice


    map_file = tom_info[op['tomogram_id']]['vol_file']
    clip_file = None

    if 'clip' in op:
        print 'cut out a sub-region for faster display'

        clip_file = os.path.join('/tmp', str(uuid.uuid1())+'.mrc')

        subprocess.call(['clip', 'resize', '-cx', str((op['clip']['start'][0]+op['clip']['end'][0])/2), '-cy', str((op['clip']['start'][1]+op['clip']['end'][1])/2), '-cz', str((op['clip']['start'][2]+op['clip']['end'][2])/2), '-ox', str(op['clip']['end'][0]-op['clip']['start'][0]), '-oy', str(op['clip']['end'][1]-op['clip']['start'][1]), '-oz', str(op['clip']['end'][2]-op['clip']['start'][2]), '-m', 'short', map_file, clip_file])

        # remove all points outside the clip
        in_clip = N.zeros(x.shape[0], dtype=N.bool)
        in_clip[:] = True
        for dim in range(3):
            in_clip[x[:,dim] < op['clip']['start'][dim]] = False
            in_clip[x[:,dim] >= op['clip']['end'][dim]] = False

        x = x[in_clip, :]

        for dim in range(3):        x[:, dim] -= op['clip']['start'][dim]



    l = PMALP.generate_lines(x_full=x, rad=rad, view_direction=op['view_direction'])

    PMALP.display_map_with_lines(l=l, map_file=(map_file if clip_file is None else clip_file), remove_intermediate_file=False)
    
    if clip_file is not None:       os.remove(clip_file)


if __name__ == '__main__':
    main()

