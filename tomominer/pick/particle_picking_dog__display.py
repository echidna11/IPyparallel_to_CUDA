#!/usr/bin/env python


# display selected peaks for visual inspection


import matplotlib
matplotlib.use('Qt4Agg')


if __name__ == '__main__':

    import json
    with open('particle_picking_dog__display__config.json') as f:       op = json.load(f)

    with open(op['particle_picking_dog__config_file']) as f:       pp_op = json.load(f)

    with open(op['particle_picking_dog__out_file']) as f:       pp = json.load(f)

    with open(pp_op['tomogram_info_file']) as f:    tom_info = json.load(f)
    tom_info = {  _['id']:_ for _ in tom_info }



    print 'loading map', tom_info[op['tomogram_id']]['vol_file']
    import tomominer.io.file as IF
    v = IF.read_mrc(tom_info[op['tomogram_id']]['vol_file'])['value']

    import numpy as N 
    if 'clip_region' in op:
        # use this option to show only a subvolume, in order for fast processing
        cr = N.array(       op['clip_region']       )
        print 'clip_region', cr

        v = v[cr[0,0]:cr[0,1], cr[1,0]:cr[1,1], cr[2,0]:cr[2,1]]

        keep_flag = [True] * len(pp)
        for peak_i, pp_t in enumerate(pp):
            for dim_i in range(cr.shape[0]):
                pp_t['x'][dim_i] -= cr[dim_i,0]
                if pp_t['x'][dim_i] < 0:        keep_flag[peak_i] = False
                if pp_t['x'][dim_i] >= v.shape[dim_i]:      keep_flag[peak_i] = False
        pp = [pp[_] for _ in range(len(pp)) if keep_flag[_]]

    # select only top number of peaks with highest values
    pp = [_ for _ in pp if (_['tomogram_id'] == op['tomogram_id'])]
    if pp_op['find_maxima']:
        pp = sorted(pp, key=lambda _:(-_['val']))      # sort pp according to decrease of peak value
    else:
        pp = sorted(pp, key=lambda _:(_['val']))      # sort pp according to increase of peak value

    if 'tom_num' in op:     pp = pp[:op['top_num']]

    sigma1 = pp_op['sigma1'] / tom_info[op['tomogram_id']]['voxel_spacing']

    rad = sigma1
    if 'radius_factor' in op:   rad *= op['radius_factor']


    if op['do_smoothing']:
        import tomominer.filter.gaussian as FG
        v = FG.dog_smooth__large_map(v, s1=sigma1, s2=sigma1*pp_op['sigma_k'])


    print 'displaying', len(pp), 'markers'
    import tomominer.plot.marker_display  as PM
    PM.display(v=v, x=N.array([_['x'] for _ in pp]), r=rad, linewidth=op['linewidth'])

