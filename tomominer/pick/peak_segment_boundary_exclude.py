#!/usr/bin/env python

# peak_segment_boundary_exclude.py


# exclude those peaks that are close to the boundary of segmentation. This is just an adhoc way to remove those peaks that corresponds to membrane.


if __name__ == '__main__':

    import json
    with open('peak_segment_boundary_exclude__config.json') as f:     op = json.load(f)

    dist_cutoff_voxel = float(op['dist_cutoff_voxel'])

    with open(op['peak_file']) as f:    pp = json.load(f)

    with open(op['tomogram_info_file']) as f:    tom_info = json.load(f)

    import pickle
    import numpy as N
    from scipy.spatial.distance import cdist
    import tomominer.geometry.surface.amira.util as GFAU
 
    pp_new = []
    for tom in tom_info:
        print 'tomogram', tom['id']

        pp_tom = [_ for _ in pp if (_['tomogram_id'] == tom['id'])]
        pp_tom = sorted(pp_tom, key=lambda _ : _['id'])         # order peaks according to ids,  assume the peak ids are ordered so that the smaller id has better value / importance

        with open(tom['header_file'], 'rb') as f:   mrc = pickle.load(f)['header']['MRC']

        with open(tom['segment_file']) as f:  surf = GFAU.surf_parse(f)
        vert_org = surf['vertices']

        # transfer vertice location from original space to the voxel grid space
        vert = N.zeros(vert_org.shape)
        vert[:,0] = (   (vert_org[:,0] - mrc['xorg']) / mrc['xlen']   ) * mrc['nx']
        vert[:,1] = (   (vert_org[:,1] - mrc['yorg']) / mrc['ylen']   ) * mrc['ny']
        vert[:,2] = (   (vert_org[:,2] - mrc['zorg']) / mrc['zlen']   ) * mrc['nz']

        pp_new_t = []
        for i, p in enumerate(pp_tom):
            print '\r', float(i)/len(pp_tom), '          ',
            x = N.array(p['x'], dtype=N.float).reshape(     (1,3)      )
            d = cdist(vert, x)
            if N.any(d < dist_cutoff_voxel):    continue

            pp_new_t.append(p)
        print 'in total', len(pp_new_t), 'peaks obtained'

        pp_new.extend(pp_new_t)
        


    with open('peak_segment_boundary_exclude__out.json', 'w') as f:     json.dump(pp_new, f, indent=2)





'''
in order for efficient calculation, we can also process this slice by slice. For each slice, we collect boundary points of segment, then remove all peaks that are within certain distance to these boundary points.
'''


