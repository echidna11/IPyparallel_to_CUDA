#!/usr/bin/env python



'''


exclude those peaks that are close to the boundary of segmentation. This is just an adhoc way to remove those peaks that corresponds to membrane.





modified according to ~/ln/tomominer/tomominer/pick/peak_segment_boundary_exclude.py


~/ln/tomominer/tomominer/pick/subtomogram/filtering/location/segment_boundary_exclude.py

'''




if __name__ == '__main__':

    import json
    with open('segment_boundary_exclude__op.json') as f:     op = json.load(f)

    dist_cutoff_voxel = float(op['dist cutoff voxel'])

    with open(op['data file in']) as f:    dj = json.load(f)

    import pickle
    import numpy as N
    from scipy.spatial.distance import cdist
 
    import tomominer.io.file as IF

    mrc = IF.read_mrc_header(op['tomogram file'])['MRC']

    import tomominer.geometry.surface.amira.amira as GFAA
    with open(op['segment file']) as f:  surf = GFAA.surf_parse(f)
    vert_org = surf['vertices']

    # transfer vertice location from original space to the voxel grid space
    vert = N.zeros(vert_org.shape)
    vert[:,0] = (   (vert_org[:,0] - mrc['xorg']) / mrc['xlen']   ) * mrc['nx']
    vert[:,1] = (   (vert_org[:,1] - mrc['yorg']) / mrc['ylen']   ) * mrc['ny']
    vert[:,2] = (   (vert_org[:,2] - mrc['zorg']) / mrc['zlen']   ) * mrc['nz']

    dj_new = []
    for i, p in enumerate(dj):
        print '\r', float(i)/len(dj), '          ',

        x = N.array(p['peak']['loc'], dtype=N.float).reshape(     (1,3)      )
        d = cdist(vert, x)
        if N.any(d < dist_cutoff_voxel):    continue

        dj_new.append(p)

    print 'in total', len(dj_new), 'peaks obtained'

    with open(op['data file out'], 'w') as f:     json.dump(dj_new, f, indent=2)





