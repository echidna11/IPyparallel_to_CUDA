#!/usr/bin/env python



'''

given each, membrane bound particle, collect membrane surface points within membrane, then use PCA to obtain surface norm direction.
then calculate the angle difference of embeded template norm and filter




~/ln/tomominer/tomominer/template/search/filter/surface_norm_difference.py

'''


import json
import numpy as N

import scipy.spatial.distance as SPD

import tomominer.geometry.surface.vrml as GSV
import tomominer.geometry.surface.amira.util as GSAU
import tomominer.geometry.ang_loc as GA
import tomominer.geometry.point_cloud.util as GPU

import tomominer.io.file as IF



def surface_norm(x, vs, neighborhood_radius):

    d = SPD.cdist(vs, x).flatten()
    vst = vs[d <= neighborhood_radius, :]

    if vst.shape[0] < 3:    return None

    return GPU.surface_norm(vst)


def main():
    with open('surface_norm_difference__op.json') as f:     op = json.load(f)

    with open(op['data file in']) as f:      dj = json.load(f)

    mrc = IF.read_mrc_header(op['reference tomogram'])['MRC']

    nt = N.array(op['template norm'])
    nt = N.reshape(nt, (1, 3))

    with open(op['surface file']) as f:         vs = GSAU.surf_parse(f)['vertices']        # vertices on the surface

    x = N.zeros((len(dj), 3))
    for i in range(len(dj)):        x[i,:] = N.array(dj[i]['peak']['loc'])

    x_c = x.mean(axis=0)


    vs = GSAU.vertice_location_transform__amira_to_vol(x=vs, mrc=mrc)

    dj_new = []
    for i,d in enumerate(dj):
        n = surface_norm(x=N.reshape(x[i,:], (1,3)), vs=vs, neighborhood_radius=op['neighborhood radius'])
        if n is None:       continue
        
        if N.dot(x[i,:].flatten() - x_c, n) < 0:    n = -n          # make sure that the normal points outwards, in future, need to modify this using inside mesh test

        n0 = GPU.rotate_translate(x0=nt, rm=GA.rotation_matrix_zyz(d['angle']).T, c0=N.zeros(3)).flatten()

        ang_dif = GA.angle_between_vectors(n, n0, take_abs=op['angle take abs'])

        d['surface_norm'] = {'n0':n0.tolist(), 'n':n.tolist(), 'ang_dif':ang_dif}
        dj_new.append(d)

    with open(op['data file out'], 'w') as f:       json.dump(dj_new, f, indent=2)
    print 'number of data entries reduced from', len(dj), 'to', len(dj_new)
    

if __name__ == '__main__':
    main()



