#!/usr/bin/env python


'''

this program inverse rotate the surface segmenation (obtained from a rotated tomogram) back to the original tomogram

pipeline:

given a tomogram v
bin v to make a smaller tomogram vb
rotate vb to vbr using the slicer program     ~/ln/tomominer/tomominer//image/vol/plot/slicer_mayavi.py
manually segment vbr to obtain a segment surface file vbr.surf

use this program to inverse rotate vbr.surf to obtain the surface vbri.surf that correspond to vb or v



~/ln/tomominer/tomominer/segmentation/manual/surface/surf_inverse_rotate.py


'''


def main():
    import json

    with open('surf_inverse_rotate.json') as f:     op = json.load(f)

    with open(op['rot_param_file']) as f:     rot = json.load(f)

    import numpy as N
    loc = N.array(rot['loc'])
    ang = N.array(rot['angle'])


    import tomominer.geometry.ang_loc as GA
    (ang_rev, loc_rev) = GA.reverse_transform_ang_loc(ang, loc)             # calculate inverse rigid transform


    import tomominer.io.file as IF
    mrc_org = IF.read_mrc_header(op['map_org_file'])['MRC']
    siz = N.array([mrc_org['nx'], mrc_org['ny'], mrc_org['nz']])

    mrc_rot = IF.read_mrc_header(op['map_rot_file'])['MRC']

    import tomominer.geometry.surface.amira.util as GFAU
    with open(op['surf_file_in']) as f:  surf = GFAU.surf_parse(f)
    vert_org = surf['vertices']         
    vert_map = GFAU.vertice_location_transform__amira_to_vol(vert_org, mrc_rot)             # transform vertice coordinates from volume space to map grid space


    import tomominer.geometry.point_cloud.util as GPU
    vert_map_inv = GPU.rotate_translate(vert_map, angle=ang_rev, loc_r=loc_rev, c0=N.array(siz)/2.0)        # rotate the vertices according to inverse transform

    vert_inv = GFAU.vertice_location_transform(vert_map_inv, mrc_org)       # transform vertice coordinates from map grid space to volume space


    import copy
    surf_inv = copy.deepcopy(surf)
    surf_inv['vertices'] = vert_inv

    with open(op['surf_file_out'], 'w') as f:       GFAU.export_surf_ascii_simple(surf_inv, f)

if __name__ == '__main__':
    main()

