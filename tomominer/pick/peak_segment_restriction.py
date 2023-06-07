#!/usr/bin/env python

# restrict the detected peaks to cell region segments defined using surface

import json, pickle

import numpy as N
import pymatlab


import tomominer.io.file as IF

import tomominer.geometry.surface.amira.util as GFAU




def process_surface_segmentation(session, segment_file, pp_tom, mrc):
    x = [_['x'] for _ in pp_tom]
    x = N.array(x, dtype=N.float)
   
    print 'loading', segment_file

    with open(segment_file) as f:  surf = GFAU.surf_parse(f)
    vert_org = surf['vertices']

    # transfer vertice location from original space to the voxel grid space
    vert = N.zeros(vert_org.shape)
    vert[:,0] = (   (vert_org[:,0] - mrc['xorg']) / mrc['xlen']   ) * mrc['nx']
    vert[:,1] = (   (vert_org[:,1] - mrc['yorg']) / mrc['ylen']   ) * mrc['ny']
    vert[:,2] = (   (vert_org[:,2] - mrc['zorg']) / mrc['zlen']   ) * mrc['nz']

    session.run('clear all;')
    session.putvalue('vert', vert)
    session.putvalue('faces', surf['faces'])
    session.putvalue('x', x)

    if False:
        session.run('is_in_t = PointCloud.inpolyhedron(faces, vert,  x);')
    else:
        session.run('is_in_t = PointCloud.InPolyedron(vert, faces, PointCloud.InPolyedron__get_surf_normal(vert, faces), x);')

    is_in_t = session.getvalue('is_in_t')
    print 'peaks found inside segment', is_in_t.sum()

    pp_tom = [_ for i,_ in enumerate(pp_tom) if (is_in_t[i])]

    return pp_tom


def process_vol_segmentation(segment_file, pp_tom, cutoff=0.5):
        
    print 'loading', segment_file

    m = IF.read_mrc_vol(segment_file)

    pp = []
    for p in pp_tom:
        x = N.round(N.array(p['x'])).astype(N.int)
        if m[x[0]][x[1]][x[2]] < cutoff:        continue
        pp.append(p)

    return pp



def main():

    with open('peak_segment_restriction__config.json') as f:     op = json.load(f)

    print 'loading', op['tomogram_info_file']
    with open(op['tomogram_info_file']) as f:    tom_info = json.load(f)

    print 'loading', op['peak_file']
    with open(op['peak_file']) as f:    pp = json.load(f)          # this file contains peaks with redundancy removed using particle_picking_dog__filter.py


    session = pymatlab.session_factory(options='-nodesktop -nodisplay')
    session.run( 'addpath(\'%s\')'%(op['matlab_code_path']) )

    pp_re = []
    for tom in tom_info:
        print 'tomogram', tom['id']

        print 'loading', tom['header_file']
        with open(tom['header_file'], 'rb') as f:   mrc = pickle.load(f)['header']['MRC']
        

        pp_tom = [_ for _ in pp if (_['tomogram_id'] == tom['id'])]
        if len(pp_tom) == 0:
            print 'no peaks'
            continue

        pp_tom = sorted(pp_tom, key=lambda _ : _['id'])         # order peaks according to ids,  assume the peak ids are ordered so that the smaller id has better value / importance

        if tom['segment_file'].lower().endswith('.surf'):
            pp_tom = process_surface_segmentation(session=session, segment_file=tom['segment_file'], pp_tom=pp_tom, mrc=mrc)
        elif tom['segment_file'].lower().endswith('.mrc') or tom['segment_file'].lower().endswith('.rec'):
            pp_tom = process_vol_segmentation(segment_file=tom['segment_file'], pp_tom=pp_tom)
        else:
            raise Exception('Segmentation file type unrecognized')

        if 'top_num' in op:        pp_tom = pp_tom[:op['top_num']]          # assume peaks are ordered by response

        pp_re.extend(pp_tom)

    print 'filtered from', len(pp), 'to', len(pp_re), 'peaks'
    with open(op['out_file'], 'w') as f:       json.dump(pp_re, f, indent=2)


if __name__ == '__main__':
    main()

