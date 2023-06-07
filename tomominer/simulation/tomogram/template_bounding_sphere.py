#!/usr/bin/env python

# given a set of templates, segment and calculate bounding spheres

"""
For UCLA Cluster make sure matlab is loaded: module load matlab
"""


import sys, json, pickle, copy, pymatlab, gc
import numpy as N
import tomominer.io.file as IF
import tomominer.segmentation.active_contour.chan_vese.segment as SAC
import oct2py

def get_segment(m, op):

    m = N.abs(m)        # sometimes a map has negative intensities (especially at boundaries, from those FFT based pdb to density converters), remove them to reduce influence to segmentation
    phi = N.sign(m - m.mean())
    phi = SAC.segment(m, phi, smooth_weight=float(op['smooth_weight']), image_weight=float(op['image_weight'])/m.var(), delta_t=1.0, n_iters=1000)

    if m[phi > 0.0].mean() < m[phi < 0.0].mean():
        # sometimes the positive phi corresponds to lower intensity region, in such case, need to reverse phi
        print 'reverse phi'
        phi = -phi

    s = N.copy(m)
    s[phi <= 0.0] = float('NaN')       # use NaN for outside mask regions
    return s

def get_contour(m,pid, op):
    s = copy.deepcopy(m)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            for k in range(m.shape[2]):
                if m[i][j][k] > op['contour']:
                    s[i][j][k] = 1.0
                else:
                    s[i][j][k] = 0.0
    
    #out_o = "/tmp/origonal-" + pid + ".mrc"
    #out_s = "/tmp/segmented-" + pid + ".mrc"
    #IF.put_mrc(m, out_o)
    #IF.put_mrc(s, out_s)
    s[s==0.0] = float('NaN')
    return s

def get_bounding_sphere(s, session):
    #session.run('clear all;')
    session.push('s', s)
    session.eval('i = find(isfinite(s));')
    session.eval('[i1, i2, i3] = ind2sub(size(s), i);')
    session.eval('[center,radius] = minboundsphere( [i1, i2, i3] );')
    mb = {'c' : session.pull('center')-1, 'r' : session.pull('radius')}     # important: indicis in matlab starts from 1, but in python starts from 0!!
    return mb


if __name__ == '__main__':

    with open('template_bounding_sphere__config.json') as f:        op = json.load(f)
    with open(op['map_file']) as f:     maps = pickle.load(f)
    # -----------------------------------------------------
    print 'segmentating'
    segments = {}
    pid_id = 0
    matlab_dict = {}
    for pid in maps:
        segments[pid_id] = {}
        segments[pid_id]['map'] = get_contour(m=maps[pid]['map'], pid=pid,op=op['contour'])
        segments[pid_id]['pdb_id'] = pid
        print '\r', pid, '\tmean intensity inside:', maps[pid]['map'][N.isfinite(segments[pid_id]['map'])].mean(), 'outside:', maps[pid]['map'][N.isnan(segments[pid_id]['map'])].mean()
        pid_id += 1
        #segments[pid] = get_segment(m=maps[pid]['map'], op=op['segmentation'])
        #print '\r', pid, '\tmean intensity inside:', maps[pid]['map'][N.isfinite(segments[pid])].mean(), 'outside:', maps[pid]['map'][N.isnan(segments[pid])].mean()
        #IF.put_mrc(N.array(maps[pid]['map'], dtype=float, order='F'), '/tmp/%s.mrc'%(pid))
        #IF.put_mrc(N.array(N.isfinite(segments[pid]), dtype=float, order='F'), '/tmp/%s-seg.mrc'%(pid))t

    #import scipy.io
    #matlab_dict = {}
    #for pid in maps:
    #    matlab_dict[segments]
    #scipy.io.savemat('./arrays-segments.mat', mdict=segments)
    #session = pymatlab.session_factory(options='-nodisplay')
    session = oct2py.Oct2Py()
    session.eval( 'addpath(\'%s\')'%(op['matlab_code_path']) )
    session.eval( 'addpath(\'%s\')'%(op['matlab_bounding_code_path']) )

    # ------------------------------------------------------
    print 'calculating bounding sphere'
    bounding_spheres = {}
    for pid in segments:
        print '\r', pid,
	pdb_id = segments[pid]["pdb_id"]
        bounding_spheres[pdb_id] = get_bounding_sphere(s=segments[pid]['map'], session=session)
        bounding_spheres[pdb_id]['pdb_id'] = segments[pid]['pdb_id']
        bounding_spheres[pdb_id]['id'] = maps[segments[pid]['pdb_id']]['id']       # record a integer id for convinence
        print bounding_spheres[pdb_id]

    with open(op['output_file'], 'wb') as f:        pickle.dump({'segments':segments, 'bounding_spheres':bounding_spheres}, f)
