#!/usr/bin/env python

'''
plot isosurfaces of instances embeded into density map, using different colors prepared using 

~/ln/tomominer/tomominer/simulation/tomogram/analysis/class_isosurface_plot.py

'''

'''
todo: consider to merge the plotting part with 
~/ln/tomominer/tomominer/simulation/tomogram/analysis/class_isosurface_plot.py

'''


import os, json, pickle, copy
import numpy as N
from mayavi import mlab
import gts
import tomominer.geometry.surface.util as GSU
import tomominer.geometry.surface.gts.isosurface as GSGI

def main():
    with open('class_isosurface_plot__op.json') as f:         op = json.load(f)

    with open(op['map file'], 'rb') as f:       maps = pickle.load(f)
    with open(op['model op file']) as f:       model_op = json.load(f)
    with open(op['model file']) as f:       models = json.load(f)
    with open(op['bounding sphere file'], 'rb') as f:       bounding_spheres = pickle.load(f)
    with open(op['color file in']) as f:        colors = json.load(f)
    bounding_spheres = bounding_spheres['bounding_spheres']

    box  = model_op['packing']['param']['box']
    map_size = N.array( [box['x'], box['y'], box['z']] )

    figure = mlab.figure(bgcolor=tuple(op['bgcolor']), size=(1000,1000))

    def picker_callback(picker):
        e = mlab.view()
        r = mlab.roll()
        view = {'azimuth':e[0], 'elevation':e[1], 'roll':r, 'distance':e[2], 'focalpoint':e[3].tolist()}
        print 'view', view

        if op['view']['record']:
            with open(op['view']['file'], 'w') as f:       json.dump(view, f, indent=2)
            print 'saved to', op['view']['file']


    picker = figure.on_mouse_pick(picker_callback)

    model_id = op['model id']
    for pid in maps:
        if pid in op["reject_pid"]:
            continue
        print 'generating instances for', pid
        model_t = [_ for _ in models[str(model_id)]['instances'] if (_['pdb_id'] == pid)]
        if len(model_t) == 0:       continue

        model_t = copy.deepcopy(model_t)
        for d in model_t:            d['mid_co'] = d['x']

        # generate and simplify isosurface
        v = maps[pid]['map']
        v_iso = GSGI.extract(v, isovalue=op['contour'], coarsen_n=op['corsen number'])

        #ss = GSU.instance_embedding_list(dj=model_t, ts=v_iso, tc=N.array(v.shape)/2.0, reverse=False)            # according to simulation.tomogram.model_to_density_map.get_map(), tc=N.array(v.shape)/2.0, not tc=bounding_spheres[pid]['c']
        ss = GSU.instance_embedding_list(dj=model_t, ts=v_iso, tc=bounding_spheres[pid]['c'], reverse=False)
        ie = GSU.concat_list(ss)

        # generate and add the corresponding isosurface
        if 'opacity' in op:
            mlab.triangular_mesh(ie['vertices'][:,0], ie['vertices'][:,1], ie['vertices'][:,2], ie['faces'], color=tuple(colors[pid]), opacity=op['opacity'])
        else:
            mlab.triangular_mesh(ie['vertices'][:,0], ie['vertices'][:,1], ie['vertices'][:,2], ie['faces'], color=tuple(colors[pid]))


    if os.path.isfile(op['view']['file']):
        with open(op['view']['file']) as f:   view = json.load(f)
        mlab.view(azimuth=view['azimuth'], elevation=view['elevation'], distance=view['distance'], focalpoint=N.array(view['focalpoint']), roll=view['roll'])

    if 'figure file out' in op:
        print 'saving', op['figure file out'], 'size', op['figure size']
        mlab.savefig(op['figure file out'], size=op['figure size'])         # notice: you can save as VRML format, and open it using chimera, or amira!!
        print 'done'

    mlab.show()


if __name__ == '__main__':
    main()



