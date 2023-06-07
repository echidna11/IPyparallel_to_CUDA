#!/usr/bin/env python

'''
plot embeded instances (generated through instance_embed.py) using mayavi
~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/instance_embed_plot_mayavi.py

code similiar to
~/ln/tomominer/tomominer/simulation/tomogram/analysis/class_isosurface_plot.py
'''

import os, json, pickle
from mayavi import mlab
import numpy as N


def main():
    with open('instance_embed_plot_mayavi__op.json') as f:     op = json.load(f)

    print 'loading', op['mesh file']
    with open(op['mesh file'], 'rb') as f:       s = pickle.load(f)

    if ('color file' in op) and os.path.isfile(op['color file']):
        print 'loading', op['color file']
        with open(op['color file']) as f:        colors = json.load(f)
    else:
        colors = None

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

    for i in s:
        ie = s[i]['fv']
        if 'opacity' in op:
            mlab.triangular_mesh(   ie['vertices'][:,0], ie['vertices'][:,1], ie['vertices'][:,2], ie['faces'], color= (tuple(colors[s[i]['subtomogram']]) if ((colors is not None) and (s[i]['subtomogram'] in colors)) else None), opacity=(op['opacity'] if 'opacity' in op else None)        )
        else:
            mlab.triangular_mesh(       ie['vertices'][:,0], ie['vertices'][:,1], ie['vertices'][:,2], ie['faces'], color=(     tuple(colors[s[i]['subtomogram']]) if ((colors is not None) and (s[i]['subtomogram'] in colors)) else None    )        )



    if ('view' in op) and ('file' in op['view']) and os.path.isfile(op['view']['file']):
        print 'loading view file', op['view']['file']
        with open(op['view']['file']) as f:   view = json.load(f)
        mlab.view(azimuth=view['azimuth'], elevation=view['elevation'], distance=view['distance'], focalpoint=N.array(view['focalpoint']), roll=view['roll'])

    if 'figure file out' in op:
        print 'saving', op['figure file out'], 'size', op['figure size'],
        mlab.savefig(op['figure file out'], size=op['figure size'])
        print 'done'

    mlab.show()



if __name__ == '__main__':
    main()

