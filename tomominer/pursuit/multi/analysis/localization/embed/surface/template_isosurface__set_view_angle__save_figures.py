#!/usr/bin/env python




'''
given the output of,        template_isosurface__set_view_angle.py
save corresponding isosurface figures


~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/template_isosurface__set_view_angle__save_figures.py

'''

import os, json, pickle
import numpy as N

from mayavi import mlab

def main():
    with open('template_isosurface__set_view_angle__save_figures__op.json') as f:     op = json.load(f)
    with open(op['set_view_angle__op file']) as f:       vop = json.load(f)

    with open(vop['view file out']) as f:     dj = json.load(f)

    with open(vop['mesh file in'], 'rb') as f:       surf = pickle.load(f)

    if 'color file' in vop:
        with open(vop['color file']) as f:       colors = json.load(f)
        colors = [tuple(colors[d['subtomogram']]['color']) for d in dj]
    else:
        colors = None


    op['out dir'] = os.path.abspath(op['out dir'])
    if not os.path.isdir(op['out dir']):        os.makedirs(op['out dir'])

    for i, d in enumerate(dj):
        fv = surf[i]['fv']
        fig = mlab.figure(bgcolor=tuple(op['bgcolor']), size=tuple(op['figure size']))
        mlab.triangular_mesh(fv['vertices'][:,0], fv['vertices'][:,1], fv['vertices'][:,2], fv['faces'], color=(colors[i] if colors is not None else tuple(op['default color'])), opacity=vop['opacity'])
        mlab.view(azimuth=d['view']['azimuth'], elevation=d['view']['elevation'], distance=d['view']['distance'], focalpoint=N.array(d['view']['focalpoint']), roll=d['view']['roll'])
        d['isosurface'] = {'file': os.path.join(op['out dir'], '%03d.png'%(i,))}
        mlab.savefig(filename=d['isosurface']['file'])
        mlab.close()


    with open(op['stat out'], 'w') as f:     dj = json.dump(dj, f, indent=2)






if __name__ == '__main__':
    main()



'''

related code


~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/template_isosurface__set_view_angle.py

~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/template_isosurface.py

'''

