#!/usr/bin/env python


"""
plot isosurface of whole tomogram

~/ln/tomominer/tomominer/simulation/tomogram/analysis/tomogram_isosurface_plot.py

"""



import json
from mayavi import mlab
import numpy as N
import tomominer.io.file as IF
import tomominer.geometry.surface.gts.isosurface as GSGI

def main():
    with open('tomogram_isosurface_plot__op.json') as f:         op = json.load(f)

    if op['tomogram'].endswith('.npy'):
        v = N.load(op['tomogram'])
    elif op['tomogram'].endswith('.mrc') or op['tomogram'].endswith('.rec'):
        IF.read_mrc_vol(op['tomogram'])
    else:
        raise Exception('unknown file type')


    if not op['density positive']:      v = -v
    print 'min', v.min(), 'max', v.max()
    ie = GSGI.extract(v, isovalue=op['isovalue'], coarsen_n=(op['corsen number'] if 'corsen number' in op else None))
    print 'plotting...'

    figure = mlab.figure(bgcolor=tuple(op['bgcolor']))
    mlab.triangular_mesh(ie['vertices'][:,0], ie['vertices'][:,1], ie['vertices'][:,2], ie['faces'], color=tuple(op["color"]))

    if os.path.isfile(op['view']['file']):
        with open(op['view']['file']) as f:   view = json.load(f)
        mlab.view(azimuth=view['azimuth'], elevation=view['elevation'], distance=view['distance'], focalpoint=N.array(view['focalpoint']), roll=view['roll'])

    print 'saving', op['figure file out'], 'size', op['figure size'],
    mlab.savefig(op['figure file out'], size=op['figure size'])
    print 'done'
    mlab.show()   

if __name__ == '__main__':
    main()


