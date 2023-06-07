#!/usr/bin/env python


"""

plot embeded instances (generated through instance_embed.py) using mayavi, save instances for each pattern seperately.

Save to vrml files directly so that they can be opened and explored using Chimera, do not plot using mavavi


~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/instance_embed_seperate_plot_mayavi.py

"""




import os, json, pickle

from mayavi import mlab
import numpy as N


def main():
    with open('instance_embed_seperate_plot_mayavi__op.json') as f:     op = json.load(f)

    op['out dir'] = os.path.abspath(op['out dir'])
    if not os.path.isdir(op['out dir']):        os.makedirs(op['out dir'])

    print 'loading', op['mesh file']
    with open(op['mesh file'], 'rb') as f:       s = pickle.load(f)

    if ('color file' in op) and os.path.isfile(op['color file']):
        print 'loading', op['color file']
        with open(op['color file']) as f:        colors = json.load(f)
    else:
        colors = None

    figure = mlab.figure(bgcolor=tuple(op['bgcolor']), size=(1000,1000))

    for i in s:
        ie = s[i]['fv']
        tm = mlab.triangular_mesh(       ie['vertices'][:,0], ie['vertices'][:,1], ie['vertices'][:,2], ie['faces'], color=(     tuple(colors[s[i]['subtomogram']]['color']) if ((colors is not None) and (s[i]['subtomogram'] in colors)) else None    )        )

        out_file = os.path.join(op['out dir'], '%04d.vrml'%(i,))
        print 'saving', out_file
        mlab.savefig(out_file)

        tm.remove()

    #mlab.show()



if __name__ == '__main__':
    main()


'''
code similiar to
~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/instance_embed_plot_mayavi.py
~/ln/tomominer/tomominer/simulation/tomogram/analysis/class_isosurface_plot.py
'''


