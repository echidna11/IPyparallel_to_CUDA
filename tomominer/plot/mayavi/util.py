'''
functions for mayavi plotting

~/n/tomominer/tomominer/plot/mayavi/util.py
'''

from mayavi import mlab
from tvtk.tools import visual
import numpy as N


'''
# example use of visual

from mayavi import mlab
from tvtk.tools import visual
fig=mlab.figure( bgcolor=(0, 0, 0))
visual.set_viewer(fig)
'''

'''
adapted from
http://stackoverflow.com/questions/20085398/simple-arrow-mayavi-tvtk-in-mlab-weird-behaviour-looks-like-a-bug
'''
def arrow_from_a_to_b(x1, y1, z1, x2, y2, z2):
    ar1=visual.arrow(x=x1, y=y1, z=z1)
    ar1.length_cone=0.4

    arrow_length=N.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    ar1.actor.scale=[arrow_length, arrow_length, arrow_length]
    ar1.pos = ar1.pos/arrow_length
    ar1.axis = [x2-x1, y2-y1, z2-z1]
    return ar1

