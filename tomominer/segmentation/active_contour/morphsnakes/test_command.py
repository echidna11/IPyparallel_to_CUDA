


# Test on a volume containing a sphere, with noise added.

dim_len = 200
shape = [dim_len, dim_len, dim_len]

import tomominer.model.util as MU
#v = MU.sphere_mask(shape, radius=dim_len/4)

import numpy as N
v = N.zeros(shape)
v[(dim_len/2 - dim_len/4) : (dim_len/2 + dim_len/4), (dim_len/2 - dim_len/4) : (dim_len/2 + dim_len/4), (dim_len/2 - dim_len/4) : (dim_len/2 + dim_len/4)] = 1      # rectangle


import tomominer.image.vol.util as GV

import matplotlib
matplotlib.use('Qt4Agg')
from matplotlib import pyplot as ppl

GV.dspcub(v)


# run one iteration, and plot the segment

draw = True

import sys
import tomominer.segmentation.active_contour.morphsnakes.morphsnakes as SM

vg = SM.gborders(v, alpha=1000, sigma=2)
if draw:    GV.dspcub(vg)

# Morphological GAC. Initialization of the level-set.
mgac = SM.MorphGAC(vg, smoothing=1, threshold=0.31, balloon=1)
mgac.levelset = MU.sphere_mask(shape, radius=dim_len/8)

# Visual evolution.

if draw:
    fig = ppl.gcf()
    fig.clf()
    ax = fig.add_subplot(1,1,1)
    ax_u = ax.imshow(GV.cub_img(mgac.levelset)['im'], cmap=ppl.cm.gray)
    ppl.pause(0.01)


import sys
import time
for i in range(100):
    start_time = time.time()
    mgac.step()
    print '\r', i, '    ', time.time() - start_time, '           ',
    sys.stdout.flush()
    if draw:
        ax_u.set_data( GV.cub_img(mgac._u)['im'] )
        fig.canvas.draw()
        ppl.pause(0.01)

