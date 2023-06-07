#!/usr/bin/env python

# given a volume and a set of markers with radius, display their locations


import sys
import numpy as N

# parameters: v: volume    x: locations      r:radius
current_slice = None
def display(v, x=None, r=None, linewidth=2.0):
    print 'display....... press left or right arrow to change slices'

    import matplotlib.pyplot as PLT
    fig = PLT.figure()

    global current_slice
    current_slice = int(v.shape[2]/2)

    ax = fig.add_subplot(111)
    ax.old_axis = None

    import matplotlib.cm as cm
    def display_slice():
        global current_slice
        print '\r', current_slice, '       ',       ;        sys.stdout.flush()

        old_axis = ax.axis()

        PLT.cla()
        ax.imshow(v[:,:,current_slice].transpose(), cmap = cm.Greys_r)            # IMPORTANT!! Transpose is needed to make the x-y coordinates of images consistant with coordinateds in x!!!

        if x is not None:
            for xt in x:
                if N.abs(xt[2] - current_slice) >= r:   continue

                rad = N.sqrt( r*r - N.square(current_slice - xt[2]) )

                if xt[2] <= current_slice:
                    ax.plot([xt[0]-rad, xt[0]+rad], [xt[1], xt[1]], 'g-', linewidth=linewidth)

                if xt[2] >= current_slice:
                    ax.plot([xt[0], xt[0]], [xt[1]-rad, xt[1]+rad], 'r-', linewidth=linewidth)

        if ax.old_axis is not None:
            ax.axis(old_axis)
        ax.old_axis = old_axis

        fig.canvas.draw()


    def on_key_press(event):
        global current_slice
        if (event.key == 'left') and (current_slice > 0):
            current_slice -= 1
            display_slice()

        elif (event.key == 'right') and (current_slice < (v.shape[2]-1)):
            current_slice += 1
            display_slice()


    fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)
    cid = fig.canvas.mpl_connect('key_press_event', on_key_press)

    display_slice()

    PLT.show()

    fig.canvas.mpl_disconnect(cid)


'''
# test code

import matplotlib
matplotlib.use('Qt4Agg')

import numpy as N
v = N.random.random((50,50,50))

sigma = 3
import tomominer.filter.gaussian as FG
v = FG.smooth(v, sigma)

import tomominer.filter.local_extrema as FL

if False:
    x = FL.local_maxima(v)
else:
    x = FL.local_minima(v)

x = N.array(x).T

import tomominer.plot.marker_display as PM
PM.display(v=v, x=x, r=sigma)


'''

