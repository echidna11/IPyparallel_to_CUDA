# find ridge regions of a given volume

import numpy as N
import scipy.ndimage.interpolation as SNI

import tomominer.filter.differential as FD
import tomominer.linalg.eigen as LE


def feature(vg, alpha=1.0):

    print '# calculate hessian matrices'
    dif = FD.diff_3d(vg)
    h = FD.hessian_3d(vg, d=dif)

    hmm = FD.hessian_3d__max_magnitude(h)
    h = FD.hessian_3d__normalize(h, hmm)

    print '# calculate eignevalues and largest eigen vectors of hessian matrices'
    e = LE.eigen_value_3_symmetric_batch(h)
    c = LE.eigen_vector_3_symmetric_batch__given_eigen_value(h, e[0])

    c = N.array(c)
    c[N.logical_not(N.isfinite(c))] = 0.0       # change invalid eigen vectors to zeros

    print '# calculate find ridge region'       #       according to first part of equation 6 of paper      Martinez-Sanchez13 A ridge-based framework for segmentation of 3D electron microscopy datasets
    flag = N.zeros(vg.shape, dtype=N.bool)
    flag[:] = True

    g = N.mgrid[0:vg.shape[0], 0:vg.shape[1], 0:vg.shape[2]]
    neighbor_compare(flag=flag, vg=vg, g=g, alpha=alpha, c=c)
    neighbor_compare(flag=flag, vg=vg, g=g, alpha=(-alpha), c=c)

    return flag


def neighbor_compare(flag, vg, g, alpha, c):
    gd = g + alpha*c

    vgi = SNI.map_coordinates(vg, gd)

    flag[vg < vgi] = False      # this indicates the corresponding position are not local maxima along the direction of principle eigen vector


