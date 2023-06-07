#!/usr/bin/env python


'''
use simple plotting to determine the color, contour level, coarsen level of a batch of templates / averages. Then save the triangle mesh data.
~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/template_isosurface.py

note: initial list of averages can be obtained using some commands like following:
grep common_frame ../../out/pass/025/data.json | sort | uniq | grep clus_vol_avg_ > averages.json
'''




import json, pickle, copy
import numpy as N
from mayavi import mlab
import tomominer.image.vol.util as IVU
import tomominer.geometry.surface.gts.isosurface as GSGI
import tomominer.geometry.surface.util as GSU
import tomominer.filter.gaussian as FG
import tomominer.segmentation.util as SU
import tomominer.io.file as IF
def extract_isosurface(v0, op):

    v0 = N.copy(v0)

    if not op['density positive']:      v0 = -v0
    if 'shift mean' in op and op['shift mean']:        v0 -= v0.mean()

    if 'smooth sigma' in op:    v0 = FG.smooth(v0, sigma=op['smooth sigma'])

    print 'v0', 'min', v0.min(), 'max', v0.max(), 'mean', v0.mean()
    assert op['iso value'] > v0.mean()
    assert op['iso value'] < v0.max()

    # sometimes at boundary of volume the isosurface is not closed. In order to make the isosurface closed, we first enlarge the volume to double size
    v = IVU.resize_center(v0, N.array(v0.shape)*2.0, cval=v0.mean())

    s0 = v > op['iso value']         ;           assert  s0.sum() > 0
    s = copy.deepcopy(s0)

    # -------------------------------------
    # calculate size of connected components
    sc = SU.connected_components(s)

    l_n_s = N.zeros(sc['num']+1)
    scl = sc['label']
    for l in range(1, sc['num']+1):     l_n_s[l] = (sc['label'] == l).sum()
    print 'connected component sizes', sorted(l_n_s, reverse=True)


    chosen_labels = [True] * len(l_n_s)
    if ('small connected component size cutoff' in op) and (op['small connected component size cutoff'] > 0):
        print '# remove particles size smaller than', op['small connected component size cutoff']
        for l in range(len(l_n_s)):
            if l_n_s[l] < op['small connected component size cutoff']:     chosen_labels[l] = False
        
    if ('largest connected component only' in op) and (op['largest connected component only'] > 0) and (op['largest connected component only'] + 1 < len(l_n_s)):
        print '# only keep', op['largest connected component only'],'largest connected component(s)'
        l_n_s_i = N.argsort(-l_n_s)
        for i in range(op['largest connected component only']+1, len(l_n_s_i)):     chosen_labels[l_n_s_i[i]] = False

    s = N.zeros(sc['label'].shape, dtype=bool)
    for l in range(len(l_n_s)):
        if chosen_labels[l]:    s[sc['label'] == l] = True

    v[N.logical_and(s0, N.logical_not(s))] = v0.min()       # mask out values at rejected segment
    assert  (v >= op['iso value']).sum() > 0
    IF.put_mrc(v, '/tmp/v.mrc')
    IF.put_mrc(v >= op['iso value'], '/tmp/v-s.mrc')

    fv = GSGI.extract(v=v, isovalue=op['iso value'], coarsen_n=(op['coarsen_n'] if 'coarsen_n' in op else None))

    # move vertices back to scale of v0
    x = fv['vertices']
    for i in range(3):      x[:,i] -= v0.shape[i] / 2.0

    return fv


def main():
    with open('template_isosurface__op.json') as f:     op = json.load(f)
    with open(op['template file']) as f:     dj = json.load(f)
    print 'loaded', len(dj), 'averages'

    if 'color file' in op:
        with open(op['color file']) as f:       colors = json.load(f)
    else:
        colors = None

    surf = {}
    for j in range(len(dj)):
        d = dj[j]['subtomogram']
        i = dj[j]['id']
        print d, i

        v = IF.read_mrc_vol(d)
        print d, 'min', v.min(), 'max', v.max()

        surf[i] = {}
        surf[i]['subtomogram'] = d
        surf[i]['fv'] = extract_isosurface(v, op['extract'])
        surf[i]['center'] = N.array(v.shape) / 2.0

        mlab.figure(figure=str(i))
        GSU.dsp_surface_mayavi(fv=surf[i]['fv'], color=colors[d], opacity=op['opacity'])
        if op['plot slices']:    IVU.dsp_orthogonal_slices__matavi(v, vmin=-op['extract']['iso value'], vmax=op['extract']['iso value'])

    print 'saving trangular mesh to', op['mesh file out']
    with open(op['mesh file out'], 'wb') as f:          pickle.dump(surf, f, protocol=-1)
    print 'done'

    mlab.show()


if __name__ == '__main__':
    main()
