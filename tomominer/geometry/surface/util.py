
"""
functions for handling surfaces


~/ln/tomominer/tomominer/geometry/surface/util.py
"""



import copy
import numpy as N

import scipy.spatial.distance as SPD


# written according to Surface.concat()
def concat(s0=None, s1=None):

    if s0 is None:      return copy.deepcopy(s1)
    if s1 is None:      return copy.deepcopy(s0)

    n0 = s0['vertices'].shape[0]
    n1 = s1['vertices'].shape[1]

    s = {}
    s['vertices'] = N.vstack(   (s0['vertices'], s1['vertices'])    )
    s['faces'] = N.vstack(  (s0['faces'], s1['faces']+n0)   )

    # add a class label field so that in future the surfaces can be decomposed
    if 'class' not in s0:   s0['class'] = N.ones(n0)
    if 'class' not in s1:   s1['class'] = N.ones(n1)

    s['class'] = N.hstack(      (s0['class'], s1['class'] + s0['class'].max() + 1)      )

    return s


def concat_list(ss):
    s = None
    for st in ss:       s = concat(s, st)
    
    return s




'''
given the isosurface of one average, and the rigid transform of instances against this average, put the isosurface to the corresponding location in the map, the composed surface can be dumped as an mat file, which can be imported by Amira
options: ts: template isosurface, tc: the center of template.    dj[i]['angle'], dj[i]['loc']: rigid transforms     dj[i]['mid_co']: center of instance subtomograms
'''
import tomominer.geometry.ang_loc as GA
import tomominer.geometry.point_cloud.util as GPC
def instance_embedding_list(dj, ts, tc, reverse=True, min_dist=None, min_dist_instance_check_num=5):

    if tc is None:  tc = N.array(ts.shape, dtype=N.float) / 2.0

    mid_co_s = N.array([_['mid_co'] for _ in dj], dtype=N.float)

    count_face = 0
    count_vertice = 0

    ss = [None] * len(dj)
    for i, djt in enumerate(dj):
        if reverse:
            # in this case dj contains subtomogram alignment against template information
            rev_rm, rev_loc_r = GA.reverse_transform(GA.rotation_matrix_zyz(djt['angle']), djt['loc'])
            sr = {'faces':ts['faces'], 'vertices':GPC.rotate_translate(x0=ts['vertices'], c0=tc, c1=mid_co_s[i,:], rm=rev_rm, loc_r=rev_loc_r)}        # mxu 150427: bug fix: should use loc_r=rev_loc_r instead of loc_r=djt['loc']
        else:
            # in this case dj contains template transform for generating simulated tomogram
            sr = {'faces':ts['faces'], 'vertices':GPC.rotate_translate(x0=ts['vertices'], c0=tc, c1=mid_co_s[i,:], angle=N.array(djt['angle']), loc_r=(N.array(djt['loc']) if 'loc' in djt else None))}

        ss[i] = sr

        count_face += len(sr['faces'])
        count_vertice += len(sr['vertices'])

        print '\r', i, '   ', 'total vertices', count_vertice, 'total faces', count_face,

    return ss




'''
find and flag non-intersecting instances,
assume dj is ordered by importance
'''
def instance_find_nonintersect_instances(dj, ss, min_dist, min_dist_instance_check_num=5):


    mid_co_s = N.array([_['mid_co'] for _ in dj], dtype=N.float)

    ss_flag = N.zeros(len(ss), dtype=N.int)

    inds = N.array(range(len(ss)), dtype=N.int)
    for i, djt in enumerate(dj):

        if ss_flag.sum() == 0:
            ss_flag[i] = 1          # we want keep at least one non-intersect instance
            continue

        # this is just an adhoc way to exclude intersecting surfaces....

        # in order to save computation, find only the most close instances for checking
        md2 = SPD.cdist(mid_co_s[ss_flag==1, :].reshape((-1,3)), mid_co_s[i,:].reshape((1,3))).flatten()
        md2_si = N.argsort(md2)

        is_intersect = False
        min_dist_instance_check_num__count = 0
        for si_t in inds[ss_flag==1][md2_si]:
            assert si_t != i
            if ss[si_t] is None:        continue

            # then for each instance, check if its surface intesects with this new instance
            cd2 = SPD.cdist(ss[si_t]['vertices'], ss[i]['vertices'])
            if cd2.min() <= min_dist:
                is_intersect = True
                break

            min_dist_instance_check_num__count += 1
            if min_dist_instance_check_num__count > min_dist_instance_check_num:        break

        if is_intersect: continue

        ss_flag[i] = 1

        print '\r', i, '   ', 

    print '\r'

    return ss_flag




# use matlab to plot isosurface
def isosurface__matalb(v, isovalue, session=None, matlab_code_dir='/home/rcf-47/mxu/ln/frequent_structure/code'):

    if session is None:
        import pymatlab
        session = pymatlab.session_factory(options='-nodesktop -nodisplay')
        session.run( 'addpath(\'%s\')'%(matlab_code_dir,))

    session.run('clear all;')


    session.putvalue('v', v)
    session.run('isovalue = %f'%(isovalue,))

    session.run('fv = Surface.vol_isosurface(v, isovalue);')
    session.run('faces = fv.faces;')
    session.run('vertices = fv.vertices;')

    return {'faces':session.getvalue('faces'), 'vertices':session.getvalue('vertices')}



'''

import tomominer.model.util as MU
v = MU.generate_toy_model()

import tomominer.geometry.surface.util as GSU
fv = GSU.isosurface(v, v.mean())
GSU.dsp_surface(fv)

'''



'''
other ways to extract isosurface
https://github.com/pmneila/PyMCubes

https://github.com/Pyroevil/CubeSurfer

https://pyscience.wordpress.com/2014/09/11/surface-extraction-creating-a-mesh-from-pixel-data-using-python-and-vtk/

http://itk.org/ITKExamples/src/Core/Mesh/ExtractIsoSurface/Documentation.html
'''




def dsp_surface__matlab(fv, session=None, matlab_code_dir='/home/rcf-47/mxu/ln/frequent_structure/code', wait_enter=True):

    if session is None:
        import pymatlab
        session = pymatlab.session_factory(options='-nodesktop')
        session.run( 'addpath(\'%s\')'%(matlab_code_dir,))

    session.run('clear all;')

    session.putvalue('faces', fv['faces'])
    session.putvalue('vertices', fv['vertices'])

    session.run('fv = struct();')
    session.run('fv.faces = faces;')
    session.run('fv.vertices = vertices;')

    session.run('op = struct();')
    session.run('op.fv = fv;')

    session.run('figure;')
    session.run('Surface.dsp_surface(op);')
    session.run('drawnow;')

    if wait_enter:      raw_input('press enter')
    
    return


# use matlab to simplify isosurface in order to reduce number of data points
def resample__matlab(fv, keepratio=0.1, session=None, matlab_code_dir='/home/rcf-47/mxu/ln/frequent_structure/code'):

    if session is None:
        import pymatlab
        session = pymatlab.session_factory(options='-nodesktop -nodisplay')
        session.run( 'addpath(\'%s\')'%(matlab_code_dir,))

    session.run('clear all;')

    session.putvalue('faces', fv['faces'])
    session.putvalue('vertices', fv['vertices'])

    session.run('fv = struct();')
    session.run('fv.faces = faces');
    session.run('fv.vertices = vertices')

    session.run('keepratio = %f'%(keepratio,))

    session.run('fv = Surface.resample(fv, keepratio);')

    session.run('faces_t = fv.faces;')
    session.run('vertices_t = fv.vertices;')

    return {'faces':session.getvalue('faces_t'), 'vertices':session.getvalue('vertices_t')}



def dsp_surface_mayavi(fv, color=None, opacity=None):
    if color is not None:   color = tuple(color)

    from mayavi import mlab
    mlab.triangular_mesh(fv['vertices'][:,0], fv['vertices'][:,1], fv['vertices'][:,2], fv['faces'], color=color, opacity=opacity)



'''
test if a set of points are inside a polygon mesh, 

for hints, see
http://stackoverflow.com/questions/18135614/querying-of-a-point-is-within-a-mesh-maya-python-api
'''

def test_if_inside_mesh():
    raise Exception('To be implemented')

