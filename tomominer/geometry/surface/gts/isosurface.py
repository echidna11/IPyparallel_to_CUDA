
"""
functions for handling surfaces and isosurfaces through PyGTS
"""


'''
~/ln/tomominer/tomominer/geometry/surface/gts/util.py
'''



import gts
import numpy as N


"""
extract isosurface into a triangle mesh,
simplify if needed using coarse_n
"""
def extract(v, isovalue, coarsen_n=None):

    s2 = N.array(v.shape) / 2.0

    iso = gts.isosurface(v, isovalue)

    if coarsen_n is not None:        iso.coarsen(coarsen_n)

    x0, x1, x2, t = gts.get_coords_and_face_indices(iso, unzip=True)
    x0 = N.array(x0).flatten() * s2[0] + s2[0] + 0.5
    x1 = N.array(x1).flatten() * s2[1] + s2[1] + 0.5
    x2 = N.array(x2).flatten() * s2[1] + s2[2] + 0.5

    return {'vertices': N.array([x0, x1, x2]).T, 'faces': N.array(t)}


    
    



'''

#useful functions

S = gts.Surface()
S.add(s0)
S.add(s1)

x,t = gts.get_coords_and_face_indices(surf)

'''



'''

# example code to use gts to get isosurface, and use gts to simplify


import tomominer.model.util as MU
v = MU.generate_toy_model()

def plot_surf(surf):
    x,y,z,t = gts.get_coords_and_face_indices(surf,True)
    from mayavi import mlab
    mlab.clf()
    mlab.triangular_mesh(x,y,z,t,color=(0.25,0.25,0.75))
    mlab.show()


import gts
s = gts.isosurface(v, 0.1)     # see http://sourceforge.net/p/pygts/code/HEAD/tree/examples/isosurface.py
s.Nvertices
s.Nfaces
plot_surf(s)


s1 = gts.Surface()
s1.copy(s)          # need to make a copy before simplify
s1.coarsen(100)     # simplify to under 100 vertices / faces
s1.Nvertices
s1.Nfaces
plot_surf(s1)


s1 = gts.Surface()
s1.copy(s)      # need to make a copy before simplify
s1.coarsen(1000)    # simplify to under 1000 vertices / faces
s1.Nvertices
s1.Nfaces
plot_surf(s1)


'''



# check how to scale and translate the gts isosurface coordinates to fit the coordinate of original volumeric data in mayavi
def scale_location_check():

    s = 13
    s2 = s/2
    s2s = s/2.0

    import numpy as N

    v = N.zeros([s,s,s])

    v[(s2-1):(s2+1), (s2-2):(s2+2), (s2-3):(s2+3)] = 1.0

    from tomominer.geometry.surface.gts.isosurface import extract
    iso = extract(v, 0.5, coarsen_n=10)


    from mayavi import mlab
    mlab.clf()


    mlab.triangular_mesh(iso['vertices'][:,0], iso['vertices'][:,1], iso['vertices'][:,2], iso['faces'], color=(0.25,0.25,0.75), opacity=0.5)


    from tomominer.image.vol.util import dsp_orthogonal_slices__matavi
    dsp_orthogonal_slices__matavi(v, c=(s2s, s2s, s2s))

    mlab.show()


    
    # export to amira format, see if the resulting surf file can be read by amira
    import copy
    iso_t = copy.deepcopy(iso)
    iso_t['faces'] += 1
    import tomominer.geometry.surface.amira.util as GSAA
    with open('/tmp/tmp.surf', 'w') as f:    GSAA.export_surf_ascii_simple(iso_t, f)

