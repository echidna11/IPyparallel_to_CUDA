


'''

functions for handling vrml



~/ln/tomominer/tomominer/geometry/surface/vrml.py

'''






def amira_surface_parse(vrml_file):
    from OpenGLContext.loaders.loader import Loader
    scenegraph = Loader.load( vrml_file, baseURL=None )
    return {'vertices':scenegraph.getDEF(scenegraph.DEF).vector}




'''

info of vrml97
http://pyopengl.sourceforge.net/context/vrml97.html


pip install --pre OpenGLContext-full



# example usage

vrml_file = '/home/rcf-47/mxu/ln/frequent_structure/data/out/analysis/jensen/ding/150503/match/segmentation/surface/Bdellovibrio_1k-labels_om_mask.wrl'

from OpenGLContext.loaders.loader import Loader
scenegraph = Loader.load( vrml_file, baseURL=None )

# to get points
scenegraph.getDEF(scenegraph.DEF).vector

# to get triangles, use scipy.spatial.Delaunay()


'''






'''
PyVRML97


to install pyvrml97

pip install --pre pyvrml97
pip install simpleparse PyDispatcher



# test code for using PyVRML97

vrml_file = '/home/rcf-47/mxu/ln/frequent_structure/data/out/analysis/jensen/ding/150503/match/segmentation/surface/Bdellovibrio_1k-labels_om_mask.wrl'
with open(vrml_file) as f:   v = f.read()

p = VVP.buildParser()
o = p.parse(v)
o[1][1]


'''



