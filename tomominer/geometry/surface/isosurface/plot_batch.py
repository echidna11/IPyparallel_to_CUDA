#!/usr/bin/env python



'''

given a list of cluster averages or templates, for each one, plot an isosurface, use GUI to control view density cutoff and view angle, then record these information

~/ln/tomominer/tomominer/geometry/surface/isosurface/plot_batch.py

'''


'''
related code, see
https://sites.google.com/site/xumin100/home/technical-topics/graphics/mayavi

'''



import os, json, copy

import tomominer.io.file as IF

def main():
    with open('geometry_surface_isosurface_plot_batch__op.json') as f:        op = json.load(f)

    if os.path.isfile(op['data file out']):
        # when we already have a record of previous results
        print 'using existing', op['data file out']
        with open(op['data file out']) as f:    dj = json.load(f)
    else:
        print 'loading', op['data file in']
        with open(op['data file in']) as f:   dj = json.load(f)

    # load averages / templates
    vs = {}
    for i, d in enumerate(dj):
        vt = IF.read_mrc_vol(d['subtomogram'])
        vs[i] = vt if op['intensity positive'] else (-vt)


    my_models = MyModels(dj, vs, out_file=op['data file out'])
    my_models.configure_traits()


import numpy as N

import traits.api as TA
import traitsui.api as TUA
import mayavi.modules.api as TMA
import mayavi.core.api as MCA
import mayavi.core.ui.api as MCUA

class MyModels(TA.HasTraits):
    scene = TA.Instance(MCUA.MlabSceneModel, ())

    plot = TA.Instance(MCA.PipelineBase)

    def __init__(self, dj, vs, out_file):
        TA.HasTraits.__init__(self)

        self.dj = copy.deepcopy(dj)
        self.vs = vs
        self.out_file = out_file

        v_max_min = []
        for i in vs:    v_max_min.extend([vs[i].max(), vs[i].min()])

        if 'view' in self.dj[0]:
            int_cutoff_default = self.dj[0]['view']['contour']
        else:
            int_cutoff_default = max(v_max_min)/10

        self.add_trait('int_cutoff', TA.Range(low=float(N.floor(min(v_max_min))), high=float(N.ceil(max(v_max_min))), value=int_cutoff_default))  # this trait's range need to be determined at initialization stage, it is not a statistic trait

        self.add_trait('vol_ind', TA.Range(low=0, high=len(vs)-1, value=0))

        self.add_trait('save', TA.Button())

        self.current_vol_ind = -1
        self.iso = None


    @TA.on_trait_change('int_cutoff, vol_ind')
    def update_plot(self):
        if self.iso is not None:    self.record_view()

        if self.current_vol_ind != int(self.vol_ind):
            self.current_vol_ind = int(self.vol_ind)
            if self.iso is not None:
                self.current_dist = self.scene.mlab.view()[2]
                self.iso.remove()
            elif 'view' in self.dj[self.current_vol_ind]:
                self.current_dist = self.dj[self.current_vol_ind]['view']['distance']
            else:
                self.current_dist = None
            print self.dj[self.current_vol_ind]['subtomogram']

            #self.scene.mlab.view(azimuth=0, elevation=0, distance=self.current_dist)
            #self.scene.mlab.roll(0.0)

            self.current_voi = self.get_voi()
            self.iso = self.scene.mlab.pipeline.iso_surface(self.current_voi, contours=[self.int_cutoff, ], colormap='Spectral')
            
            self.set_view()

        else:
            self.iso.contour.contours = [self.int_cutoff, ]
            self.iso.update_pipeline()
            #self.scene.save_png()


    def get_voi(self):
        src = self.scene.mlab.pipeline.scalar_field(self.vs[self.vol_ind])
        src.spacing = [1, 1, 1]
        return src
 
    def set_view(self):
        if ('view' in self.dj[self.current_vol_ind]) and ('azimuth' in self.dj[self.current_vol_ind]['view']):
            # only reset view angles, use a common reset distance and contour value
            if False:
                self.scene.mlab.view(azimuth=self.dj[self.current_vol_ind]['view']['azimuth'], elevation=self.dj[self.current_vol_ind]['view']['elevation'], distance=self.current_dist, focalpoint=N.array(self.dj[self.current_vol_ind]['view']['focalpoint']), roll=self.dj[self.current_vol_ind]['view']['roll'])          # IMPORTANT: distance must be set together with angles and focalpoint
                #self.scene.mlab.roll(self.dj[self.current_vol_ind]['view']['roll'])
            else:
                # somehow distance must set after change angle, otherswise some structures cannot be correctly displayed, don't know why
                self.scene.mlab.view(azimuth=self.dj[self.current_vol_ind]['view']['azimuth'], elevation=self.dj[self.current_vol_ind]['view']['elevation'], roll=self.dj[self.current_vol_ind]['view']['roll'])
                self.scene.mlab.view(distance=self.current_dist, focalpoint=N.array(self.dj[self.current_vol_ind]['view']['focalpoint']))


            print 'set_view', self.dj[self.current_vol_ind]['view']
        

    def record_view(self):
        if 'view' not in self.dj[self.current_vol_ind]:      self.dj[self.current_vol_ind]['view'] = {}

        e = self.scene.mlab.view()
        self.dj[self.current_vol_ind]['view']['azimuth'] = e[0]
        self.dj[self.current_vol_ind]['view']['elevation'] = e[1]

        #self.dj[self.current_vol_ind]['view']['distance'] = e[2]
        self.current_dist = e[2]

        self.dj[self.current_vol_ind]['view']['focalpoint'] = e[3].tolist()
        
        self.dj[self.current_vol_ind]['view']['roll'] = self.scene.mlab.roll()

        #self.dj[self.current_vol_ind]['view']['contour'] = self.int_cutoff

        for i in range(len(self.dj)):
            if 'view' not in self.dj[i]:    self.dj[i]['view'] = {}
            self.dj[i]['view']['distance'] = self.current_dist
            self.dj[i]['view']['contour'] = self.int_cutoff

        #print 'record_view', self.dj[self.current_vol_ind]['view']


    def _save_fired(self):
        self.record_view()
        with open(self.out_file, 'w') as f:     json.dump(self.dj, f, indent=2)
        print 'views saved'


    # The layout of the dialog created
    view = TUA.View(TUA.Item('scene', editor=MCUA.SceneEditor(scene_class = MCUA.MayaviScene), height=400, width=400, show_label=False), TUA.Group('_', 'int_cutoff', 'vol_ind', '_', 'save'), resizable=True)


if __name__ == '__main__':
    main()





'''
# some test code on isosurface


import tomominer.model.util as MU
v = MU.generate_toy_model()

from mayavi import mlab

mlab.figure(bgcolor=(0, 0, 0), size=(400, 400))

src = mlab.pipeline.scalar_field(v)
src.spacing = [1, 1, 1]
src.update_image_data = True
voi = mlab.pipeline.extract_grid(src)
#voi.set(x_min=125, x_max=193, y_min=92, y_max=125, z_min=34, z_max=75)
p = mlab.pipeline.iso_surface(voi, contours=[v.max()/10,], colormap='Spectral')

p.get()
p.mlab_source.get().keys()

def picker_callback(picker):
    print mlab.view(), mlab.roll()

figure = mlab.gcf()
figure.on_mouse_pick(picker_callback)

view_org = mlab.view()
focal_point = view_org[3]
mlab.view(azimuth=-125, elevation=54, distance=100, focalpoint=focal_point)
mlab.roll(roll=-175)

mlab.draw()

p.contour.contours = [self.int_cutoff, ]
p.update_pipeline()

mlab.show()

p.remove()      # remove the isosurface from the engine,    see   http://docs.enthought.com/mayavi/mayavi/advanced_scripting.html


mlab.savefig(filename='/tmp/t.pdf')



'''





'''

# use skimage to get isosurface, and use mlab to simplify


import tomominer.model.util as MU
v = MU.generate_toy_model()

from skimage import measure
verts, faces = measure.marching_cubes(v, 0.1)

from mayavi import mlab
mesh = mlab.pipeline.triangular_mesh_source(verts[:,0], verts[:,1], verts[:,2], faces, figure=None)          # see http://nullege.com/codes/show/src@p@y@pysurfer-0.4@surfer@viz.py/2247/enthought.mayavi.mlab.pipeline.surface


# simplify surfaces
dec = mlab.pipeline.decimate_pro(mesh, figure=None)          # see       http://docs.enthought.com/mayavi/mayavi/auto/example_julia_set_decimation.html
dec.filter.feature_angle = 1
dec.filter.target_reduction = 0.95

mlab.pipeline.surface(dec)

mlab.show()


# extract mesh data
verts[0,:]
faces[0,:]

mesh.data.points
mesh.data.polys.to_array()

dec.data.points


'''


