#!/usr/bin/env python



'''
this is used to interactivelly set view angles, focal points, and common distance for plotting individual isosurfaces, 
once this is set, we can have another program to plot all isosurfaces


~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/template_isosurface__set_view_angle.py
'''


import os, json, copy
import cPickle as pickle


def main():
    with open('template_isosurface__set_view_angle__op.json') as f:       op = json.load(f)

    if os.path.isfile(op['view file out']):
        # when we already have a record of previous results
        print 'using existing', op['view file out']
        with open(op['view file out']) as f:    dj = json.load(f)
    else:
        print 'loading', op['data file in']
        with open(op['data file in']) as f:   dj = json.load(f)


    print 'loading', op['mesh file in']
    with open(op['mesh file in'], 'rb') as f:       surf = pickle.load(f)

    if 'color file' in op:
        print 'loading', op['color file']
        with open(op['color file']) as f:       colors = json.load(f)
        colors = [tuple(colors[d['subtomogram']]['color']) for d in dj]
    else:
        colors = None

    
    my_models = MyModels(dj=dj, surf=surf, colors=colors, opacity=op['opacity'], view_file=op['view file out'])
    my_models.configure_traits()



import numpy as N

import traits.api as TA
import traitsui.api as TUA
import mayavi.modules.api as TMA
import mayavi.core.api as MCA
import mayavi.core.ui.api as MCUA
from mayavi import mlab

class MyModels(TA.HasTraits):
    scene = TA.Instance(MCUA.MlabSceneModel, ())

    plot = TA.Instance(MCA.PipelineBase)

    def __init__(self, dj=None, surf=None, colors=None, opacity=None, view_file=None):
        TA.HasTraits.__init__(self)

        self.dj = copy.deepcopy(dj)
        self.surf = surf
        self.colors = colors
        self.opacity = opacity
        self.view_file = view_file

        self.add_trait('vol_ind', TA.Range(low=0, high=len(self.dj)-1, value=0))

        self.add_trait('save', TA.Button())

        self.current_dist = None
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

            i = self.current_vol_ind
            surf = self.surf
            fv = surf[i]['fv']
            self.iso = mlab.triangular_mesh(fv['vertices'][:,0], fv['vertices'][:,1], fv['vertices'][:,2], fv['faces'], color=(self.colors[i] if self.colors is not None else (0.8, 0.8, 0.8)), opacity=self.opacity)

            self.set_view()


 
    def set_view(self):
        if 'view' in self.dj[self.current_vol_ind]:
            # only reset view angles, use a common reset distance and contour value
            self.scene.mlab.view(azimuth=self.dj[self.current_vol_ind]['view']['azimuth'], elevation=self.dj[self.current_vol_ind]['view']['elevation'], distance=self.current_dist, focalpoint=N.array(self.dj[self.current_vol_ind]['view']['focalpoint']), roll=self.dj[self.current_vol_ind]['view']['roll'])          # IMPORTANT: distance must be set together with angles and focalpoint
            #self.scene.mlab.roll(self.dj[self.current_vol_ind]['view']['roll'])

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
            if 'view' not in self.dj[i]:    continue
            self.dj[i]['view']['distance'] = self.current_dist

        #print 'record_view', self.dj[self.current_vol_ind]['view']


    def _save_fired(self):
        self.record_view()
        with open(self.view_file, 'w') as f:     json.dump(self.dj, f, indent=2)
        print 'views saved'


    # The layout of the dialog created
    view = TUA.View(TUA.Item('scene', editor=MCUA.SceneEditor(scene_class = MCUA.MayaviScene), height=400, width=400, show_label=False), TUA.Group('_', 'vol_ind', '_', 'save'), resizable=True)




if __name__ == '__main__':
    main()



'''
related code

~/ln/tomominer/tomominer/geometry/surface/isosurface/plot_batch.py

~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_high_fsc_consistency_with_ground_truth/ground_truth_average_isosurface.py

~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/template_isosurface.py

'''

