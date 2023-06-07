#!/usr/bin/env python




'''

interactively translate and rotate a volume (template, or average), and use mayavi to plot central slices.
record rotation parameters


code similiar to
~/ln/tomominer/tomominer/geometry/surface/isosurface/plot_batch.py



~/ln/tomominer/tomominer/image/vol/plot/slicer_mayavi.py


todo: also include a scale bar

'''




import os, json, copy

import tomominer.io.file as IF
import tomominer.geometry.rotate as GR


def main():
    with open('slicer_mayavi__op.json') as f:        op = json.load(f)

    if not os.path.isdir(op['out dir']):        os.makedirs(op['out dir'])

    v = IF.read_mrc_vol(op['map file'])
    
    my_models = MyModels(v=v, rotation_file=op['rotation file out'], pad_mode=op['pad mode'], out_dir=op['out dir'])
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

    def __init__(self, v, rotation_file, out_dir, pad_mode):
        TA.HasTraits.__init__(self)

        self.v = v
        self.v_min = v.min()
        self.v_max = v.max()

        self.rotation_file = rotation_file
        self.pad_mode = pad_mode
        self.out_dir = out_dir

        if os.path.isfile(self.rotation_file):
            print 'load existing rotation from', self.rotation_file
            with open(self.rotation_file) as f:      self.rotation = json.load(f)
        else:
            self.rotation = {'angle':[0., 0., 0.], 'loc':[0., 0., 0.]}

        pi = 3.14
        self.add_trait('rot_xi', TA.Range(low=-pi, high=pi, value=self.rotation['angle'][0]))
        self.add_trait('rot_eta', TA.Range(low=-pi, high=pi, value=self.rotation['angle'][1]))
        self.add_trait('rot_omega', TA.Range(low=-pi, high=pi, value=self.rotation['angle'][2]))

        self.add_trait('trans_x', TA.Range(low=-self.v.shape[0]/2.0, high=self.v.shape[0]/2.0, value=self.rotation['loc'][0]))
        self.add_trait('trans_y', TA.Range(low=-self.v.shape[0]/2.0, high=self.v.shape[0]/2.0, value=self.rotation['loc'][1]))
        self.add_trait('trans_z', TA.Range(low=-self.v.shape[0]/2.0, high=self.v.shape[0]/2.0, value=self.rotation['loc'][2]))
 
        self.add_trait('save', TA.Button())

        self.draw_init()


    @TA.on_trait_change('rot_xi, rot_eta, rot_omega, trans_x, trans_y, trans_z')
    def update_plot(self):
        self.rotation = {'angle':[self.rot_xi, self.rot_eta, self.rot_omega], 'loc':[self.trans_x, self.trans_y, self.trans_z]}

        self.draw()


    def draw_init(self):
        vr = self.rotate()

        s = vr.shape

        m0 = N.squeeze(vr[s[0]/2,:,:])
        self.o0 = self.scene.mlab.imshow(m0, colormap='gray', vmin=self.v_min, vmax=self.v_max)
        self.o0.mlab_source.m_data.origin = [0,0,0]
        #self.o0.mlab_source.dataset.origin = [0,0,0]

        m1 = N.squeeze(vr[:,s[1]/2,:])
        self.o1 = self.scene.mlab.imshow(m1, colormap='gray', vmin=self.v_min, vmax=self.v_max)
        self.o1.mlab_source.m_data.origin = [m0.shape[0],0,0]

        m2 = N.squeeze(vr[:,:,s[2]/2])
        self.o2 = self.scene.mlab.imshow(m2, colormap='gray', vmin=self.v_min, vmax=self.v_max)
        self.o2.mlab_source.m_data.origin = [m0.shape[0]+m1.shape[0],0,0]


        #self.scene.mlab.view(azimuth=0, elevation=0, roll=0)



    def draw(self):
        vr = self.rotate()

        s = vr.shape

        m0 = N.squeeze(vr[s[0]/2,:,:])
        #self.o0.mlab_source.dataset.point_data.scalars = m0.T.flatten()
        self.o0.mlab_source.scalars = m0
        #self.o0.mlab_source.dataset.scalar_range = (m.min(), m.max())
        self.o0.update_pipeline()

        m1 = N.squeeze(vr[:,s[1]/2,:])
        self.o1.mlab_source.scalars = m1
        self.o1.update_pipeline()

        m2 = N.squeeze(vr[:,:,s[2]/2])
        self.o2.mlab_source.scalars = m2
        self.o2.update_pipeline()


    def rotate(self):

        if self.pad_mode == 'mean':
            vr = GR.rotate_pad_mean(self.v, angle=self.rotation['angle'], loc_r=self.rotation['loc'])
        elif self.pad_mode == 'zero':
            vr = GR.rotate_pad_zero(self.v, angle=self.rotation['angle'], loc_r=self.rotation['loc'])
        else:
            raise Exception('pad_mode')
        return vr


       
    def _save_fired(self):
        with open(self.rotation_file, 'w') as f:     json.dump(self.rotation, f, indent=2)
        vr = self.rotate()
        save_slice(N.squeeze(vr[vr.shape[0]/2,:,:]), os.path.join(self.out_dir, 'slice-0.png'))
        save_slice(N.squeeze(vr[:,vr.shape[1]/2,:]), os.path.join(self.out_dir, 'slice-1.png'))
        save_slice(N.squeeze(vr[:,:,vr.shape[2]/2]), os.path.join(self.out_dir, 'slice-2.png'))

        IF.put_mrc(vr, os.path.join(self.out_dir, 'vol.mrc'))

        print 'rotation parameters, rotated volume, and slices saved'


    # The layout of the dialog created
    view = TUA.View(TUA.Item('scene', editor=MCUA.SceneEditor(scene_class = MCUA.MayaviScene), height=400, width=400, show_label=False), TUA.Group('_', 'rot_xi', 'rot_eta', 'rot_omega', 'trans_x', 'trans_y', 'trans_z', '_', 'save'), resizable=True)



import matplotlib.pyplot as PLT
import matplotlib.cm as MCM


def save_slice(s, out_file):

    f, sp = PLT.subplots(nrows=1, ncols=1)

    sp.imshow( s, cmap = MCM.Greys_r )
    sp.axis('off') # clear x- and y-axes

    PLT.draw()

    PLT.savefig(out_file, bbox_inches='tight')
    
    PLT.close("all")




if __name__ == '__main__':
    main()


