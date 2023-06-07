#!/usr/bin/env python





'''
plot and save all resulting isosurfaces, 
and also save central slice of density maps using mayavi (for simplification, we do not plot them, but just save them)

There will be an adjust bar to adjust isovalue, and a button to save picture and go to next average

todo: in future, consider to decompose into several programs


~/ln/tomominer/tomominer/average/genetic_algorithm/config/analysis/from_multi_pursuit__batch_plot.py
'''


import os, json, pickle, copy

import tomominer.io.file as IF


def main():
    with open('from_multi_pursuit__batch_plot__op.json') as f:        op = json.load(f)

    if not os.path.isdir(op['out dir']):    os.makedirs(op['out dir'])

    with open(op['config file']) as f:       conf = json.load(f)

    print 'collecting related data',
    data = []
    for clus in conf['clusters']:
        for c in conf['clusters'][clus]:
            d = collect(c)
            for dd in d:        dd['cluster'] = int(clus)

            data.extend(d)
            print len(data),

    export_statistics(data, op['stat file out'])

    data = data_reorganize(data)
    save_central_slices(data=data, out_dir=op['out dir'])

    if not op['intensity positive']:
        for d in data:      d['vol'] = -d['vol']

    my_models = MyModels(data=data, opacity=op['opacity'], out_dir=op['out dir'])
    my_models.configure_traits()



'''
collect the average from iteration 0, and the best iteration (if different than 0)
'''
def collect(c):

    with open(os.path.join(c['out_dir'], 'file_stat.json')) as f:       fs = json.load(f)
    
    p = fs['passes']
    p = {int(_):p[_] for _ in p}
    for i in p:     p[i]['pass_i'] = i


    data = []
    d = copy.deepcopy(c)
    d['p'] = collect_single_pass(p[0])
    data.append(d)

    if fs['best_pass_i'] > 0:
        d = copy.deepcopy(c)
        d['p'] = collect_single_pass(p[fs['best_pass_i']])
        data.append(d)

    return data


def collect_single_pass(p):
    p = copy.deepcopy(p)
    with open(p['cluster_averaging_file'], 'rb') as f:      ca = pickle.load(f)
    p['subtomogram'] = IF.read_mrc_vol(ca['template_keys'][0]['subtomogram'])

    with open(p['fsc_stat_file']) as f:     fsc_stat = json.load(f)
    fsc_stat = {int(_):fsc_stat[_] for _ in fsc_stat}
    

    p['fsc_stat'] = fsc_stat[p['pass_i']]

    return p


import csv
def export_statistics(data, out_file):
    print 'export_statistics()'
    s = []
    for i, d in enumerate(data):
        s.append([i, 'cluster', d['cluster'], 
            'aligned', d['aligned'], 
            'smooth', d['smooth'], 
            'repeat', d['repeat_i'], 
            'pass', d['p']['pass_i'], 
            'size', d['p']['fsc_stat']['cluster_size'], 
            'fsc', d['p']['fsc_stat']['fsc_sum']])

    with open(out_file, 'w') as f:
        sw = csv.writer(f, delimiter='\t')
        sw.writerows(s)



'''
only keep information needed for plotting isosurface and information
'''
def data_reorganize(data):
    print 'data_reorganize()'

    data_new = []
    for d in data:
        desc = 'clus-%03d'%(d['cluster'], )
        desc += '--'
        desc += 'aligned' if d['aligned'] else 'noalign'
        desc += '--'
        desc += 'smooth' if d['smooth'] else 'nosmooth'
        desc += '--'
        desc += 'repeat-%03d'%(d['repeat_i'], )
        desc += '--'
        desc += 'pass-%03d'%(d['p']['pass_i'],)

        dt = {}
        dt['desc'] = desc
        dt['vol'] = d['p']['subtomogram']
        
        data_new.append(dt)

    return data_new



import tomominer.pursuit.multi.analysis.instance.instance_plotting as TPMAII

def save_central_slices(data, out_dir):
    print 'save_central_slices()'

    for i, d in enumerate(data):
        slice_out_file = os.path.join(out_dir, '%04d--%s--slice.png'%(i,d['desc']))
        TPMAII.save_slices(d['vol'], slice_out_file)


import numpy as N

import traits.api as TA
import traitsui.api as TUA
import mayavi.modules.api as TMA
import mayavi.core.api as MCA
import mayavi.core.ui.api as MCUA


class MyModels(TA.HasTraits):
    scene = TA.Instance(MCUA.MlabSceneModel, ())

    plot = TA.Instance(MCA.PipelineBase)

    def __init__(self, data, opacity, out_dir):
        TA.HasTraits.__init__(self)

        self.d = copy.deepcopy(data)
        self.opacity = opacity
        self.out_dir = out_dir

        v_max_min = []
        for d in data:    v_max_min.extend([d['vol'].max(), d['vol'].min()])

        int_cutoff_default = float(N.mean(v_max_min))

        self.add_trait('int_cutoff', TA.Range(low=float(N.floor(min(v_max_min))), high=float(N.ceil(max(v_max_min))), value=int_cutoff_default))  # this trait's range need to be determined at initialization stage, it is not a statistic trait
        
        self.add_trait('save', TA.Button())
        self.add_trait('next', TA.Button())


        self.vol_ind = 0
        self.current_vol_ind = -1
        self.iso = None
        self.iso_count = 0

        self.update_plot()


    @TA.on_trait_change('int_cutoff')
    def update_plot(self):
        if self.current_vol_ind != self.vol_ind:
            self.current_vol_ind = self.vol_ind
            if self.iso is not None:                self.iso.remove()

            self.current_voi = self.get_voi()
            self.iso = self.scene.mlab.pipeline.iso_surface(self.current_voi, contours=[self.int_cutoff, ], colormap='Spectral', opacity=self.opacity)
            
        else:
            self.iso.contour.contours = [self.int_cutoff, ]
            self.iso.update_pipeline()


    def get_voi(self):
        v = self.d[self.vol_ind]['vol']

        src = self.scene.mlab.pipeline.scalar_field(v)
        src.spacing = [1, 1, 1]
        return src
 

    def _save_fired(self):
        iso_out_file = os.path.join(self.out_dir, '%04d--%s--iso--%03d.png'%(self.vol_ind, self.d[self.vol_ind]['desc'], self.iso_count))            # we save multiple iso files to display different aspects
        self.scene.save_png(iso_out_file)
        print 'saved', iso_out_file
        self.iso_count += 1

    def _next_fired(self):
        if self.vol_ind >= (len(self.d) - 1):       return

        self.vol_ind += 1
        self.iso_count = 0
        self.update_plot()


    # The layout of the dialog created
    view = TUA.View(TUA.Item('scene', editor=MCUA.SceneEditor(scene_class = MCUA.MayaviScene), height=400, width=400, show_label=False), TUA.Group('_', 'int_cutoff', '_', 'save', 'next'), resizable=True)



if __name__ == '__main__':
    main()




'''
for similiar code, see

~/ln/tomominer/tomominer/geometry/surface/isosurface/plot_batch.py
~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_high_fsc_consistency_with_ground_truth/ground_truth_average_isosurface.py
/home/rcf-47/mxu/ln/tomominer/tomominer/image/vol/plot/slicer_mayavi.py

~/ln/tomominer/tomominer/simulation/tomogram/analysis/class_isosurface_plot.py

'''

