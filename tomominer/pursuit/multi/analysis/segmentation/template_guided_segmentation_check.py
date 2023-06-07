#!/usr/bin/env python


'''
given a pursuit parameter and results, perform segmentation to templates, and use template to guide subtomogram segmentations
this is done by calling exising routine in pursuit.multi.util.template_guided_segmentation() and template_segmentation__single()


~/ln/tomominer/tomominer/pursuit/multi/analysis/segmentation/template_guided_segmentation_check.py
'''


import os, shutil, sys, copy, json
import cPickle as pickle
import tomominer.pursuit.multi.util as PMU
import tomominer.common.obj as CO
import tomominer.io.file as IV

def main():
    
    self = CO.Object()
    self.pool = None

    with open('template_guided_segmentation_check__op.json') as f:      op = json.load(f)
    with open(op['pursuit_op_file']) as f:      pop = json.load(f)


    segmentation_op = copy.deepcopy(pop['segmentation'])
    segmentation_op['density_positive'] = pop['density_positive']


    template_segmentation_op = copy.deepcopy(segmentation_op)
    if ('segmentation' in pop['template']) and ('normalize_and_take_abs' in pop['template']['segmentation']) and (pop['template']['segmentation']['normalize_and_take_abs']):      template_segmentation_op['normalize_and_take_abs'] = True
    print 'template_segmentation_op', template_segmentation_op


    segmentation_tg_op = copy.deepcopy(segmentation_op)     # options for segmenting a subtomogram, guided by template
    if ('guided_segmentation' in pop['template']) and ('gaussian_smooth_sigma' in pop['template']['guided_segmentation']):            segmentation_tg_op['gaussian_smooth_sigma'] =  pop['template']['guided_segmentation']['gaussian_smooth_sigma']
    print 'segmentation_tg_op', segmentation_tg_op



    with open(op['file_stat_file']) as f:      fs = json.load(f)
    fs['passes'] = {int(_):fs['passes'][_] for _ in fs['passes']}


    if 'pass_i_seleted' not in op:      op['pass_i_seleted'] = fs['pass_i_current']
    fsp = fs['passes'][op['pass_i_seleted']]

    with open(fsp['data_json_file']) as f:  dj = json.load(f)


    if os.path.isdir(op['out_dir']):    shutil.rmtree(op['out_dir'])

    #-----------------------------------------------------------------
    # first copy templates to a temporary dir
    with open(fsp['cluster_average_select_file'], 'rb') as f:       cas = pickle.load(f)

    tem_dir = os.path.join(op['out_dir'], 't')
    os.makedirs(tem_dir)

    tk = {}
    for c in cas['selected_templates_common_frame']:
        tk[c] = {}
        tk[c]['subtomogram'] = os.path.join(tem_dir, os.path.basename(cas['selected_templates_common_frame'][c]['subtomogram']))
        shutil.copyfile(cas['selected_templates_common_frame'][c]['subtomogram'], tk[c]['subtomogram'])
    del c

    # segment templates
    PMU.template_segmentation(self=self, tk=tk, op=template_segmentation_op)

    #--------------------------------------------------------------
    # template guided segmentation
    subt_dir = os.path.join(op['out_dir'], 's')
    os.makedirs(subt_dir)

    dj = sorted(dj, key=lambda _:_['score'], reverse=True)      # order dj according to alignment scores

    cls_cnt = {}            # use this to restrict the number of subtomograms in each class
    for i,d in enumerate(dj):
        print '\r', i,          ;       sys.stdout.flush()

        c = d['template']['id']
        if c not in cls_cnt:        cls_cnt[c] = 0
        if cls_cnt[c] > op['max_num_each_class']:       continue
        cls_cnt[c] += 1

        file_name_root = os.path.join(subt_dir, '%03d-%05d'%(c, i))
        v = IV.read_mrc_vol(d['subtomogram'])
        IV.put_mrc(v, file_name_root + '.mrc')

        dt = copy.deepcopy(d)
        dt['template'] = tk[c]
        s = PMU.align_to_templates__segment(rec=dt, v=v, segmentation_op=template_segmentation_op)
        IV.put_mrc(s['v'], file_name_root + '-seg.mrc')
        IV.put_mrc(s['phi_mr'], file_name_root + '-phim.mrc')









if __name__ == '__main__':
    main()

