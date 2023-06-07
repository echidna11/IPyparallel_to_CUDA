#!/usr/bin/env python



'''
export groud truth, predicted pattern, and contegency tables in excel format, insert figures into the tables if needed,
for EPS format figures, need to convert to PNG format

~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_high_fsc_consistency_with_ground_truth/table_excel_export.py
'''


import os, json, pickle, copy, uuid, csv
import xlsxwriter
import numpy as N

import tomominer.image.format_convert as TIF
import tomominer.doc.table.excel.insert_image as TDTEI
import tomominer.pursuit.multi.analysis.cluster_avg_high_fsc_consistency_with_ground_truth.cluster_avg_high_fsc_consistency_with_ground_truth as PMACC
import tomominer.io.file as IF
from tomominer.image.io import save_image_matplotlib

   

def main():
    with open('table_excel_export__op.json') as f:     op = json.load(f)

    with open(op['cluster_avg_high_fsc_consistency_with_ground_truth__op']) as f:       st_op = json.load(f)

    with open(op['cluster_avg_high_fsc_consistency_with_ground_truth__stat'], 'rb') as f:     stc = pickle.load(f)

    if "plot slice" in op:
        # load example subtomograms
        vs_count = {}       # we use this to pick up instances of a particular count, the purpose is to plot better center slices
        vs = {}
        for k in stc['true_lbl']:
            pid = stc['pid'][stc['true_lbl'][k]]
            if pid in vs:   continue
            if pid not in vs_count:     vs_count[pid] = 0
            vs_count[pid] += 1
            if ('selected instances' in op['plot slice']) and (pid in op['plot slice']['selected instances']) and (vs_count[pid] < op['plot slice']['selected instances'][pid]):
                print   'skipping instance', vs_count[pid], 'of', pid
                continue
            v = IF.read_mrc_vol(k)
            vs[pid] = v
    else:
        vs = None

    with open(op['ground_truth_average_isosurface__stat']) as f:     st_iso = json.load(f)
    st_iso = {int(_):st_iso[_] for _ in st_iso}

    if 'model_stat_file' in op:
        with open(op['model_stat_file']) as f:   model_stat = json.load(f)
        model_instance_count_t = model_stat['instance_count']
        model_instance_count_t = {int(_):model_instance_count_t[_] for _ in model_instance_count_t}

        with open(op['model_tomogram_subtomogram_file']) as f:       mts_op = json.load(f)
        if 'selected_models' in mts_op:
            tomogram_ids = set(mts_op['selected_models'])
        else:
            tomogram_ids = set(model_instance_count_t.keys())
 
        model_instance_count = {}
        for tid_t in tomogram_ids:
            for pid_t in model_instance_count_t[tid_t]:
                if pid_t not in model_instance_count:     model_instance_count[pid_t] = 0
                model_instance_count[pid_t] += model_instance_count_t[tid_t][pid_t]
    else:
        model_instance_count = None

    with open(op['cluster_avg_ssnr_fsc__plot__stat']) as f:     st_fsc = json.load(f)
    st_fsc = {int(_):st_fsc[_] for _ in st_fsc}


    with open(op['pattern truth correspondance']) as f:
        cr = csv.reader(f)
        ptc = {int(_[0]):str(_[1]) for _ in cr if (len(_) == 2)}         # manually defined pattern to truth correspondance


    tab = PMACC.stat_table(stc)
    PMACC.stat_table__arrange(tab)

    wb = xlsxwriter.Workbook(op['table out'])
    formats = {'bold': wb.add_format({'bold': True}), 'bg_red':wb.add_format({'bg_color': 'red'}), 'bg_green':wb.add_format({'bg_color': 'green'})}

    et_op = {'figure_width':op['figure width'], 'figure_resolution':op['figure resolution'], 'formats':formats}
    if ('plot slice' in op) and ('v_threshold_std_cutoff' in op['plot slice']):          et_op['v_threshold_std_cutoff'] = op['plot slice']['v_threshold_std_cutoff']
    export_truth(ws=wb.add_worksheet('Groud truth'), tab=tab, iso=st_iso, vs=vs, model_instance_count=model_instance_count,op=et_op)

    export_pred(tab=tab, iso=st_iso, ptc=ptc, grey_pattern_ids=(set(op['grey pattern ids']) if 'grey pattern ids' in op else None), fsc=st_fsc, figure_width=op['figure width'], figure_resolution=op['figure resolution'], ws=wb.add_worksheet('Prediction'), formats=formats, include_central_slice=(op['include central slice'] if 'include central slice' in op else False))
    export_contingency(tab=tab, ws=wb.add_worksheet('Contingency'), formats=formats, intersect_size_cutoff=op['contingency__intersect_size_cutoff'], resolution_cutoff=op['contingency__resolution_cutoff'])
    export_fp_fn(calculate_fp_fn(tab=tab), ws=wb.add_worksheet('FP FN'), formats=formats)

    wb.close()



'''
columns:    pdb id,     instance number,     isosurface,        one central slice of one subtomogram
parameters:     tab:table       iso:isosurface      vs: example subtomograms
'''
def export_truth(ws, tab, iso, vs, op, model_instance_count=None):


    if vs is not None:
        # use a vector to collect all intensities
        vv_t = []
        for pid in vs:      vv_t.extend(vs[pid].flatten().tolist())
        vv_t = N.array(vv_t)

        if 'v_threshold_std_cutoff' in op:
            v_threshold_min = vv_t.mean() - op['v_threshold_std_cutoff'] * vv_t.std()
            v_threshold_max = vv_t.mean() + op['v_threshold_std_cutoff'] * vv_t.std()
        else:
            v_threshold_min = vv_t.min()
            v_threshold_max = vv_t.max()


    ts = tab['template_stat']

    col_i = 0
    ws.write(0, col_i, 'PDB ID', op['formats']['bold']);            col_i += 1

    if model_instance_count is not None:
        ws.write(0, col_i, 'Initial copy number before particle picking', op['formats']['bold']);          col_i += 1
        ws.write(0, col_i, 'MPP input copy number after particle picking', op['formats']['bold']);       col_i += 1
    else:
        ws.write(0, col_i, 'Copy number', op['formats']['bold']);       col_i += 1

    ws.write(0, col_i, 'Isosurface', op['formats']['bold']);        col_i += 1

    if vs is not None:      
        ws.write(0, col_i, 'Subtomogram central slice example (X-Z)', op['formats']['bold'])
        col_i += 1

    imgs = []

    row = 1
    for ti in ts:
        assert ts[ti]['id'] == ti

        col_i = 0
        ws.write(row, col_i, ts[ti]['pid']);        col_i += 1

        if model_instance_count is not None:
            ws.write(row, col_i, model_instance_count[ts[ti]['pid']]);       col_i += 1

        ws.write(row, col_i, ts[ti]['instance_num']);       col_i += 1

        imgs.append({'row':row, 'col':col_i, 'file':iso[ti]['plot'], 'size_scaled':[op['figure_width'], None], 'x_offset':0, 'y_offset':3});        col_i += 1

        if vs is not None:
            slice_out = '/tmp/tmp--%s--1.png'%(str(uuid.uuid4()), )
            v = vs[ts[ti]['pid']]
            save_image_matplotlib(m=N.squeeze(v[:,v.shape[1]/2,:]), out_file=slice_out, vmin=v_threshold_min, vmax=v_threshold_max)           # we plot only the x-z slice
            imgs.append({'row':row, 'col':col_i, 'file':slice_out, 'size_scaled':[op['figure_width'], None], 'x_offset':0, 'y_offset':3});          col_i += 1
       
        row += 1

    
    TDTEI.insert_images(ws, imgs, op={'row_margin':10.0, 'col_margin':0.0})


    

'''
parameters:     tab:    iso: iso surfaces   ptc: prediction true correspondance
columns:    pattern id,     subtomogram number,     isosurface,     fsc curve
'''
def export_pred(tab, iso, ptc, grey_pattern_ids, fsc, figure_width, figure_resolution, ws, formats, include_central_slice=False):

    # find corresponding isosurface file for each pattern
    ptc_i = {ptc[_]:_ for _ in ptc }        # pdb id : pattern id

    ts = tab['template_stat']
    pid_ti = {ts[_]['pid']:_ for _ in ts}           # mapping between pid and truth id
    ti_pid = {_:ts[_]['pid'] for _ in ts}

    ci_ti = {_:pid_ti[ptc[_]] for _ in ptc}

    cs = tab['clus_stat']

    # get and process averages for plotting center slices
    vs = {}
    for ci in cs:   vs[ci] = IF.read_mrc_vol(fsc[ci]['subtomogram'])

    vs_v = []
    for ci in vs:       vs_v.extend(vs[ci].flatten().tolist())
    vs_max = N.max(vs_v)
    vs_min = N.min(vs_v)


    col_i = 0
    ws.write(0, col_i, 'Pattern ID', formats['bold'])       ;       col_i += 1
    ws.write(0, col_i, 'Copy number', formats['bold'])       ;       col_i += 1
    if include_central_slice:        ws.write(0, col_i, 'Center slice', formats['bold'])       ;       col_i += 1
    ws.write(0, col_i, 'Isosurface', formats['bold'])       ;       col_i += 1
    ws.write(0, col_i, 'FSC sum', formats['bold'])       ;       col_i += 1
    ws.write(0, col_i, 'FSC curve', formats['bold'])       ;       col_i += 1


    imgs = []
    row = 1
    for ci in cs:
        col_i = 0
        ws.write(row, col_i, cs[ci]['id'])       ;       col_i += 1
        ws.write(row, col_i, cs[ci]['size'])       ;       col_i += 1

        if include_central_slice:            imgs.append({'row':row, 'col':col_i, 'file':image_convert(generate_center_slice(v=vs[ci], v_max=vs_max, v_min=vs_min), resolution=figure_resolution), 'size_scaled':[figure_width, None], 'x_offset':0, 'y_offset':3})       ;       col_i += 1

        if (grey_pattern_ids is not None) and (cs[ci]['id'] in grey_pattern_ids):
            iso_file = iso[ci_ti[ci]]['pred'][unicode(ci)]['plot-grey']     # insert grey version of the isosurface, because the pattern is a highly mixed version of ground truth
        else:
            iso_file = iso[ci_ti[ci]]['pred'][unicode(ci)]['plot']          # insert colored version of the isosurface, color correspond to ground truth

        imgs.append({'row':row, 'col':col_i, 'file':image_convert(fi=iso_file, resolution=figure_resolution), 'size_scaled':[figure_width, None], 'x_offset':0, 'y_offset':10})       ;       col_i += 1

        ws.write(row, col_i, '%.3f'%(fsc[ci]['fsc_sum'],))       ;       col_i += 1

        imgs.append({'row':row, 'col':col_i, 'file':image_convert(fi=fsc[ci]['file'], resolution=figure_resolution), 'size_scaled':[figure_width, None], 'x_offset':0, 'y_offset':3})       ;       col_i += 1

        row += 1

    TDTEI.insert_images(ws, imgs, op={'row_margin':10.0, 'col_margin':0.0})

 


def image_convert(fi, resolution):
    if fi.endswith('.png'):
        fo = fi
    else:
        fo = os.path.join('/tmp', 'tmp--%s.png'%(str(uuid.uuid4()),))
        TIF.convert_a2ping(fi, fo, resolution=resolution)

    return fo


import tomominer.image.io as TII
def generate_center_slice(v, v_max, v_min):
    fo = os.path.join('/tmp', 'tmp--%s.png'%(str(uuid.uuid4()),))
    TII.save_image_matplotlib(N.squeeze(v[:,:,v.shape[2]/2]), out_file=fo, vmin=v_min, vmax=v_max)
    return fo


def export_contingency(tab, ws, formats, intersect_size_cutoff=None, resolution_cutoff=None):

    clus_stat = tab['clus_stat']
    template_stat = tab['template_stat']
    tab_stat = tab['tab_stat']


    # print a table of statistics
    col = 1
    for cp_i, cp in enumerate(tab['arrange']['cps']):
        ws.write(0, col, clus_stat[cp]['id'], formats['bold'])
        col += 1

    ws.write(0, col, '(Pattern ID)', formats['bold'])


    row = 1
    for ct_i, ct in enumerate(tab['arrange']['cts']):

        ws.write(row, 0, template_stat[ct]['pid'], formats['bold'])

        col = 1
        for cp_i, cp in enumerate(tab['arrange']['cps']):
            cel_str = '%d / %.1f'%(tab_stat[ct][cp]['count_cluster'], tab_stat[ct][cp]['resolution'])
            if (intersect_size_cutoff is None) or (tab_stat[ct][cp]['count_cluster'] < intersect_size_cutoff):
                ws.write(row, col, cel_str)
            else:
                if tab_stat[ct][cp]['resolution'] >= resolution_cutoff:
                    ws.write(row, col, cel_str, formats['bg_red'])
                else:
                    ws.write(row, col, cel_str, formats['bg_green'])

            col += 1

        row += 1

    ws.write(row, 0, '(PDB ID)', formats['bold'])



'''

code for calculating false positive (FP) and false negative (FN) from stat.csv.

Suppose we have one to one match of a particular complex and particular pattern, this is indicated in the diagonal of re-arranged contingency table. Suppose a complex C matches a pattern P. The FN of C is the number of instances of C that are not included into P. The FP of P is the number of instances of pattern P that does not belong to C.
'''

def calculate_fp_fn(tab):
    clus_stat = tab['clus_stat']
    template_stat = tab['template_stat']
    tab_stat = tab['tab_stat']


    pred_size = {}
    pred_ids = {}
    for cp_i, cp in enumerate(tab['arrange']['cps']):
        pred_ids[cp_i] = clus_stat[cp]['id']
        pred_size[cp_i] = clus_stat[cp]['size']


    true_size = {}
    true_ids = {}
    for ct_i, ct in enumerate(tab['arrange']['cts']):
        true_ids[ct_i] = template_stat[ct]['pid']
        true_size[ct_i] = template_stat[ct]['instance_num']

    t = N.zeros((len(tab['arrange']['cts']), len(tab['arrange']['cps'])))
    for ct_i, ct in enumerate(tab['arrange']['cts']):
        for cp_i, cp in enumerate(tab['arrange']['cps']):
            t[ct_i, cp_i] = tab_stat[ct][cp]['count_cluster']


    tp = {}         # true positive
    tpr = {}        # true positive rate
    fp = {}         # false positive
    fdr = {}        # False Discovery Rate
    fn = {}         # false negative
    fnr = {}        # False Negative Rate

    for i in range(min(t.shape)):
        # for a match to be valid, we require the diagonal value to be maximum
        if t[i,i] < max(t[i,:]):        continue
        if t[i,i] < max(t[:,i]):        continue

        tp[i] = t[i,i]

        fp[i] = pred_size[i] - tp[i]
        fn[i] = true_size[i] - tp[i]

        tpr[i] = float(tp[i])/true_size[i]
        fnr[i] = float(fn[i])/true_size[i]
        fdr[i] = float(fp[i])/pred_size[i]

    return {'fp':fp, 'fn':fn, 'tpr':tpr, 'fdr':fdr, 'fnr':fnr, 'pred_size':pred_size, 'true_size':true_size, 'true_ids':true_ids, 'pred_ids':pred_ids}



def export_fp_fn(fp_fn, ws, formats):

    fp = fp_fn['fp']
    fn = fp_fn['fn']
    fdr = fp_fn['fdr']
    fnr = fp_fn['fnr']

    true_ids = fp_fn['true_ids']
    pred_ids = fp_fn['pred_ids']

    ws.write_row(0, 0, ['PDB ID', 'False Negative', 'False Negative Rate', 'Pattern ID', 'False Positive', 'False Discovery Rate'], formats['bold'])

    row = 1
    for i in fp:
        ws.write_row(row, 0, [true_ids[i], fn[i], '%0.2f'%(fnr[i],), pred_ids[i], fp[i], '%0.2f'%(fdr[i],)])
        row += 1


if __name__ == '__main__':
    main()




