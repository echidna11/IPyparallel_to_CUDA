#!/usr/bin/env python




'''
export a excel table that contains following information
pattern id, whether the pattern is refined,  pattern subtomogram number, pattern isosurface, pattern SSNR based FSC


~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/template_stat_table_excel_export.py

'''

import os, json, pickle, csv, xlsxwriter
import numpy as N

from tomominer.pursuit.multi.analysis.cluster_avg_high_fsc_consistency_with_ground_truth.table_excel_export import image_convert
import tomominer.doc.table.excel.insert_image as TDTEI

def main():

    with open('template_stat_table_excel_export__op.json') as f:      op = json.load(f)

    # get comments
    comments = {}
    if 'comments' in op:
        with open(op['comments'], 'rb') as f:
            cr = csv.reader(f)
            for row in cr:
                if len(row) != 2:       continue
                comments[row[0]] = row[1]


    with open(op['fsc collect'], 'rb') as f:        co = pickle.load(f)         # use this to get the number of subtomograms for each pattern, and whether pattern is generated through pursuit or refinment

    cod = {}
    for fs_fn in co:
        for c in co[fs_fn]:
            t = co[fs_fn][c]
            cod[t['template']['subtomogram']] = t

    with open(op['fsc plot stat']) as f:      fp = json.load(f)

    with open(op['gold standard fsc config stat']) as f:      gs_st = json.load(f)
    with open(op['gold standard fsc plot stat']) as f:      gs_fp = json.load(f)
    gs_fp = {int(_):gs_fp[_] for _ in gs_fp}


    # reformat gold standard FSC information
    gs_plot_info = {}
    for gs_st_t in gs_st:
        k = gs_st_t['template']['subtomogram']
        c = gs_st_t['count']
        
        # format resolution string
        ress = gs_fp[c]['plot_res']['resolutions']
        ress = sorted(ress, key=lambda _:(-_['fsc_cutoff']))
        res_str = ' / '.join(['%.01f'%(_['resolution'],) for _ in ress])        # construct a combined resolution string, seperated using /

        gs_plot_info[k] = {'file':gs_fp[c]['file'], 'resolution_str':res_str}


    with open(op['iso plot stat']) as f:     ip = json.load(f)


    # exporting table
    wb = xlsxwriter.Workbook(op['table out'])
    formats = {'bold': wb.add_format({'bold': True})}
    ws = wb.add_worksheet('Predicted patterns')

    ws.write_row(0, 0, ['Pattern ID', 'Refined using single pattern pursuit', 'Isosurface', 'Comment', 'Subtomogram number', 'SSNR-FSC sum', 'SSNR-FSC curve', 'Resolution (nm)', 'Gold standard FSC curve'], formats['bold'])

    row = 1
    imgs = []
    for i, ipt in enumerate(ip):
        k = ipt['subtomogram']

        ws.write(row, 0, i)
        ws.write(row, 1, 'Y' if (cod[k]['mode'] == 'averaging') else ' ')
        imgs.append({'row':row, 'col':2, 'file':image_convert(fi=ipt['isosurface']['file'], resolution=op['figure resolution']), 'size_scaled':[op['figure width'], None], 'x_offset':0, 'y_offset':10})
        ws.write(row, 3, comments[k] if k in comments else ' ')
        ws.write(row, 4, len(cod[k]['dj']))
        ws.write(row, 5, '%.2f'%(N.sum(cod[k]['fsc']),))
        imgs.append({'row':row, 'col':6, 'file':image_convert(fi=fp[k]['file'], resolution=op['figure resolution']), 'size_scaled':[op['figure width'], None], 'x_offset':0, 'y_offset':10})
        ws.write(row, 7, gs_plot_info[k]['resolution_str'])
        imgs.append({'row':row, 'col':8, 'file':image_convert(fi=gs_plot_info[k]['file'], resolution=op['figure resolution']), 'size_scaled':[op['figure width'], None], 'x_offset':0, 'y_offset':10})

        row += 1

    TDTEI.insert_images(ws, imgs, op={'row_margin':10.0, 'col_margin':0.0})
    wb.close()

    for img in imgs:        os.remove(img['file'])



if __name__ == '__main__':
    main()



'''
related code

~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_high_fsc_consistency_with_ground_truth/table_excel_export.py

~/ln/tomominer/tomominer/pursuit/multi/analysis/statistics/collect_fsc_across_multiple_experiments.py
~/ln/tomominer/tomominer/pursuit/multi/analysis/statistics/collect_fsc_across_multiple_experiments__plot.py

'''

