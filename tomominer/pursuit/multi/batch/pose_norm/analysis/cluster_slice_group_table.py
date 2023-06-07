#!/usr/bin/env python




'''
generate an excel table 
one sheet containing cluster id, cluster size, three center slices cluster average, 
another sheet containing grouping of clusters

~/ln/tomominer/tomominer/pursuit/multi/batch/pose_norm/analysis/cluster_slice_group_table.py
'''


import os, sys, json, xlsxwriter, uuid
import numpy as N

import tomominer.io.file as IF
from tomominer.image.io import save_image_matplotlib
from tomominer.pursuit.multi.analysis.cluster_avg_high_fsc_consistency_with_ground_truth.table_excel_export import image_convert
import tomominer.doc.table.excel.insert_image as TDTEI



def main():
    with open('cluster_slice_group_table__op.json') as f:       op = json.load(f)

    with open(op['cluster stat']) as f:     cs = json.load(f)
    cs = {_['cluster_label']:_ for _ in cs}

    with open(op['cluster groups']) as f:     cg = json.load(f)
    cg = {_['id']:_ for _ in cg}

    cgd = {}        # mapping from cluster to gourp id
    for g in cg:
        cgt = cg[g]
        for c in cgt['clusters']:
            assert c not in cgd
            cgd[c] = cgt['id']


    # re-order cs accourding to groups
    cs_new = []
    for g in cg:
        cgt = cg[g]
        for c in sorted(cgt['clusters']):
            cs_new.append(cs[c])

    for c in cs:
        if c in cgd:    continue
        cs_new.append(cs[c])


    wb = xlsxwriter.Workbook(op['table out'])
    formats = {'bold': wb.add_format({'bold': True})}
    ws = wb.add_worksheet('Pose normalization and filtering')

    ws.write_row(0, 0, ['Group ID', 'Cluster ID', 'Cluster size', 'x-y center slice', 'x-z center slice', 'y-z center slice'], formats['bold'])

    row = 1
    imgs = []
    for i, cst in enumerate(cs_new):
        print '\r', i, '          ',        ;       sys.stdout.flush()

        c = cst['cluster_label']

        ws.write(row, 0, str(cgd[c]) if c in cgd else ' ')
        ws.write(row, 1, c)
        ws.write(row, 2, cst['size'])

        sf = save_tmp_slices(s=cst['subtomogram'], v_threshold=float(op['v_threshold']))

        imgs.append({'row':row, 'col':3, 'file':image_convert(fi=sf[2], resolution=op['figure resolution']), 'size_scaled':[op['figure width'], None], 'x_offset':0, 'y_offset':10})
        imgs.append({'row':row, 'col':4, 'file':image_convert(fi=sf[1], resolution=op['figure resolution']), 'size_scaled':[op['figure width'], None], 'x_offset':0, 'y_offset':10})
        imgs.append({'row':row, 'col':5, 'file':image_convert(fi=sf[0], resolution=op['figure resolution']), 'size_scaled':[op['figure width'], None], 'x_offset':0, 'y_offset':10})

        row += 1


    TDTEI.insert_images(ws, imgs, op={'row_margin':10.0, 'col_margin':0.0})
    wb.close()

    for img in imgs:        os.remove(img['file'])



def save_tmp_slices(s, v_threshold):
    v = IF.read_mrc_vol(s)
    v -= v.mean()
    v /= v.std()

    sf = {}

    slice_tmp_root = os.path.join('/tmp', 'tmp--' + str(uuid.uuid4()))
    for d in range(3):      sf[d] = slice_tmp_root + '--%d.png'%(d,)

    save_image_matplotlib(m=N.squeeze(v[v.shape[0]/2,:,:]).T, out_file=sf[0], vmin=-v_threshold, vmax=v_threshold)
    save_image_matplotlib(m=N.squeeze(v[:,v.shape[1]/2,:]).T, out_file=sf[1], vmin=-v_threshold, vmax=v_threshold)
    save_image_matplotlib(m=N.squeeze(v[:,:,v.shape[2]/2]).T, out_file=sf[2], vmin=-v_threshold, vmax=v_threshold)

    return sf



if __name__ == '__main__':
    main()



'''

related code

~/ln/tomominer/tomominer/image/vol/plot/mrc_center_slice.py

~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/template_stat_table_excel_export.py

~/ln/tomominer/tomominer/pursuit/multi/batch/pose_norm/config_prepare.py

'''

