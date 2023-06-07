#!/usr/bin/env python




'''
plot isosurfaces, their instance numbers, and resolution into a single figure


~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/template_isosurface_stat_combined_plot.py
'''



import os, json, pickle, csv
import numpy as N


def main():

    with open('template_isosurface_stat_combined_plot__op.json') as f:      op = json.load(f)

    with open(op['template_stat_table_excel_export__op.json']) as f:      tab_op = json.load(f)


    with open(tab_op['fsc collect'], 'rb') as f:        co = pickle.load(f)         # use this to get the number of subtomograms for each pattern, and whether pattern is generated through pursuit or refinment

    cod = {}
    for fs_fn in co:
        for c in co[fs_fn]:
            t = co[fs_fn][c]
            cod[t['template']['subtomogram']] = t


    with open(tab_op['gold standard fsc config stat']) as f:      gs_st = json.load(f)
    with open(tab_op['gold standard fsc plot stat']) as f:      gs_fp = json.load(f)
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


    with open(os.path.join(os.path.dirname(op['template_stat_table_excel_export__op.json']), tab_op['iso plot stat'])) as f:     ip = json.load(f)


    # plotting
    dpi = float(op['dpi'])
    subplot_size = op['subplot_size']

    import matplotlib.pyplot as PLT
    import matplotlib.image as MI
    
    fig = PLT.figure(figsize=(len(ip)*subplot_size[0]/dpi, subplot_size[1]/dpi), dpi=dpi)

    for i, ipt in enumerate(ip):
        k = ipt['subtomogram']
        
        img = MI.imread(ipt['isosurface']['file'])

        ax = PLT.subplot2grid((1, len(ip)+1), (0, i))
        PLT.imshow(img, origin='upper')
        PLT.text(img.shape[0]/2, 0, str(i))
        PLT.text(img.shape[0]/2, img.shape[1], str(len(cod[k]['dj'])) + '\n' + gs_plot_info[k]['resolution_str'])
        PLT.axis('off')

    PLT.savefig(op['out_file'], bbox_inches='tight')
    PLT.close('all')



if __name__ == '__main__':
    main()





'''

modified from
~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/template_stat_table_excel_export.py

and 

~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_high_fsc_consistency_with_ground_truth/isosurface_stat_combine_plot.py

'''


