#!/usr/bin/env python



'''
use this program to plot isusurfaces of ground truth and predicted patterns and statistics into combined figure(s)


run following using proot

~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_high_fsc_consistency_with_ground_truth/isosurface_stat_combine_plot.py

'''


import os, json, pickle, csv


import tomominer.pursuit.multi.analysis.cluster_avg_high_fsc_consistency_with_ground_truth.cluster_avg_high_fsc_consistency_with_ground_truth as PMACC
import tomominer.io.file as IF


def main():
    with open('isosurface_stat_combine_plot__op.json') as f:     op = json.load(f)

    with open(op['table_excel_export__op.json']) as f:     tab_op = json.load(f)

    with open(tab_op['cluster_avg_high_fsc_consistency_with_ground_truth__op']) as f:       st_op = json.load(f)

    with open(tab_op['cluster_avg_high_fsc_consistency_with_ground_truth__stat'], 'rb') as f:     stc = pickle.load(f)

    with open(tab_op['ground_truth_average_isosurface__stat']) as f:     st_iso = json.load(f)
    st_iso = {int(_):st_iso[_] for _ in st_iso}


    with open(os.path.join(os.path.dirname(op['table_excel_export__op.json']), tab_op['pattern truth correspondance'])) as f:
        cr = csv.reader(f)
        ptc = {int(_[0]):str(_[1]) for _ in cr if (len(_) == 2)}         # manually defined pattern to truth correspondance

    tab = PMACC.stat_table(stc)
    PMACC.stat_table__arrange(tab)


    from tomominer.pursuit.multi.analysis.cluster_avg_high_fsc_consistency_with_ground_truth.table_excel_export import calculate_fp_fn
    pid_hits = export_pred_true(tab=tab, iso=st_iso, ptc=ptc, fp_fn=calculate_fp_fn(tab), intersect_size_cutoff=tab_op['contingency__intersect_size_cutoff'], grey_pattern_ids=(set(tab_op['grey pattern ids']) if 'grey pattern ids' in op else None), out_file=op['pred_true_file'], dpi=op['dpi'], subplot_size=op['subplot_size'])

    export_miss_true(tab=tab, iso=st_iso, pid_hits=pid_hits, out_file=op['miss_true_file'], dpi=op['dpi'], subplot_size=op['subplot_size'])


'''
export figure on missed true complexes
'''
def export_miss_true(tab, iso, pid_hits, out_file, dpi, subplot_size):
    dpi = float(dpi)

    ts = tab['template_stat']

    ti_miss = [_ for _ in ts if ts[_]['pid'] not in pid_hits]

    import matplotlib.pyplot as PLT
    import matplotlib.image as MI
    
    fig = PLT.figure(figsize=(subplot_size[0]/dpi, len(ti_miss)*subplot_size[1]/dpi), dpi=dpi)
    
    for row_i, ti in enumerate(ti_miss):

        img = MI.imread(iso[ti]['plot'])

        ax = PLT.subplot2grid((len(ti_miss), 1), (row_i, 0))
        PLT.imshow(img, origin='upper')
        PLT.text(img.shape[0]/2, img.shape[1], ts[ti]['pid'] + '\n' + str(ts[ti]['instance_num']))
        PLT.axis('off')

    PLT.savefig(out_file, bbox_inches='tight')
    PLT.close('all')


'''
export predicted patterns, then followed by the true complexes that matches the predicted pattern
'''
def export_pred_true(tab, iso, fp_fn, ptc, intersect_size_cutoff, grey_pattern_ids, out_file, subplot_size, dpi):
    cs = tab['clus_stat']

    tpr = fp_fn['tpr']
    fdr = fp_fn['fdr']
    true_ids = fp_fn['true_ids']
    pred_ids = fp_fn['pred_ids']
    tpr_d = {pred_ids[_]:tpr[_] for _ in range(len(tpr))}
    fdr_d = {pred_ids[_]:fdr[_] for _ in range(len(fdr))}   
     


    # find corresponding isosurface file for each pattern
    ptc_i = {ptc[_]:_ for _ in ptc }        # pdb id : pattern id

    ts = tab['template_stat']
    pid_ti = {ts[_]['pid']:_ for _ in ts}           # mapping between pid and truth id
    ti_pid = {_:ts[_]['pid'] for _ in ts}

    ci_ti = {_:pid_ti[ptc[_]] for _ in ptc}

    import matplotlib.pyplot as PLT
    import matplotlib.image as MI
   
    pid_hits = set()

    tab_stat = tab['tab_stat']
    pid_intersection_count = {}
    for ci in cs:
        # plot true complex with decreasing order of intersection
        pid_intersection_count_t = []
        for ct_i, ct in enumerate(tab['arrange']['cts']):
            count_cluster_t = tab_stat[ct][ci]['count_cluster']
            if count_cluster_t < intersect_size_cutoff:     continue
            pid_intersection_count_t.append( {'ct_i':ct_i, 'ct':ct, 'pid':ts[ct]['pid'], 'count':count_cluster_t} )
        if len(pid_intersection_count_t) == 0:    continue

        pid_intersection_count[ci] = sorted(pid_intersection_count_t, key=lambda _:_['count'], reverse=True)

  
    def plot_set(ci_set, row_i):
        for ci in ci_set:
            if (grey_pattern_ids is not None) and (cs[ci]['id'] in grey_pattern_ids):
                iso_file = iso[ci_ti[ci]]['pred'][unicode(ci)]['plot-grey']
            else:
                iso_file = iso[ci_ti[ci]]['pred'][unicode(ci)]['plot']


            # plot predicted pattern
            img = MI.imread(iso_file)

            ax = PLT.subplot2grid((len(cs), len(ts)+1), (row_i, 0))
            PLT.imshow(img, origin='upper')
            if ci in fdr_d:
                PLT.text(img.shape[0]/2, img.shape[1], 'pattern %d\n%d (%d%%, %d%%)'%(ci, cs[ci]['size'], int(round(fdr_d[ci]*100)), int(round(tpr_d[ci]*100))))
            else:
                PLT.text(img.shape[0]/2, img.shape[1], 'pattern %d\n%d'%(ci, cs[ci]['size']))
            PLT.axis('off')

            col_i = 1
            for pid_intersection_count_t in pid_intersection_count[ci]:
                img = MI.imread(iso[pid_intersection_count_t['ct']]['plot'])
                ax = PLT.subplot2grid((len(cs), len(ts)+1), (row_i, col_i))
                PLT.imshow(img, origin='upper')
                PLT.text(img.shape[0]/2, img.shape[1], '%s\n%d'%(pid_intersection_count_t['pid'], ts[pid_intersection_count_t['ct']]['instance_num']))
                PLT.axis('off')

                pid_hits.add(pid_intersection_count_t['pid'])

                col_i += 1

            row_i += 1
        return row_i


    fig = PLT.figure(figsize=(len(ts)*subplot_size[0]/dpi, len(cs)*subplot_size[1]/dpi), dpi=dpi)

    row_i = 0
    row_i = plot_set([_ for _ in cs if len(pid_intersection_count[_]) == 1], row_i)
    plot_set([_ for _ in cs if len(pid_intersection_count[_]) > 1], row_i)

    PLT.savefig(out_file, bbox_inches='tight')
    PLT.close('all')


    return pid_hits



if __name__ == '__main__':
    main()




'''

code modified from

~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_high_fsc_consistency_with_ground_truth/table_excel_export.py

'''

