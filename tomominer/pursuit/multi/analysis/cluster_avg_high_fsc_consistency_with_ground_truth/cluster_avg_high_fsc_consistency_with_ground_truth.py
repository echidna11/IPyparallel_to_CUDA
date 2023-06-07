#!/usr/bin/env python





'''
calculate contegincy table for each cluster
find its largest overlaping true cluster, 
then see if the true cluster is fully covered by this cluster

see if the predicted cluster also can correctly seperate other clusters

perform pairwise alignment and calculate structual consistency, in terms of cross correlation, and fourier shell correlation, list together with contegincy table

calculate and report overall cluster statistics using scikit-learn

~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_high_fsc_consistency_with_ground_truth/cluster_avg_high_fsc_consistency_with_ground_truth.py

'''



import sys, os, shutil, copy, json, pickle, csv
import numpy as N

from collections import defaultdict

import tomominer.io.file as IV
import tomominer.statistics.vol as SV
import tomominer.core as core



def fsc_stat(pa, voxel_spacing, out_dir=None, count_cluster=None):
    print 'fsc_stat()', out_dir

    if out_dir is not None:
        if not os.path.isdir(out_dir):      os.makedirs(out_dir)

    count = 0

    record = defaultdict(dict)
    for a in pa:
        v1k = a['v1_key']['subtomogram']
        v2k = a['v2_key']['subtomogram']

        print v1k
        print v2k
        v1 = IV.get_mrc(v1k)
        v2 = IV.get_mrc(v2k)

        # my code for test
        print "Before rotation of 2nd subtomogram:", v1.shape, v2.shape
        v2r = core.rotate_vol_pad_mean(v2, a['angle'], a['loc'])
        print "After rotation of 2nd subtomogram:", v1.shape, v2r.shape

        id1 = a['v1_key']['cluster_id']
        id2 = a['v2_key']['cluster_id']


        if out_dir is not None: 
            if (count_cluster is None) or (count_cluster[id1][id2] > 0):
                v1ko = os.path.join(out_dir, '%03d.mrc'%(id1,))
                v2ko = os.path.join(out_dir, '%03d-%03d.mrc'%(id1, id2))
                IV.put_mrc( v2r, v2ko )         # if there are not shared members between a true and predicted cluster, we do not write aligned volumes, in order to same storage
                if not os.path.isfile(v1ko):    IV.put_mrc( v1, v1ko )         # also save ground truth for reference

        resolution = SV.resolution(v1, v2r, voxel_spacing=voxel_spacing, cutoff=0.5)

        assert id2 not in record[id1]
        record[id1][id2] = {}

        record[id1][id2]['v1k'] = v1k
        record[id1][id2]['v2k'] = v2k
        record[id1][id2]['align'] = a
        record[id1][id2]['resolution'] = resolution
        record[id1][id2]['fsc'] = SV.fsc(v1, v2r)
        record[id1][id2]['cor'] = N.corrcoef(v1.flatten(), v2r.flatten())[0,1]

        print '\r', count, float(count) / len(pa), '                        ',           ;       sys.stdout.flush()
        count += 1


    return record



def stat(op, out_dir):
    with open(op['true_member']) as f:          true_member = json.load(f)
    with open(op['true_avg']) as f:             true_avg = json.load(f)

    pdb_lbl = {_['pdb_id']:_['cluster_label'] for _ in true_avg}

    true_lbl = {}
    for _ in true_member:
        # pdb_id should be checked before cluster_label because in case of analysis of ground truth both are present and we want pdb_id to be picked up in that case.
        if 'pdb_id' in _:
            if _['pdb_id'] is None:     continue        # this means that the extracted subtomogram is junk and does not contain a valid complex
            true_lbl[ _['subtomogram'] ] = pdb_lbl[_['pdb_id']]
            continue
        
        if 'cluster_label' in _:
            true_lbl[ _['subtomogram'] ] = int(_['cluster_label'])
            continue
        
        continue        # this means that the extracted subtomogram is junk and does not contain a valid complex
        #raise Exception('cannot identify true label')



    with open(op['pred_member_file']) as f:  pred_member = json.load(f)         # the data json file in pass folder and with best aligned templates
    pred_lbl_match = {_['subtomogram'] : _['template']['id'] for _ in pred_member if (('template' in _) and ('id' in _['template']))}     # predicted labels (template id)  defined by best template match after alignment


    #------------------------------------
    # get selected templates
    matched_template_id__set = set( _['template']['id'] for _ in pred_member if (('template' in _) and ('id' in _['template'])) )       # the collection of mached template ids

    with open(op['cluster_average_select_file']) as f:  cas = pickle.load(f)

    selected_templates = cas['selected_templates']
    for _ in selected_templates.keys():
        if _ not in matched_template_id__set:
            del selected_templates[_]       #  we restrict only to those selcted templates that has best matches

    selected_templates_inv = {selected_templates[_]['subtomogram']:int(_) for _ in selected_templates}      # map from selected template's subtomogram key to its id

    # in order to save storage, clean up cas to keep only needed entries
    unused_templates = [_ for _ in cas['tk_info'] if _ not in selected_templates_inv]
    for _ in unused_templates:        del cas['tk_info'][_]




    # -----------------------------------
    # get pred_lbl_cluster
    cluster_info = {}
    for pass_i in op['file_stat']['passes']:
        if pass_i > op['selected_pass_i']:      continue
        with open(op['file_stat']['passes'][pass_i]['cluster_info_file'], 'rb') as f:    cluster_info[pass_i] = pickle.load(f)

    with open(op['cluster_info_stat_file'], 'rb') as f:    cluster_info_stat = pickle.load(f)
    for pass_i in cluster_info:
        for c in cluster_info[pass_i]:
            cluster_info[pass_i][c]['is_specific'] = cluster_info_stat[pass_i][c]['is_specific']


    pred_lbl_cluster = {}           # predicted labels defined by clusters
    for pass_i in cluster_info:
        for c in cluster_info[pass_i]:
            if not cluster_info[pass_i][c]['template_key']['subtomogram'] in selected_templates_inv:    continue
            
            assert      (cluster_info[pass_i][c]['is_specific'] is not True) and (cluster_info[pass_i][c]['is_specific'] is not False)      # we do not work on old version of boolean assignment any more.

            if cluster_info[pass_i][c]['is_specific'] is not None:  continue            # we ignore those non-specific clusters!!!      Note in new format, is_specific==(a record)

            clus_lbl_t = selected_templates_inv[cluster_info[pass_i][c]['template_key']['subtomogram']]

            for s in cluster_info[pass_i][c]['data_json']:            pred_lbl_cluster[s['subtomogram']] = clus_lbl_t


    true_label_max = max(max([pdb_lbl[_] for _ in pdb_lbl]), max([true_lbl[_] for _ in true_lbl]))
    pred_label_max = max([pred_lbl_match[_] for _ in pred_lbl_match])        # number of subtomograms in a true cluster that best match the average of a predicted cluster


    # ---------------------------------------
    # counting contegency tables
    count_match = N.zeros([true_label_max+1, pred_label_max+1], dtype=N.int)        # number of subtomograms in a true cluster that best match the average of a predicted cluster
    for s in pred_lbl_match:
        if s not in true_lbl:   continue
        count_match[true_lbl[s], pred_lbl_match[s]] += 1

    count_cluster = N.zeros(count_match.shape, dtype=N.int)         # number of subtomograms that belongs to a true and a predicted cluster at same time
    for s in pred_lbl_cluster:
        if s not in true_lbl:       continue

        count_cluster[true_lbl[s], pred_lbl_cluster[s]] += 1



    # -----------------------------
    # align and calculate FSC
    v1ks_v = true_avg

    v1ks = []
    for r in v1ks_v:  
        r['subtomogram'] = r['ground_truth']
        r['cluster_id'] = int(r['cluster_label'])
        v1ks.append(r)

    pid = {_['cluster_id']:str(_['pdb_id']) for _ in v1ks }


    v2ks_d = selected_templates

    v2ks = []
    for c in v2ks_d:
        r = v2ks_d[c]
        r['cluster_id'] = c
        v2ks.append(r)


    align_op = {'with_missing_wedge':True, 'L':36}
    import tomominer.pursuit.multi.util as CU
    pa = CU.pairwise_align_keys__multiprocessing(self=None, v1ks=v1ks, v2ks=v2ks, op=align_op)

    fs = fsc_stat(pa, voxel_spacing=op['voxel_spacing'], out_dir=(os.path.join(out_dir, 'aligned-averages') if op['output aligned averages'] else None), count_cluster=count_cluster)

    # ------------------------------------------------
    clus_stat = {}
    for cp in range(pred_label_max + 1):
        if cp not in selected_templates:    continue
        clus_stat[cp] = {'id':cp, 'size':len([_ for _ in pred_lbl_cluster if pred_lbl_cluster[_] == cp]), 'fsc':cas['tk_info'][selected_templates[cp]['subtomogram']]['fsc'].sum(), 'pass_i':cas['tk_info'][selected_templates[cp]['subtomogram']]['pass_i'], 'cluster':cas['tk_info'][selected_templates[cp]['subtomogram']]['cluster']}


    return {'selected_templates':selected_templates, 'pred_lbl_cluster':pred_lbl_cluster, 'pid':pid, 'pdb_lbl':pdb_lbl, 'true_lbl':true_lbl, 'pred_lbl_match':pred_lbl_match, 'pred_lbl_cluster':pred_lbl_cluster, 'fs':fs, 'clus_stat':clus_stat, 'true_label_max':true_label_max, 'pred_label_max':pred_label_max, 'count_match':count_match, 'count_cluster':count_cluster}






def stat__count(st):
    true_lbl = st['true_lbl']
    pred_lbl_match = st['pred_lbl_match']
    pred_lbl_cluster = st['pred_lbl_cluster']

    st['count_match'] = count_match
    st['count_cluster'] = count_cluster



# form a table of statistics
def stat_table(st):

    pid = st['pid']
    fs = st['fs']
    clus_stat = st['clus_stat']


    template_stat = {}
    for ct in range(st['true_label_max'] + 1):
        if ct not in pid:   continue
        template_stat[ct] = {'id':ct, 'pid':pid[ct], 'instance_num':st['count_match'][ct, :].sum()}


    # a table of statistics
    # IMPORTANT: because generated using MATLAB, true class label starts from 1
    tab = defaultdict(dict)
    for cp in range(st['count_match'].shape[1]):
        if cp not in clus_stat: continue

        for ct in range(st['count_match'].shape[0]):
            if ct not in template_stat: continue

            tab[ct][cp] = {'count_match':st['count_match'][ct][cp], 'count_cluster':st['count_cluster'][ct][cp], 'resolution':fs[ct][cp]['resolution'], 'cor':fs[ct][cp]['cor']}

    return {'clus_stat':clus_stat, 'template_stat':template_stat, 'tab_stat':tab}



# prepare an ideal table for comparison purpose
def stat_table_ideal(tab):
    tab = copy.deepcopy(tab)

    template_stat = tab['template_stat']
    clus_stat = {_:{'size':template_stat[_]['instance_num']} for _ in template_stat}

    tab_stat = defaultdict(dict)
    for cp in template_stat:
        for ct in clus_stat:
            if cp == ct:
                tab_stat[ct][cp] = {'count_cluster':template_stat[ct]['instance_num'], 'resolution':tab['resolution_min']}
            else:
                tab_stat[ct][cp] = {'count_cluster':0, 'resolution':tab['resolution_max']}

    return  {'template_stat':template_stat, 'clus_stat':clus_stat, 'tab_stat':tab_stat, 'resolution_min':tab['resolution_min'], 'resolution_max':tab['resolution_max']}






# use Hungarian algorithm to arrange the table
def stat_table__arrange(tab):

    template_stat = tab['template_stat']
    clus_stat = tab['clus_stat']
    tab_stat = tab['tab_stat']

    cts = template_stat.keys()
    cps = clus_stat.keys()
    
    cnt = N.zeros(  (len(cts), len(cps)), dtype=N.int  )
    for ct_i, ct in enumerate(cts):
        for cp_i, cp in enumerate(cps):
            cnt[ct_i, cp_i] = tab_stat[ct][cp]['count_cluster']

    ai = stat_table__arrange__ind(cnt=cnt, cts=cts, cps=cps)

    tab['arrange'] = {'cts':ai['cts'], 'cps':ai['cps'], 'cnt':ai['cnt_v']}



def stat_table__arrange__ind(cnt, cts, cps):
    # use Hungarian algorithm to match the counts
    cnt_t = cnt.max() - cnt

    from munkres import Munkres
    m = Munkres()
    ind = m.compute(cnt_t.tolist())
    cnt_v = [None] * len(ind)       # vector of corresponding entries in cnt
    ct_ind = [None] * len(ind)
    cp_ind = [None] * len(ind)


    # order the matched counts
    for i, (ct_i, cp_i) in enumerate(ind):
        cnt_v[i] = cnt[ct_i, cp_i]
        ct_ind[i] = ct_i
        cp_ind[i] = cp_i

    cnt_v = N.array(cnt_v, dtype=N.int)
    i = N.argsort(-cnt_v)

    cnt_v_s = [cnt_v[_] for _ in i]

    cts_s = [cts[ct_ind[_]] for _ in i]
    cps_s = [cps[cp_ind[_]] for _ in i]


    # add the rest items if any
    cts_s.extend(set(cts).difference(cts_s))
    cps_s.extend(set(cps).difference(cps_s))

    return {'cts':cts_s, 'cps':cps_s, 'cnt_v':cnt_v_s}





def stat_table_out(tab, out_dir):

    clus_stat = tab['clus_stat']
    template_stat = tab['template_stat']
    tab_stat = tab['tab_stat']


    

    # print a table of statistics
    t = []
    r = ['', '', '']
    for cp_i, cp in enumerate(tab['arrange']['cps']):
        r.append('%d, %d, %.2f, %d, %d'%(clus_stat[cp]['id'], clus_stat[cp]['size'], clus_stat[cp]['fsc'], clus_stat[cp]['pass_i'], clus_stat[cp]['cluster'] ))        # output cluster id, cluster size, cluster average's FSC

    t.append(r)

    for ct_i, ct in enumerate(tab['arrange']['cts']):
        r = [template_stat[ct]['id'], template_stat[ct]['pid'], template_stat[ct]['instance_num']]

        for cp_i, cp in enumerate(tab['arrange']['cps']):
            #r.append( '%d, %d, %.1f, %.2f'%(tab_stat[ct][cp]['count_match'], tab_stat[ct][cp]['count_cluster'],tab_stat[ct][cp]['resolution'], tab_stat[ct][cp]['cor']) )
            r.append( '%d, %.1f'%(tab_stat[ct][cp]['count_cluster'], tab_stat[ct][cp]['resolution']) )

        t.append(r)


    out_file = os.path.join(out_dir, 'stat.csv')
    with open(out_file, 'w') as f:
        cw = csv.writer(f, quotechar='\"')
        cw.writerows(t)

    print 'written:', '\t', out_file



# plot the count of memberships, color the bars according to FSC!!
def stat_table__plot_count(tab, out_file, label_font_size=24):
    clus_stat = tab['clus_stat']
    template_stat = tab['template_stat']
    tab_stat = tab['tab_stat']

    cts = tab['arrange']['cts']
    cps = tab['arrange']['cps']


    import matplotlib
    matplotlib.use('Qt4Agg')

    import matplotlib.cm as CM



    # coding the resolution using colors

    resolution_tab = N.zeros(   (len(cts), len(cps))  )
    for ct_i, ct in enumerate(cts):
        for cp_i, cp in enumerate(cps):
            resolution_tab[ct_i, cp_i] = tab_stat[ct][cp]['resolution']
    
    rainbow_color_num = 1000
    rainbow_color_lsp = N.linspace(0, 1, (rainbow_color_num+1))
    rainbow_color = CM.rainbow(rainbow_color_lsp)
    resolution_tab_norm = N.ceil(rainbow_color_num * (resolution_tab - tab['resolution_min']) / (tab['resolution_max'] - tab['resolution_min'])).astype(N.int)
    resolution_tab_norm[resolution_tab_norm < 0] = 0        # somehow the zero value gives purple color, which is strange
    resolution_tab_norm[resolution_tab_norm > rainbow_color_num] = rainbow_color_num


    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt



    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(111, projection='3d')
    #ax.pbaspect = [1.0, 2.0, 5]


    max_cnt = 0

    # plot cluster sizes
    xs = []
    ys = []
    for cp_i, cp in enumerate(cps):
        xs.append(cp_i)
        ys.append(clus_stat[cp]['size'])
    ax.bar(xs, ys, zs=-1, zdir='y', color='y', alpha=0.4, linewidth=0)          #   see https://github.com/matplotlib/matplotlib/blob/master/lib/mpl_toolkits/mplot3d/axes3d.py

    ax.xaxis.set_ticks(N.array(range(len(cps))))
    ax.xaxis.set_ticklabels([_ for _ in cps], size=label_font_size, verticalalignment='center', horizontalalignment='right')
    ax.xaxis._axinfo['ticklabel']['space_factor'] = 0.35

    max_cnt = max(max_cnt, max(ys))

    # plot ground truth instance number
    zs = []
    ys = []
    for ct_i, ct in enumerate(cts):
        ys.append(template_stat[ct]['instance_num'])
        zs.append(ct_i)
    ax.bar(zs, ys, zs=[-1]*len(ys), zdir='x', color='y', alpha=0.4, linewidth=0)
    max_cnt = max(max_cnt, max(ys))

    ax.yaxis.set_ticks(N.array(range(len(cts))))
    y_lbl = [template_stat[_]['pid'] for _ in cts]
    ax.yaxis.set_ticklabels(y_lbl, size='small', rotation=-35, verticalalignment='baseline', horizontalalignment='left')
    ax.yaxis._axinfo['ticklabel']['space_factor'] = 0.2


    '''
    print str(ax.yaxis.get_major_ticks()[0].label)
    for t in ax.yaxis.get_major_ticks():       t.set_pad(0)
    '''



    # plot match counts
    for ct_i, ct in enumerate(cts):
        xs = []
        ys = []

        for cp_i, cp in enumerate(cps):
            xs.append(cp_i)
            ys.append(tab_stat[ct][cp]['count_cluster'])

        # see http://matplotlib.org/mpl_toolkits/mplot3d/api.html
        ax.bar(xs, ys, zs=ct_i, zdir='y', bottom=-30, color=rainbow_color[resolution_tab_norm[ct_i, :].flatten()], alpha=0.4, linewidth=0, align='center')      # set bottom to a negative value so that the bar can stand on the xy plane which actually has a negative axis below z=0

    #ax.set_xlim(-1, len(cps), auto=False)
    #ax.set_ylim(-1, len(cts), auto=False)

    z_ticks = ax.get_zticks()
    z_ticks = [_ for _ in z_ticks if _ >= 0]        # ignore negative values
    ax.set_zticks(z_ticks)
    ax.set_zticklabels([str(int(_)) for _ in z_ticks], size=label_font_size)


    ax.set_xlabel('Pattern id', size=label_font_size)
    ax.xaxis._axinfo['label']['space_factor'] = 1         # use this to adjust the distance of label to axis
    
    ax.set_ylabel('PDB id', size=label_font_size)
    ax.yaxis._axinfo['label']['space_factor'] = 1         # use this to adjust the distance of label to axis
    
    ax.set_zlabel('Member count', size=label_font_size)


    ax.invert_yaxis()
    #ax.set_aspect(1)
    #ax.axis('equal')
    max_clus_template_num = max(len(cps), len(cts))
    #ax.auto_scale_xyz([0, max_clus_template_num], [0, max_clus_template_num], [0, max_cnt+100])

    ax.set_zlim3d([0, max_cnt+100])

    # set transparent background
    ax.patch.set_facecolor('None')
    ax.patch.set_visible(False)

    plt.draw()

    if True:
        plt.savefig(out_file, bbox_inches='tight', transparent=True)
    else:
        print 'select a suitable aspect and size, then save the figure as PDF file, then close'
        plt.show()


    plt.close('all')



'''

todo: plot multiple bars
see
http://matplotlib.org/examples/api/barchart_demo.html

'''

def stat_table__plot_colorbar(tab, out_dir, label_font_size=24):

    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.cm as CM


    # plot and save color bar
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.1, 0.8])

    
    cbar = matplotlib.colorbar.ColorbarBase(ax, cmap=CM.rainbow, norm=matplotlib.colors.Normalize(vmin=tab['resolution_min'], vmax=tab['resolution_max']), orientation='vertical')
    cbar_ticks = [round(_*10.0)/10.0 for _ in N.linspace(tab['resolution_min'], tab['resolution_max'], 10)]
    cbar.set_ticks(cbar_ticks)

    cbar_ticklabels = ['1/%.1f'%(_,) for _ in cbar_ticks]
    cbar_ticklabels[0] = '>' + cbar_ticklabels[0]
    cbar_ticklabels[-1] = '<' + cbar_ticklabels[-1]
    cbar.set_ticklabels(cbar_ticklabels)

    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(label_font_size)


    cbar.set_label('Frequency (FSC cutoff=0.5, unit 1/nm)', size=label_font_size)


    # set transparent background
    ax.patch.set_facecolor('None')
    ax.patch.set_visible(False)


    plt.savefig(os.path.join(out_dir, 'cluster_count_colorbar.pdf'), bbox_inches='tight')

    #plt.show()

    plt.close('all')






# for each matched template - cluster center pair, rotate the cluster center according to alignment, and save them
def save_best_matched_averages(cts, cps, fs, out_dir, density_positive=True):

    print 'save_best_matched_averages()'

    if not os.path.isdir(out_dir):      os.makedirs(out_dir)

    for i in range(min(len(cts), len(cps))):
        ct = cts[i]
        cp = cps[i]

        v1 = IV.get_mrc(fs[ct][cp]['v1k'])
        if not density_positive:    v1 = -v1
        IV.put_mrc(v1, os.path.join(out_dir, '%03d--%03d.mrc'%(i, ct,)), overwrite=True)        # save ground truth

        v2 = IV.get_mrc(fs[ct][cp]['v2k'])
        if not density_positive:    v2 = -v2
        v2r = core.rotate_vol_pad_mean(v2, fs[ct][cp]['align']['angle'], fs[ct][cp]['align']['loc'])

        IV.put_mrc(v2r, os.path.join(out_dir, '%03d--%03d-%03d.mrc'%(i, ct,cp)), overwrite=True)    # save predicted cluster average




def main(op_org):

    # open file_stat generated in same folder as pass directory
    with open(op_org['file_stat_file']) as f:       file_stat = json.load(f)
    # read all info about each pass and store it under each pass number where pass number are keys of file_stat
    file_stat['passes'] = {int(_):file_stat['passes'][_] for _ in file_stat['passes']}

    # read paths to files of particular pass or last pass
    if 'selected_pass_i' not in op_org:         op_org['selected_pass_i'] = file_stat['pass_i_current']
    file_stat_p = file_stat['passes'][op_org['selected_pass_i']]


    # create output directory inside particular pass and check whether directory already exist and need to create again depending on parameter re-generate-all
    out_dir = os.path.abspath(os.path.join(file_stat_p['pass_dir'], 'cluster_avg_high_fsc_consistency_with_ground_truth'))
    if op_org['re-generate-all'] and os.path.isdir(out_dir):                shutil.rmtree(out_dir)
    
    if not os.path.isdir(out_dir):     os.makedirs(out_dir)


    if not os.path.isabs(op_org['stat_file_out']):      op_org['stat_file_out'] = os.path.join(out_dir, op_org['stat_file_out'])

    if os.path.isfile(op_org['stat_file_out']):
        print 'loading existing', op_org['stat_file_out']
        with open(op_org['stat_file_out'], 'rb') as f:     st = pickle.load(f)
    else:
        op = copy.deepcopy(op_org)
        op['file_stat'] = file_stat
        op['pred_member_file'] = file_stat_p['data_json_file']
        op['cluster_average_select_file'] = file_stat_p['cluster_average_select_file']
        op['cluster_info_stat_file'] = file_stat_p['cluster_info_stat_file']
        
        # function-call "stat", arguements are out_dir and op which contains parameters from input json file and information of pass we are processing
        st = stat(op=op, out_dir=out_dir)
        with open(op_org['stat_file_out'], 'wb') as f:    pickle.dump(st, f, protocol=-1)


    if not os.path.isabs(op_org['tab_file_out']):      op_org['tab_file_out'] = os.path.join(out_dir, op_org['tab_file_out'])

    if os.path.isfile(op_org['tab_file_out']):
        print 'loading', op_org['tab_file_out']
        with open(op_org['tab_file_out'], 'rb') as f:         tab = pickle.load(f)
    else:

        tab = stat_table(st)

        # set lower and upper bound of resolutions for generating color maps
        tab['resolution_min'] = op_org['resolution_ideal']
        if False:
            v = IV.get_mrc(st['fs'][0][0]['v1k'])
            tab['resolution_max'] = N.max(v.shape)
            del v
        elif False:
            tab['resolution_max'] = N.max([N.max([tab['tab_stat'][_0][_1]['resolution'] for _1 in tab['tab_stat'][_0]]) for _0 in tab['tab_stat']])

            tab['resolution_max'] *= op_org['resolution_max_factor']

        else:

            tab['resolution_max'] = tab['resolution_min'] * op_org['resolution_max_factor']

        assert      tab['resolution_max'] > tab['resolution_min']
        print 'resolution_max', tab['resolution_max']

        stat_table__arrange(tab)

        with open(op_org['tab_file_out'], 'wb') as f:        pickle.dump(tab, f, protocol=-1)


    if op_org['output best aligned averages']:      save_best_matched_averages(cts=tab['arrange']['cts'], cps=tab['arrange']['cps'], fs=st['fs'], out_dir=os.path.join(out_dir, 'best-matched-averages'), density_positive=op_org['density_positive'])
    stat_table_out(tab, out_dir=out_dir)
    stat_table__plot_count(tab, out_file=os.path.join(out_dir, 'cluster_count.pdf'))
    stat_table__plot_colorbar(tab, out_dir=out_dir)

    if op_org['plot ideal']:
        tab_ideal = stat_table_ideal(tab)
        stat_table__arrange(tab_ideal)
        stat_table__plot_count(tab_ideal, out_file=os.path.join(out_dir, 'cluster_count-ideal.pdf'))


if __name__ == '__main__':
    with open('cluster_avg_high_fsc_consistency_with_ground_truth__op.json') as f:    op = json.load(f)
    main(op)







