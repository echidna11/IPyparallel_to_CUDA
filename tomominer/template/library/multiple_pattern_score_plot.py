#!/usr/bin/env python




'''

given the alignment result against each template, this code is used for plotting the box plot of score distribution for multiple patterns, 

this code uses the output from programs in ~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/ in order to number patterns and collect subtomograms of different patterns


~/ln/tomominer/tomominer/template/library/multiple_pattern_score_plot.py

'''



import os, json, pickle, csv
import numpy as N

from tomominer.template.library.score_export import score_matrix, export_tab


def main():

    with open('multiple_pattern_score_plot__op.json') as f:      op = json.load(f)
    op['out dir'] = os.path.abspath(op['out dir'])
    if not os.path.isdir(op['out dir']):        os.makedirs(op['out dir'])


    alop_dir = os.path.dirname(os.path.abspath(op['reference multi op']))
    with open(op['reference multi op']) as f:      alop = json.load(f)
    with open(os.path.join(alop_dir, alop['data_json_file'])) as f:      dj = json.load(f)
    with open(os.path.join(alop_dir, alop['out_file'])) as f:      al = json.load(f)
    with open(os.path.join(alop_dir, alop['reference_json_file'])) as f:    ref = json.load(f)

    s = score_matrix(al)

    inds = collect_subtomogram_inds(fsc_collect=op['fsc collect'], iso_plot_stat=op['iso plot stat'], dj=dj)
    pids = [_['pid'] for _ in ref]


    for pid_i, pid in enumerate(pids):
        plot_file = os.path.join(op['out dir'], pid + '.eps')
        if not os.path.isfile(plot_file):        plot(s={_:s[inds[_], pid_i] for _ in inds}, pid=pid, out_file=plot_file)

    with open(op['wilcox test out file'], 'w') as fh:
        for pid_i, pid in enumerate(pids):
            print >>fh, 'wilcox test for', pid,
            wilcox_test(s={_:s[inds[_], pid_i] for _ in inds}, fh=fh)
            print >>fh


 


def collect_subtomogram_inds(fsc_collect, iso_plot_stat, dj):
    dji = {dj[_]['subtomogram']:_ for _ in range(len(dj))}

    with open(fsc_collect, 'rb') as f:        co = pickle.load(f)         # use this to get the list of subtomograms for each pattern, and whether pattern is generated through pursuit or refinment

    cod = {}
    for fs_fn in co:
        for c in co[fs_fn]:
            t = co[fs_fn][c]
            cod[t['template']['subtomogram']] = t

    with open(iso_plot_stat) as f:     ip = json.load(f)
    inds = {}
    for i, ipt in enumerate(ip):
        k = ipt['subtomogram']
        inds[i] = [dji[_['subtomogram']] for _ in cod[k]['dj']]

    return inds
    


def plot(s, pid, out_file):
    import rpy2.robjects as RR
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()

    # order scores and references in descending of median
    sm = [None] * len(s)
    for i in s:
        sm[i] = N.median(s[i])
    sm_i = N.argsort(-N.array(sm))

    RR.r('s = list()')
    for i in range(len(sm_i)):
        RR.r.assign('t', s[sm_i[i]])
        RR.r('s[[%d]] = t'%(i+1,))

    RR.r('names_t = vector()')
    for i in range(len(sm_i)):
        RR.r('names_t[[%d]] = \'%d\''%(i+1, sm_i[i]))

    RR.r('setEPS()')
    RR.r('postscript(\'%s\')'%(out_file,))
    RR.r('boxplot(s, names=names_t, xlab=\'Pattern ID\', ylab=\'Score\', main=\'Alignment against template of %s\')' % (pid, ))
    RR.r('dev.off()')


'''
first order the patterns in decreasing mean score order, 
then devide the pattern into two groups, 
and calculate pairwise wilcox test between patterns inside these two groups
then get maximum p-value
'''
def wilcox_test(s, fh):
    import rpy2.robjects as RR
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()

    # order scores and references in descending of median
    sm = [None] * len(s)
    for i in s:
        sm[i] = N.median(s[i])
    sm_i = N.argsort(-N.array(sm))
    sm_i = sm_i.tolist()

    for sm_i_i, sm_i_t in enumerate(sm_i):
        sm_i_i0s = sm_i[:sm_i_i]
        sm_i_i1s = sm_i[sm_i_i:]

        if not sm_i_i0s:    continue
        if not sm_i_i1s:    continue
        
        p = []
        for sm_i0 in sm_i_i0s:
            for sm_i1 in sm_i_i1s:
                RR.r.assign('s0', s[sm_i0])
                RR.r.assign('s1', s[sm_i1])
                RR.r('w = wilcox.test(x=s0, y=s1, alternative=\'greater\');' )
                p.append(RR.r('w$p.value')[0])

        print >>fh, '%d,%g\t'%(sm_i_i, max(p)),
    print >>fh



if __name__ == '__main__':
    main()




'''
related code

~/ln/tomominer/tomominer/template/library/score_export.py

~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/template_stat_table_export.py

'''

