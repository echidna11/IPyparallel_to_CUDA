#!/usr/bin/env python

'''
over iterations of a single MPP run, scatter plot of SFSC scores of generated patterns and selected patterns
~/ln/tomominer/tomominer/pursuit/multi/analysis/iteration/sfsc_score_scatter_plot.py
'''


import sys, json, os
import cPickle as pickle

def collect_stat(fs, pass_i_current=None):

    ps = fs['passes']
    ps = {int(_):ps[_] for _ in ps}
    if pass_i_current is None:    pass_i_current = fs['pass_i_current']
    pst = ps[pass_i_current]

    # selected patterns of different iterations
    with open(pst['fsc_stat_file']) as f:        fsc_sel = json.load(f)
    fsc_sel = {int(_):fsc_sel[_] for _ in fsc_sel}
    for pass_i in fsc_sel:      fsc_sel[pass_i] = {int(_):{'fsc':fsc_sel[pass_i][_]['fsc'], 'cluster_mode':fsc_sel[pass_i][_]['cluster_mode'], 'is_specific':(fsc_sel[pass_i][_]['is_specific'] is None)} for _ in fsc_sel[pass_i]}
    
    # generated patterns for different iterations
    fsc_gen = {}
    for pass_i in range(pass_i_current):
        print '\r', pass_i,             ;       sys.stdout.flush()
        with open(ps[pass_i]['cluster_info_file'], 'rb') as f:      cit = pickle.load(f)
        fsc_gen[pass_i] = {}
        for c in cit:       fsc_gen[pass_i][c] = {'fsc':cit[c]['fsc'].sum(), 'cluster_mode':cit[c]['cluster_mode']}

    return {'gen':fsc_gen, 'sel':fsc_sel}


def plot(fsc, figure_size, out_file):

    gen = fsc['gen']
    sel = fsc['sel']

    import matplotlib.pyplot as PLT
    fig = PLT.figure(figsize=tuple(figure_size))
    ax1 = fig.add_subplot(111)

    point_size = 20
     
    x, y = fsc_collect_mode(gen, cluster_mode='kmeans')
    ax1.scatter(x, y, s=point_size, edgecolors='lightgrey', facecolors='none', marker="D", label='Generated kmeans')

    x, y = fsc_collect_mode(gen, cluster_mode='kmeans-adaptave')
    ax1.scatter(x, y, s=point_size, edgecolors='lightgrey', facecolors='none', marker="s", label='Generated kmeans adaptive')

    x, y = fsc_collect_mode(gen, cluster_mode='sequential')
    ax1.scatter(x, y, s=point_size, edgecolors='lightgrey', facecolors='none', marker="^", label='Generated sequential expansion')

    point_size = 40
    x, y = fsc_collect_mode(sel, cluster_mode='kmeans')
    ax1.scatter(x, y, s=point_size, edgecolors='green', facecolors='none', marker="o", label='Selected kmeans')

    x, y = fsc_collect_mode(sel, cluster_mode='kmeans-adaptave')
    ax1.scatter(x, y, s=point_size, edgecolors='red', facecolors='none', marker="o", label='Selected kmeans adaptive')

    x, y = fsc_collect_mode(sel, cluster_mode='sequential')
    ax1.scatter(x, y, s=point_size, edgecolors='blue', facecolors='none', marker="o", label='Selected sequential expansion')
 
    point_size = 40
    x, y = fsc_collect_mode(sel, is_specific=False)
    ax1.scatter(x, y, s=point_size, color='black', marker="x", label='Redundant')

    PLT.xlim(-1, max(list(gen.keys()))+2)
    PLT.legend(loc='lower right')
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('SFSC')

    fig.savefig(out_file)
    #PLT.show()

def fsc_collect_mode(f, cluster_mode=None, is_specific=None):
    
    x = []
    y = []
    for p in f:
        for c in f[p]:
            ft = f[p][c]
            if cluster_mode is not None:
                if ft['cluster_mode'] != cluster_mode:      continue
            
            if is_specific is not None:
                if ft['is_specific'] != is_specific:        continue

            x.append(p)
            y.append(ft['fsc'])

    return (x, y)


def main():
    with open('sfsc_score_scatter_plot__op.json') as f:
        op = json.load(f)
    with open(op['rep_file']) as f:
        rep_file = json.load(f)
    
    op["rep"] = rep_file['rep']
    for i in op["rep"]:
        file_stat = os.path.abspath(os.path.join("../../out/group", str(op["group"]).zfill(4), "rep", str(i).zfill(4), "out/file_stat.json"))
        with open(file_stat) as f:
            fs = json.load(f)
        fsc = collect_stat(fs, pass_i_current=(op['pass_i_current'] if 'pass_i_current' in op else None))
        plot(fsc, figure_size = op['figure size'], out_file = (op['figure out']+str(i)+".eps") )

if __name__ == '__main__':
    main()
