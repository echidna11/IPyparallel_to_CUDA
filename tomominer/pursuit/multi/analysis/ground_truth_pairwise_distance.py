#!/usr/bin/env python



# at a certain iteration, after dimension reduction, calculate the pairwise distance between subtomograms, and plot a distance matrix, where subtomograms are ordered according to ground truth

import os
import sys
from collections import defaultdict
import pickle
import json
import scipy.spatial.distance as SPD
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt



# parameters: dt: data_json with ground truth
def calculate_and_plot_dist_mat(order, root_dir, pass_i):

    pass_dir_last = os.path.join(root_dir, 'pass_%03d' % (pass_i-1))
    pass_dir = os.path.join(root_dir, 'pass_%03d' % (pass_i))

    print '# load dimension reduction results, and calculate distance matrix'
    with open(os.path.join(pass_dir, 'dimension_reduction.pickle')) as f:       dr = pickle.load(f)

    m = SPD.squareform(   SPD.pdist(dr['cfp_re']['red'])   )

    print '# load data_json file of last pass, which has an order corresponding to the dimension reduction results'
    with open(os.path.join(pass_dir_last, 'data_config.json')) as f:        dj = json.load(f)



    print '# reorder the distance matrix according to ground truth and abundance'
    ind = [None]*len(dj)
    for i,d in enumerate(dj):
        if d['subtomogram'] not in order:   continue
        ind[order[d['subtomogram']]] = i
    ind = [_ for _ in ind if _ is not None]


    m = m[ind,:][:,ind]         # matrix reordering of both columns and rows


    print '# save plots'
    out_dir = os.path.join(pass_dir, 'ground_truth_pairwise_distance')
    if not os.path.exists(out_dir) :            os.mkdir(out_dir)
    print out_dir

    pp = PdfPages(os.path.join(out_dir, 'dist_mat.pdf'))

    plt.imshow(m, cmap='Greys_r')
    #plt.pcolor(m)
    
    pp.savefig()
    pp.close()







if __name__ == '__main__':
    
    data_truth_file = sys.argv[1]

    with open(data_truth_file) as f:      dt = json.load(f)     # load ground truth

    print '# calculate the abundance of ground truth'
    dl = defaultdict(list)      # a dictionary where keys are pdb_ids and content are list of subtomograms
    for d in dt:
        if d['pdb_id'] is None:    continue        # we ignore junk subtomograms
        dl[d['pdb_id']].append(d['subtomogram'])

    pdb_ids = sorted( dl.keys(), key=lambda _: (-len(dl[_])) )        # order pdb_ids according to abundance
    for pid in pdb_ids:     print pid, len(dl[pid]), '\t',

    print '# assign a order to each subtomogram'
    order = {}
    ind_t = 0
    for pid in pdb_ids:
        for k in dl[pid]:
            order[k] = ind_t
            ind_t += 1


    import path_util as PU
    selected_folders = PU.last_pass_dirs()


    for max_t in selected_folders:

        calculate_and_plot_dist_mat(order=order, root_dir=max_t['root_path'], pass_i=max_t['pass_i'])






