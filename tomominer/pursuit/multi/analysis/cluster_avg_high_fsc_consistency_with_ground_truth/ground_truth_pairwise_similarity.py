#!/usr/bin/env python

'''

calculate pairwise structural consistancy between ground truth complexes, in terms of FSC at 0.5 cutoff


~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_high_fsc_consistency_with_ground_truth/ground_truth_pairwise_similarity.py

'''


import os, json, copy, csv, pickle


def export_table(fs, out_file):


    pdb_ids = {}
    for l1 in fs[0]:
        v2k = fs[0][l1]['align']['v2_key']
        assert  l1 == v2k['cluster_label']
        pdb_ids[l1] = v2k['pdb_id']

    lbls = sorted(list(pdb_ids.keys()))

    with open(out_file, 'w') as f:
        cw = csv.writer(f)
        
        r = ['']
        r.extend([pdb_ids[_] for _ in lbls])
        cw.writerow(r)

        for l0 in lbls:
            r = [pdb_ids[l0]]

            for l1 in lbls:
                r.append('%.2f'%fs[l0][l1]['resolution'])
            cw.writerow(r)








def main():
    with open('ground_truth_pairwise_similarity__op.json') as f:     op = json.load(f)

    with open(op['ground truth']) as f:       true_avg = json.load(f)

    v1ks = true_avg
    for r in v1ks:  
        r['subtomogram'] = r['ground_truth']
        r['cluster_id'] = int(r['cluster_label'])


    if not os.path.isfile(op['pairwise align file out']):
        align_op = {'with_missing_wedge':True, 'L':36}
        from tomominer.pursuit.multi.util import pairwise_align_keys__multiprocessing
        pa = pairwise_align_keys__multiprocessing(self=None, v1ks=v1ks, v2ks=copy.deepcopy(v1ks), op=align_op)
        with open(op['pairwise align file out'], 'wb') as f:       pickle.dump(pa, f, protocol=-1)

    else:
        print 'load existing', op['pairwise align file out']
        with open(op['pairwise align file out'], 'rb') as f:       pa = pickle.load(f)


    if not os.path.isfile(op['fsc stat file out']):
        from tomominer.pursuit.multi.analysis.cluster_avg_high_fsc_consistency_with_ground_truth.cluster_avg_high_fsc_consistency_with_ground_truth import fsc_stat
        fs = fsc_stat(pa, voxel_spacing=op['voxel_spacing'])
        with open(op['fsc stat file out'], 'wb') as f:       pickle.dump(fs, f, protocol=-1)
    else:
        print 'load existing', op['fsc stat file out']
        with open(op['fsc stat file out'], 'rb') as f:       fs = pickle.load(f)


    export_table(fs=fs, out_file=op['table out'])


if __name__ == '__main__':
    main()






'''

# R code for drawing hierarchical cluster
d0 = read.csv('table--out.csv', header=TRUE)

d = d0[,-1]
row.names(d) = d0[,1]

h = hclust(as.dist(d))


pdf('table--out--hcluster.pdf', width=8, height=6, paper='special')
plot(h, main='Cluster Dendrogram', xlab='PDB ID', ylab='Dissimilarity (nm)')
dev.off()

'''

