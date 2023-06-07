#!/usr/bin/env python

'''
given subtomograms extracted from tomogram, and ground truth of the subtomogram pdb_id and rotation
we randomly select pairs of subtomograms of same structure, then rotate back according to ground truth so that the two subtomograms are aligned.

then calculate SNR through Pearson correlation according to Joachim Frank's estimation method
IMPORTANT: make sure following are similiar to experimenal data:  1) the subtomogram's size, 2) CTF and MTF, 3) tilt angle range, 4) the macromolecular complex of interest

~/ln/tomominer/tomominer/simulation/tomogram/instance_pair_snr.py

'''
import json, itertools
import tomominer.io.file as IF
import tomominer.geometry.ang_loc as GAL
import tomominer.geometry.rotate as GR
import numpy as N
import scipy.stats as SS
import tomominer.image.vol.util as CV

if __name__ == "__main__":

    with open('instance_pair_snr__op.json') as f:   op = json.load(f)
    with open(op['model_file']) as f: models = json.load(f)
    with open(op['model_stat_file']) as f: model_stat = json.load(f)
    tomogram_ids = op['selected_tomogram_ids']
    tom_box = model_stat["packing"]["param"]["box"]
    op['tomogram_size'] = [tom_box['x'], tom_box['y'], tom_box['z']]
    subtomogram_list = {}
    subtomograms = {}
    for pid in op['pdb_id']:
        subtomogram_list[pid] = []
        subtomograms[pid] = []
    for i in tomogram_ids:
        i = str(i)
        for pid in op['pdb_id']:
            v_list = [_ for _ in models[i]['instances'] if _['pdb_id']==pid]
            for j,vs in enumerate(v_list):
                tmp = N.array([(vs['x'][_] - int(N.ceil(vs['radius']))) for _ in range(3)])
                if (tmp < 0).sum() > 0: continue
                tmp = N.array([(vs['x'][_] + int(N.ceil(vs['radius']))) for _ in range(3)])
                if (tmp[0] > op['tomogram_size'][0]) or (tmp[1] > op['tomogram_size'][1]) or (tmp[2] > op['tomogram_size'][2]): continue
                subtomogram_list[pid].append(vs)
        # Now we have list of indexes that need to be extracted we will extract all the subtomograms
        tom = op['tomogram_files']
        tmp = tom.find('I')
        tom = tom[:tmp]+i+tom[tmp+1:]
        tmp = tom.find('I')
        tom = tom[:tmp]+i+tom[tmp+1:]
        print "Loading tomogram", i
        v = N.load(tom)
        siz = N.array(v.shape)
        for pid in op['pdb_id']:
            for j,vs in enumerate(subtomogram_list[pid]):
                vc = CV.cut_from_whole_map(whole_map=v, c=vs['x'], siz=N.array(op['subtomogram_size']))
                if vc is None:  continue
                subtomograms[pid].append(vc)
                # we will calculate snr for each tomogram separately and for each pdb_id. Juast make sure we have 10,000 pairs, i.e. atleast 142 subtomograms
    for pid in op['pdb_id']:
        if len(subtomograms[pid]) < op['len']:
            print "Not enough subtomograms for pdb_id", pid
            continue
        print 'sampling', op['pair_num'], 'pairs from', len(subtomograms[pid]), 'aligned subtomograms'
        out_file = op['out_file'] + pid + ".txt"
        combination_list = list(itertools.combinations(list(range(len(subtomograms[pid]))), 2))
        pairs = N.random.choice(list(range(len(combination_list))), op['pair_num'])
        s_list = []
        for j in range(op['pair_num']):
            j0 = combination_list[pairs[j]][0]
            j1 = combination_list[pairs[j]][1]
            c, _ = SS.pearsonr(subtomograms[pid][j0].flatten(), subtomograms[pid][j1].flatten())        # pearson correlation between two aligned subtomograms
            if c!=1:
                s = c / (1.0 - c)         # SNR according to Equation 3.61 of book     frank Three-dimensional Electron Microscopy of Macromolecular Assemblies
                s_list.append(s)
        outstr = "Average: " + str(N.mean(s_list)) + "\n"
        print(outstr)
        for s in (s_list):
            outstr += str(s) + "\n"
        f = open(out_file, "w")
        f.write(outstr)
        f.close()

