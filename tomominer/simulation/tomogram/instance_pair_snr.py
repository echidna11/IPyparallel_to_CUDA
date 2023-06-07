#!/usr/bin/env python
'''
given subtomograms extracted from tomogram, and ground truth of the subtomogram pdb_id and rotation
we randomly select pairs of subtomograms of same structure, then rotate back according to ground truth so that the two subtomograms are aligned. Then calculate SNR through Pearson correlation according to Joachim Frank's estimation method

IMPORTANT: make sure following are similiar to experimenal data:  1) the subtomogram's size, 2) CTF and MTF, 3) tilt angle range, 4) the macromolecular complex of interest
'''
import json, os
import tomominer.io.file as IF
import tomominer.geometry.ang_loc as GAL
import tomominer.geometry.rotate as GR
import numpy as N
import scipy.stats as SS

if __name__ == "__main__":
    with open('instance_pair_snr__op.json') as f:   op = json.load(f)
    with open(op['model_file']) as f:    dj = json.load(f)
   
    # read into memory subtomograms of a particular pdb_id, and rotate according to ground truth
    for pdb_id in op["pdb_id"]:
	dj_pdb = []
	for tom_id in dj.keys():
	    tmp_dj = [_ for _ in dj[tom_id]["instances"] if ('pdb_id' in _) and (_['pdb_id'] == pdb_id)]
	    for i in tmp_dj:	i["tomogram_id"] = int(tom_id)
	    dj_pdb.extend(tmp_dj)
        assert len(dj_pdb) > 0
        vrs = []
	k = op["subtomogram_size"][0]/2
        tom_file = op["tomogram_files"]
	v_tom = N.load(tom_file.replace("I", tom_id))
        for i,d in enumerate(dj_pdb):
	    if (d["x"][0]-k)<0 or (d["x"][1]-k)<0 or (d["x"][2]-k)<0 or (d["x"][0]+k)>v_tom.shape[0] or (d["x"][1]+k)>v_tom.shape[1] or (d["x"][2]+k)>v_tom.shape[2]:
		continue
	    v = v_tom[d["x"][0]-k:d["x"][0]+k, d["x"][1]-k:d["x"][1]+k, d["x"][2]-k:d["x"][2]+k]
            vr = GR.rotate(v, rm = GAL.rotation_matrix_zyz(d['angle']).T, default_val=v.mean())       # use .T to transpose rotation matrix in order to rotate back
	    print(v.shape, vr.shape)
            vrs.append(vr)
       
            if i < 20:
                # this is only for debugging through visual inspection
                IF.put_mrc(v, '/tmp/v-%05d.mrc'%(i,), overwrite=True)
                IF.put_mrc(vr, '/tmp/vr-%05d.mrc'%(i,), overwrite=True)

    	print 'sampling', op['pair_num'], 'pairs from', len(vrs), 'aligned subtomograms'
    	with open(os.path.join("./snr_estimate", op['out_file']+pdb_id+".json"), 'w') as f:   
	    snr = []
	    print >>f, pdb_id
            for i in range(op['pair_num']):
                while True:
                    i0 = N.random.randint(len(vrs))
                    i1 = N.random.randint(len(vrs))
                    if i0 != i1:        break

                c, _ = SS.pearsonr(vrs[i0].flatten(), vrs[i1].flatten())        # pearson correlation between two aligned subtomograms
                s = c / (1 - c)         # SNR according to Equation 3.61 of book     frank Three-dimensional Electron Microscopy of Macromolecular Assemblies
		snr.append(s)
                print >>f, s
	    print >>f, "Average:", N.mean(snr)

'''
related code and analysis
~/ln/tomominer/tomominer/pursuit/multi/analysis/instance_pair_snr.py
'''
