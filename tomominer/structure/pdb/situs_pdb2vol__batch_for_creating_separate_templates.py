#!/usr/bin/env python

import os, sys, json, copy
from multiprocessing.pool import Pool
import tomominer.structure.pdb.situs_pdb2vol as SP
import tomominer.image.vol.util as IVU
import tomominer.io.file as IF

# automatically scan pdb files and convert them to density maps, and save as mrc file. This is done in parallel
    
def batch_processing(op):
    # walk through every subdir, find all pdb files
    extension = '.pdb'
    pdb_path = {}
    for root, sub_folders, files in os.walk(op['pdb_dir']):
        for file_t in files:
            if not file_t.endswith(extension):     continue
            pdb_id = file_t[: len(file_t) - len(extension)]
            assert      (pdb_id + extension) == file_t
            assert      pdb_id not in pdb_path          # the pdb_id must be unique
            pdb_path[pdb_id] = os.path.join(root, file_t)

    pool = Pool()
    pool_results = []
    for pdb_id in pdb_path:
        op_t = copy.deepcopy(op)
        op_t['pdb_id'] = pdb_id
        op_t['pdb_file'] = pdb_path[pdb_id]
        pool_results.append( pool.apply_async(func=SP.convert, kwds={'op':op_t}) )

    re = {}
    for r in pool_results:
        cre = r.get()
        pdb_id = cre['pdb_id']
        if pdb_id not in re:
            re[pdb_id] = cre

    for pdb_id in pdb_path:
        v = re[pdb_id]['map']
        if 'out_map_size' in op:
            v = IVU.resize_center(v=v, s=op['out_map_size'], cval=0.0)
        if 'intensity normalize' in op:
            if op['intensity normalize']['mode'] == 'max': v_n = v / v.max()
            elif op['intensity normalize']['mode'] == 'mean':  v_n = v / v.mean()
            elif op['intensity normalize']['mode'] == 'std':   v_n = v / v.std()
            else:   raise Exception('normalize mode')
            out_file_norm = pdb_id + "-" + str(op['spacing']) + "-" + str(op['resolution']) + "-map-normalized.mrc"
            IF.put_mrc(v_n, out_file_norm)

        out_file = pdb_id + "-" + str(op['spacing']) + "-" + str(op['resolution']) + "-map.mrc"
        IF.put_mrc(v, out_file)


if __name__ == '__main__':
    with open("situs_pdb2vol__batch__op.json") as f:    op = json.load(f)
    batch_processing(op)
