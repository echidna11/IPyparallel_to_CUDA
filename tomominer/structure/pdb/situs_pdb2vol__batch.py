#!/usr/bin/env python


import situs_pdb2vol as SP

# automatically scan pdb files and convert them to density maps, and save in to a pickle file.  This is done in parallel
# Alternatively, save in matlab format that is same as Bsoft.pdb2em_batch_convert_test() .

def batch_processing(op):

    # walk through every subdir, find all pdb files
    import os
    extension = '.pdb'
    pdb_path = {}
    for root, sub_folders, files in os.walk(op['pdb_dir']):
        for file_t in files:
            if not file_t.endswith(extension):     continue

            pdb_id = file_t[: len(file_t) - len(extension)]

            assert      (pdb_id + extension) == file_t
            assert      pdb_id not in pdb_path          # the pdb_id must be unique

            pdb_path[pdb_id] = os.path.join(root, file_t)


    from multiprocessing.pool import Pool
    pool = Pool()

    import copy
    pool_results = []
    for pdb_id in pdb_path:
        for spacing in op['spacing_s']:
            for resolution in op['resolution_S']:

                op_t = copy.deepcopy(op)
                op_t['pdb_id'] = pdb_id
                op_t['pdb_file'] = pdb_path[pdb_id]

                assert      'resolution' not in op_t;       op_t['resolution'] = resolution
                assert      'spacing' not in op_t;       op_t['spacing'] = spacing

                pool_results.append(      pool.apply_async(func=SP.convert, kwds={'op':op_t})         )


    re = {}
    for r in pool_results:
        cre = r.get()

        pdb_id = cre['pdb_id']
        resolution = cre['resolution']
        spacing = cre['spacing']

        if pdb_id not in re:        re[pdb_id] = {}

        if spacing not in re[pdb_id]:        re[pdb_id][spacing] = {}

        assert resolution not in re[pdb_id][spacing]

        re[pdb_id][spacing][resolution] = cre

    import pickle
    with open(op['out_file'], 'wb') as f:   pickle.dump(re, f, protocol=-1)



if __name__ == '__main__':
    import sys
    op_file = 'situs_pdb2vol__batch__op.json'
    import json
    with open(op_file) as f:    op = json.load(f)
    batch_processing(op)

