#!/usr/bin/env python


'''
align a list of subtomograms against multiple references.

IMPORTANT: in order to keep information concise, we do not keep the keys in data list and template list, but just indicis!!





~/ln/tomominer/tomominer/template/search/fast_align/reference/multi.py

'''


import numpy as N

import tomominer.io.file as IOF
import tomominer.pursuit.multi.util as PMU

def align_single_batch(self, op):

    inds = op['inds']
    rks = op['reference_keys']

    dj_re = {}
    for ti in rks:
        t = N.load(rks[ti]['subtomogram'])
        tm = N.load(rks[ti]['mask'])

        for di, d in enumerate(op['dj']):

            v = IOF.read_mrc_vol(d['subtomogram']).astype(N.float)
            vm = IOF.read_mrc_vol(d['mask']).astype(N.float)

            a = PMU.align_vols_with_wedge(v1=t, m1=tm, v2=v, m2=vm, op=op['align'])

            if a['err'] is not None:
                # if the alignment fails, do not record result
                print   a['err']
                continue

            d['angle'] = [float(_) for _ in a['angle']]
            d['loc'] = [float(_) for _ in a['loc']]
            d['score'] = float(a['score'])

            if inds[di] not in dj_re:                dj_re[inds[di]] = {}

            dj_re[inds[di]][ti] = {'angle':d['angle'], 'loc':d['loc'], 'score':d['score']}


    return dj_re



def main():

    import os, sys, json, uuid, shutil, random

    with open('reference_search_fast_alignment_reference_multi__op.json') as f:      op = json.load(f)

    proj_id = str(uuid.uuid4())

    if not os.path.isdir(op['tmp_dir']):        os.makedirs(op['tmp_dir'])

    with open(op['data_json_file']) as f:     dj = json.load(f)
    with open(op['reference_json_file']) as f:     djr = json.load(f)

    if 'test' in op:
        if ('sample_num' in op['test']) and (op['test']['sample_num'] > 0) and (len(dj) > op['test']['sample_num']):
            print 'testing the procedure using a subsample of %d subtomograms'%(op['test']['sample_num'])
            dj = random.sample(dj, op['test']['sample_num'])


    siz = IOF.read_mrc(dj[0]['subtomogram'])['value'].shape
    djr_new = {}
    for i, d in enumerate(djr):
        # load and resize reference, reverse intensity if needed. Store into a temporary folder. Also make and store a mask
        t = IOF.read_mrc_vol(d['subtomogram']).astype(N.float);      t = t - t.mean()
        if op['reference_reverse_intensity']:    t = -t

        assert t.shape == siz

        reference_file = os.path.join(op['tmp_dir'], 'reference--%s--%03d--vol.npy'%(proj_id,i))
        N.save(reference_file, t)

        m = IOF.read_mrc_vol(d['mask']).astype(N.float)
        mask_file = os.path.join(op['tmp_dir'], 'reference--%s--%03d--mask.npy'%(proj_id,i))
        N.save(mask_file, m)

        djr_new[i] = {'subtomogram':reference_file, 'mask':mask_file}


    # distribute task chunks
    from tomominer.parallel.queue_master import QueueMaster
    qhost = os.getenv('HOSTNAME')
    qport = 5011
    runner  = QueueMaster(qhost, qport)
    
    n_chunk = op['n_chunk']
    tasks = []
    inds = list(range(len(dj)))
    while inds:
        inds_t = inds[:n_chunk]

        op_t = {'reference_keys':djr_new, 'inds':inds_t, 'dj':[dj[_] for _ in inds_t], 'align':op['align']}
        tasks.append(runner.task( module='tomominer.template.search.fast_align.reference.multi', method='align_single_batch', kwargs={'op':op_t} ))

        inds = inds[n_chunk:]

    
    dj_re = [None] * len(dj)
    for re in runner.run__except(tasks):
        re = re.result
        for i in re:
            assert  dj_re[i] is None
            dj_re[i] = re[i]

    for d in dj_re:     assert d is not None

    with open(op['out_file'], 'w') as f:    json.dump(dj_re, f, indent=2)

    for i in djr_new:
        os.remove(djr_new[i]['subtomogram'])
        os.remove(djr_new[i]['mask'])



if __name__ == '__main__':
    main()


'''

code modified from

~/ln/tomominer/tomominer/template/search/fast_align/run.py

'''


