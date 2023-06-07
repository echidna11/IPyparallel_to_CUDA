#!/usr/bin/env python


'''
template search through performing fast alignment between the template and the  extracted subtomograms

the main issue is about how to reduce the false positives
'''

import numpy as N

import tomominer.io.file as IOF
import tomominer.pursuit.multi.util as PMU

def align_single_batch(self, op):

    t = N.load(op['template_key']['subtomogram'])
    tm = N.load(op['template_key']['mask'])

    dj_re = []
    for d in op['dj']:
        v = IOF.read_mrc(d['subtomogram'])['value'].astype(N.float)
        vm = IOF.read_mrc(d['mask'])['value'].astype(N.float)

        a = PMU.align_vols_with_wedge(v1=t, m1=tm, v2=v, m2=vm, op=op['align'])

        if a['err'] is not None:
            # if the alignment fails, do not record result
            print   a['err']
            continue

        d['angle'] = [float(_) for _ in a['angle']]
        d['loc'] = [float(_) for _ in a['loc']]
        d['score'] = float(a['score'])

        dj_re.append(d)


    return dj_re




def main():

    import os, sys, json, uuid, shutil

    with open('template_search_fast_alignment__op.json') as f:      op = json.load(f)

    proj_id = str(uuid.uuid4())

    if not os.path.isdir(op['tmp_dir']):        os.makedirs(op['tmp_dir'])

    with open(op['data_json_file']) as f:     dj = json.load(f)

    # load and resize template, reverse intensity if needed. Store into a temporary folder. Also make and store a mask
    siz = IOF.read_mrc(dj[0]['subtomogram'])['value'].shape
    t = IOF.read_mrc(op['template_file'])['value'].astype(N.float);      t = t - t.mean()
    if op['template_reverse_intensity']:    t = -t

    import tomominer.geometry.rotate as GR
    if t.shape != siz:  t = GR.rotate(t, siz2=siz, default_val=t[0,0,0])        # we assume t has a constant intensity boarder, so we choose t[0,0,0] as default value when enlarging it



    template_file = os.path.join(op['tmp_dir'], 'template--%s--vol.npy'%(proj_id,))
    N.save(template_file, t)

    import tomominer.model.util as MU
    template_mask_file = os.path.join(op['tmp_dir'], 'template--%s--mask.npy'%(proj_id,))
    N.save(template_mask_file, MU.sphere_mask(siz).astype(N.float))
    


    # distribute task chunks
    from tomominer.parallel.queue_master import QueueMaster
    qhost = os.getenv('HOSTNAME')
    qport = 5011
    runner  = QueueMaster(qhost, qport)
    
    n_chunk = op['n_chunk']
    tasks = []
    dj_t = dj
    while dj_t:
        op_t = {'template_key':{'subtomogram':template_file, 'mask':template_mask_file}, 'dj':dj_t[:n_chunk], 'align':op['align']}
        tasks.append(runner.task( module='tomominer.template.search.fast_align.run', method='align_single_batch', kwargs={'op':op_t} ))

        dj_t = dj_t[n_chunk:]

    
    dj_re = []
    for r in runner.run__except(tasks):
        r = r.result
        dj_re.extend(r)

    dj_re = sorted(dj_re, reverse=True, key=lambda _:_['score'])

    with open(op['out_file'], 'w') as f:    json.dump(dj_re, f, indent=2)

    
    os.remove(template_file)
    os.remove(template_mask_file)



if __name__ == '__main__':
    main()

