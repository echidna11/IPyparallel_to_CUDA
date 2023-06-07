
import os
import json

import numpy as np
import random

import tomominer.io.file as iv
import general_util.vol as gvol
import tomominer.segmentation.thresholding as segt
import classify.config as config
from work_queue.queue_master import QueueMaster


# use one template to guide segmentation of a given set of subtomograms, parallel batch processing

def segment(op, data, out_dir, runner):

    cache_dir = os.getenv('CACHE_DIR')
    t = iv.get_mrc_cache_fs(str(op['template']['vol']), cache_dir)
    if 'mask' in op['template']:
        tm = iv.get_mrc_cache_fs(str(op['template']['mask']), cache_dir)
    else:
        tm = np.ones(t.shape, order='F')

    v = iv.get_mrc_cache_fs(data[0][0], cache_dir)
    if t.shape != v.shape:
        # sometimes the size of a subtomogram is different to the size of the template, in such case, resize the template
        if int(op['template']['resize']['default_val_mode']) == 0:
            t = gvol.resize_center(v=t, s=np.array(v.shape), cval=0.0)
        elif int(op['template']['resize']['default_val_mode']) == 1:
            t = gvol.resize_center(v=t, s=np.array(v.shape), cval=t.mean())
        else:
            raise AttributeError(op['template']['resize']['default_val_mode'])

        t = np.array(t, order='F')
        tm = gvol.resize_center(v=tm, s=np.array(v.shape), cval=0.0);      tm = np.array(tm, order='F')


    mmb_cutoffs = np.linspace(max(0, t.min()), t.max(), int(op['template']['segmentation']['cutoff_num']))
    mmb_cutoffs = mmb_cutoffs[1:];      mmb_cutoffs = mmb_cutoffs[:(len(mmb_cutoffs)-1)]
    mmb = segt.ms_mo_best( segt.ms_cutoffs(v=t, cutoffs=mmb_cutoffs) )
    t_seg = np.array( (t > mmb['best']['cutoff']), dtype=np.float, order='F' ) 



    if 'test' in op:
        if ('sample_num' in op['test']) and (op['test']['sample_num'] > 0) and (len(data) > op['test']['sample_num']):
            print 'testing the procedure using a subsample of %d subtomograms'%(op['test']['sample_num'])
            data = random.sample(data, op['test']['sample_num'])






    op['out_dir'] = out_dir

    tasks = []
    for v,m,a,l in data:
        tasks.append(runner.task('segment', vol_key=(v,m), t=t, tm=tm, t_seg=t_seg, op=op))

    # record results           
    data_json = []
    for res in runner.run__except(tasks):
        res = res.result

        json_rec = {}
        json_rec['subtomogram'] = res['vol_key']
        json_rec['mask'] = res['mask_key']

        # numpy arrays are not pickle-able.  Convert to a list of numbers.
        json_rec['angle'] = [_ for _ in (np.random.random(3) * (np.pi * 2))]
        json_rec['loc'] = [_ for _ in np.zeros(3)]
           
        data_json.append(json_rec)

    # write new data to disk
    data_json_file =  os.path.join(out_dir, 'data_config.json')
    with open(data_json_file, 'w') as f:
        json.dump(data_json, f, indent=2)




if __name__ == '__main__':


    import sys
    qhost = sys.argv[1]
    qport = 5011

    param_file = sys.argv[2]
    data_file  = sys.argv[3]
    out_dir = sys.argv[4]

    with open(param_file) as f:
        op = json.load(f)

    with open(data_file) as f:
        data_json = json.load(f)
    
    data = config.parse_data(data_json)
    
    runner = QueueMaster(qhost, qport)

    segment(op=op, data=data, out_dir=out_dir, runner=runner)


