#!/usr/bin/env python




'''
use pose.normalize.segmentation.batch.run_json
'''

import os
import json
import copy
from multiprocessing.pool import Pool

import tomominer.pose.normalize.segmentation.batch.run_json as PNSBR

if __name__ == '__main__':

    with open('pose_normalize_segmentation_batch_run_json__op.json') as f:    op = json.load(f)
    with open(op['segmentation_op_file']) as f:           op['segmentation'] = json.load(f)
    del op['segmentation_op_file']

    with open('subtomogram_extraction__out.json') as f:     se = json.load(f)

    data_dir_out = os.path.abspath(op['data_dir_out'])

    pool = Pool()

    for tom_i in se:
        for sigma1 in se[tom_i]:
            for size_ratio in se[tom_i][sigma1]:
                op_t = copy.deepcopy(op)
                del op_t['data_dir_out']


                op_t['out_dir'] = os.path.join(data_dir_out, sigma1, size_ratio, tom_i)

                op_t['data_config_in'] = se[tom_i][sigma1][size_ratio]['info_file']
                op_t['data_config_out'] = os.path.join(op_t['out_dir'], 'info.json')

                PNSBR.main(op=op_t, pool=pool)



