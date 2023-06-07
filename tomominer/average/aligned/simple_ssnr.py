#!/usr/bin/env python


'''
simple average of a set of already aligned subtomograms. 

Modified from  
~/ln/tomominer/tomominer/average/genetic_algorithm/main.py



~/ln/tomominer/tomominer/average/aligned/simple_ssnr.py

'''


import os, json, copy

import numpy as N

from tomominer.parallel.queue_master import QueueMaster
import tomominer.common.obj as CO
from tomominer.io.cache import Cache


import tomominer.pursuit.multi.util as PMU


def get_ssnr(self, op, data_json):

    with open(op['data_file']) as f:     data_json = json.load(f)


    n_chunk = op['options']['n_chunk']

    #-----------------------------------------------------------
    # calculate SSNR based FSC

    cluster_ssnr_fsc__op = {}
    cluster_ssnr_fsc__op['ssnr'] = copy.deepcopy(op['ssnr'])
    cluster_ssnr_fsc = PMU.cluster_ssnr_fsc(self=self, clusters={0:data_json}, n_chunk=n_chunk, op=cluster_ssnr_fsc__op)
    print 'fsc_sum', N.sum(cluster_ssnr_fsc['fsc'][0])

    return {'fsc':cluster_ssnr_fsc['fsc'][0].tolist(), 'fsc_sum':N.sum(cluster_ssnr_fsc['fsc'][0]), 'ssnr':cluster_ssnr_fsc['ssnr'][0].tolist()}



def main():

    with open('ssnr__op.json') as f:     op = json.load(f)

    # just get a generic object that can contain a cache
    self = CO.Object()

    self.cache = Cache(tmp_dir=op['options']['tmp_dir'])
    self.runner = QueueMaster(op['options']['network']['qhost'], op['options']['network']['qport'])
    self.pool = None

    with open(op['data_file']) as f:    dj = json.load(f)

    st = get_ssnr(self=self, op=op, data_json=dj)

    with open(op['stat_out'], 'w') as f:        json.dump(st, f, indent=2)



if __name__ == '__main__':
    main()




'''
related code

~/ln/tomominer/tomominer/average/aligned/simple_average.py
'''


