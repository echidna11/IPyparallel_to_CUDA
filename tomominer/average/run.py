#!/usr/bin/env python



import os, json


from tomominer.parallel.queue_master import QueueMaster
import tomominer.common.obj as CO
from tomominer.io.cache import Cache


import tomominer.average.main_routine as AM

if __name__ == '__main__':

    qhost = 'localhost'
    qport = 5011


    # just get a generic object that can contain a cache
    self = CO.Object()

    tmp_dir = os.getenv('TMP_DIR');         assert      tmp_dir is not None
    self.cache = Cache(tmp_dir=tmp_dir)
    self.runner = QueueMaster(qhost, qport)
 
    with open('average__op.json') as f:     op = json.load(f)
    with open(op['data_json']) as f:    dj = json.load(f)


    AM.do_average(self=self, op=op, data_json=dj, out_dir=op['out_dir'])

