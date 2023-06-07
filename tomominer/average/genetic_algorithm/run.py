#!/usr/bin/env python



import os, json


from tomominer.parallel.queue_master import QueueMaster
import tomominer.common.obj as CO
from tomominer.io.cache import Cache
from tomominer.io.redis_proxy import Redis

import tomominer.average.genetic_algorithm.main as AGM

if __name__ == '__main__':

    with open('average__op.json') as f:     op = json.load(f)

    # just get a generic object that can contain a cache
    self = CO.Object()


    self.cache = Cache(tmp_dir=op['options']['tmp_dir'])
    self.runner = QueueMaster(op['options']['network']['qhost'], op['options']['network']['qport'])
    self.pool = None
 

    with open(op['data_file']) as f:    dj = json.load(f)


    AGM.do_average(self=self, op=op, data_json=dj)

