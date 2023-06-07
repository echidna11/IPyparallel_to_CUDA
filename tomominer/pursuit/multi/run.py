#!/usr/bin/env python

import os, sys, json, warnings
from tomominer.parallel.queue_master import QueueMaster

import tomominer.common.obj as CO
from tomominer.io.cache import Cache


def main():

    warnings.filterwarnings('error')

    with open('./pursuit-op.json') as f:       op = json.load(f)
    
    # just get a generic object that can contain a cache
    self = CO.Object()

    self.pool = None
  
    self.cache = Cache(tmp_dir=op['options']['tmp_dir'])

    self.runner = QueueMaster(op['options']['network']['qhost'], op['options']['network']['qport'])

    print 'loading ', op['data_file']
    with open(op['data_file']) as f:      data_json = json.load(f)

    # convert relative path to absolute path, if needed
    for d in data_json:
        if not os.path.isabs(d['subtomogram']):     d['subtomogram'] = os.path.abspath(os.path.join(os.path.dirname(op['data_file']), d['subtomogram']))
        if not os.path.isabs(d['mask']):     d['mask'] = os.path.abspath(os.path.join(os.path.dirname(op['data_file']), d['mask']))

    del op['data_file']

    from tomominer.pursuit.multi.main import pursuit
    pursuit(self=self, op=op, data_json=data_json)



if __name__ == '__main__':
    main()


