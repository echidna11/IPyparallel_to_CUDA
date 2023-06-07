#!/usr/bin/env python

'''
given one field name (can be subfield names, seperated by :), and a top number, and a number of json record files,
for each file, select top number of records, and combine them together and write to stdout

usage example:
~/ln/tomominer/tomominer/json_/list/select_top.py highest peak:val 500 `find . -name "info.json"` > data_config_top_500.json


~/ln/tomominer/tomominer/json_/list_/select_top.py

'''


import json
import numpy as np
from tomominer.json.list.util import get_field_val


if __name__ == '__main__':

    with open('json_list_select_top__op.json') as f:   op = json.load(f)

    with open(op['file in']) as f:     rec_json = json.load(f)

    rec_merge = []
    if ('top num' in op) and (len(rec_json) > op['top num']):
        # extract the value of fields in record
        field_val = [None] * len(rec_json)
        for ri, rec in enumerate(rec_json):
            field_val[ri] = float(get_field_val(rec, op['fields']))

        if op['mode'] == 'highest':
            inds = np.argsort( -np.array(field_val) )
        elif op['mode'] == 'lowest':
            inds = np.argsort( np.array(field_val) )
        else:
            raise Exception('mode = ' + op['mode'])


        for ind_t in inds[:op['top num']]:      rec_merge.append(rec_json[ind_t])

    else:
        rec_merge.extend(rec_json)

    with open(op['file out'], 'w') as f:   json.dump(rec_merge, f, indent=2)

    print len(rec_merge), 'entries selected'


