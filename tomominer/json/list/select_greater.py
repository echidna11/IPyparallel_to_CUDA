#!/usr/bin/env python

'''
given one field name (can be subfield names), and a json list file 
select the records that is equal to only a particular value

~/ln/tomominer/tomominer/json_/list_/select_greater.py

'''



import json
from tomominer.json_.list_.util import get_field_val

import numpy as N


def main():

    with open('json_list_select_greater__op.json') as f:   op = json.load(f)

    if 'take abs' in op:        print 'take abs', op['take abs']

    with open(op['file in']) as f:     rec_json = json.load(f)


    rec_new = []

    for ri, rec in enumerate(rec_json):
        val = get_field_val(rec, op['fields'])

        if ('take abs' in op) and op['take abs']:   val = N.abs(val)

        if (op['mode'] == 'greater') and (val < op['value']):     continue
        if (op['mode'] == 'less') and (val > op['value']):     continue

        rec_new.append(rec)

    with open(op['file out'], 'w') as f:       json.dump(rec_new, f, indent=2)

    print len(rec_new), 'entries selected'


if __name__ == '__main__':
    main()

