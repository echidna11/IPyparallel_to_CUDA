#!/usr/bin/env python

'''
given one field name (can be subfield names), and a json list file 
select the records that is equal to only a particular value

~/ln/tomominer/tomominer/json/list/select_range.py

'''



import json
from tomominer.json_.list_.util import get_field_val

import numpy as N


def main():

    with open('json_list_select_range__op.json') as f:   op = json.load(f)


    with open(op['file in']) as f:     rec_json = json.load(f)

    rng = sorted(op['range'])


    rec_new = []

    for ri, rec in enumerate(rec_json):
        val = get_field_val(rec, op['fields'])

        if val < rng[0]:        continue
        if val > rng[1]:        continue

        rec_new.append(rec)

    with open(op['file out'], 'w') as f:       json.dump(rec_new, f, indent=2)

    print len(rec_new), 'entries selected'


if __name__ == '__main__':
    main()

