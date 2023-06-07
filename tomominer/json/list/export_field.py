#!/usr/bin/env python



'''
export a particular field, print one value for each line

~/ln/tomominer/tomominer/json/list/export_field.py
'''


import json
from tomominer.json_.list_.util import get_field_val


def main():
    with open('json_list_export_field__op.json') as f:   op = json.load(f)

    with open(op['list in']) as f:     dj = json.load(f)

    vs = []
    for d in enumerate(dj):
        v = get_field_val(d, op['fields'])
        vs.append(v)

    with open(op['list out'], 'w') as f:
        for v in vs:        print >>f, v


if __name__ == '__main__':
    main()


