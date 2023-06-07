#!/usr/bin/env python



'''
filter the subtomograms using a reference list

~/ln/tomominer/tomominer/template/search/filter/reference_set.py
'''

import json

def main():
    with open('reference_set__op.json') as f:     op = json.load(f)

    with open(op['main']) as f:        dj = json.load(f)

    with open(op['reference']) as f:        djr = json.load(f)

    s = set(_['subtomogram'] for _ in djr)


    with open(op['out'], 'w') as f:     json.dump([_ for _ in dj if _['subtomogram'] in s], f, indent=2)


if __name__ == '__main__':
    main()


