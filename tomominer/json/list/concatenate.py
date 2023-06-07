#!/usr/bin/env python

'''
concatenate all JSON lists specified in the parameters

~/ln/tomominer/tomominer/json/list/concatenate.py

'''

if __name__ == '__main__':

    import sys

    import json
    data = []
    for i, fn in enumerate(sys.argv):
        if i == 0: continue

        with open(fn) as f:     d = json.load(f)
        data.extend(d)

    json.dump(data, sys.stdout, indent=2)


