#!/usr/bin/env python


'''

select entries of particular range in a particular axis

~/ln/tomominer/tomominer/template/search/filter/select_axis_range.py

'''


import json

def main():
    with open('select_axis_range__op.json') as f:      op = json.load(f)

    with open(op['file in']) as f:      dj = json.load(f)

    dj_new = []
    for d in dj:
        if d['peak']['loc'][op['axis']] < op['range'][0]:        continue
        if d['peak']['loc'][op['axis']] > op['range'][1]:        continue
        dj_new.append(d)

    with open(op['file out'], 'w') as f:      json.dump(dj_new, f, indent=2)

    print len(dj_new), 'out of', len(dj), 'entries selected'

if __name__ == '__main__':
    main()

