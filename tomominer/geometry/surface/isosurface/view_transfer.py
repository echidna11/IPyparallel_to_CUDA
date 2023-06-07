#!/usr/bin/env python



'''

transfer views of one record to another batch of ground truth


~/ln/tomominer/tomominer/geometry/surface/isosurface/view_transfer.py
'''


import json

def main():
    with open('view_transfer__op.json') as f:        op = json.load(f)

    with open(op['view in']) as f:      view = json.load(f)
    view = {_[op['key field']]:_['view'] for _ in view}

    with open(op['data in']) as f:      dj = json.load(f)

    for d in dj:        d['view'] = view[d[op['key field']]]

    with open(op['view out'], 'w') as f:        json.dump(dj, f, indent=2)




if __name__ == '__main__':
    main()


