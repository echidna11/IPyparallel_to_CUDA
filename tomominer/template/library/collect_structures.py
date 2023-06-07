#!/usr/bin/env python


'''
simply go through a directory, and collect all pdb structures

~/ln/tomominer/tomominer/template/library/collect_structures.py
'''


import os, json

def main():
    with open('collect_structures__op.json') as f:      op = json.load(f)


    fs = {}
    
    for root_dir in op['root dirs']:
        for root, dirs, files in os.walk(root_dir):
            for n in files:
                pid, ext = os.path.splitext(n)
                if ext != '.pdb' :      continue
                fs[pid] = os.path.abspath(os.path.join(root, n))

    with open(op['stat out'], 'w') as f:     json.dump(fs, f, indent=2)


if __name__ == '__main__':
    main()


