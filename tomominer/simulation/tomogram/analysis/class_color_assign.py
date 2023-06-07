#!/usr/bin/env python



'''
given a list of complexes, automatically assign different colors to different complexes



~/ln/tomominer/tomominer/simulation/tomogram/analysis/class_color_assign.py

'''


import os, json, pickle

import tomominer.plot.color as PC


if __name__ == '__main__':
    
    with open('class_color_assign__op.json') as f:    op = json.load(f)

    with open(op['map file'], 'rb') as f:     maps = pickle.load(f)         # usually template_select.pickle

    
    cm = PC.get_colormap(len(maps) + 2)     # we remove the colors at two ends
    colors = {}
    for i, pid in enumerate(maps.keys()):        colors[pid] = cm(i+1)[:3]

    with open(op['color file out'], 'w') as f:       json.dump(colors, f, indent=2)

