#!/usr/bin/env python




'''

automatically assign rainbow colors to selected averages


code similiar to
~/ln/tomominer/tomominer/simulation/tomogram/analysis/class_color_assign.py



~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/average_color_auto_assign.py

'''

import os, json

import tomominer.plot.color as PC


if __name__ == '__main__':
    
    with open('average_color_auto_assign__op.json') as f:    op = json.load(f)

    with open(op['average file']) as f:     avg = json.load(f)

   
    if False:
        n_avg = len([_ for _ in avg if 'subtomogram' in _])
        cm = PC.get_colormap(n_avg + 2)     # we remove the colors at two ends
    else:
        id_max = max([_['id'] for _ in avg])
        cm = PC.get_colormap(id_max + 1 + 2)     # we remove the colors at two ends

    colors = {}
    for a in avg:
        if not 'subtomogram' in a:        continue

        colors[a['subtomogram']] = {'color':list(        cm(a['id']+1)[:3]       ), 'id':a['id']}

    with open(op['color file out'], 'w') as f:       json.dump(colors, f, indent=2)

    print 'assigned colors to', len(colors), 'templates'


