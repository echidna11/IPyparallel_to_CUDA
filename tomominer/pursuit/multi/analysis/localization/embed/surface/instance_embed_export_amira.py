#!/usr/bin/env python



"""
export embeded instances (generated using instance_embed.py) in Amira format


~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/instance_embed_export_amira.py

"""


import os, json, pickle, shutil

import tomominer.geometry.surface.amira.util as GSAA

def main():
    with open('instance_embed_export_amira__op.json') as f:     op = json.load(f)

    with open(op['mesh file'], 'rb') as f:      s = pickle.load(f)

    if os.path.isdir(op['out dir']):        shutil.rmtree(op['out dir'])

    os.makedirs(op['out dir'])

    for i in s:
        ie = s[i]['fv']
        with open(os.path.join(op['out dir'], '%03d.surf'%(i,)), 'w') as f:    GSAA.export_surf_ascii_simple(ie, f)


if __name__ == '__main__':
    main()

