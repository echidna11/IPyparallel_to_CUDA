#!/usr/bin/env python

'''
sometimes two copies of same tomogram have different origins, in such case, we can adjust the embedding obtained from one tomogram to fit another.

~/ln/tomominer/tomominer/template/search/analysis/instance_embedding_amira__origin_shift.py

'''

import os, json

import tomominer.io.file as IF
import tomominer.geometry.surface.amira.util as GSAA


def main():
    with open('instance_embedding_amira__origin_shift__op.json') as f:      op = json.load(f)

    if not os.path.isdir(op['surface new out dir']):    os.makedirs(op['surface new out dir'])
    

    mrc_org = IF.read_mrc_header(op['tomogram original'])['MRC']
    mrc_new = IF.read_mrc_header(op['tomogram new'])['MRC']

    for fn in os.listdir(op['surface original dir']):
        fnp = os.path.join(op['surface original dir'], fn)
        if not os.path.isfile(fnp):     continue
        if not fn.endswith('.surf'):        continue

        print 'loading', fnp
        with open(fnp) as f:       s = GSAA.surf_parse(f)

        x = s['vertices']
        x[:,0] -= mrc_org['xorg']       ;       x[:,0] += mrc_new['xorg']
        x[:,1] -= mrc_org['yorg']       ;       x[:,1] += mrc_new['yorg']
        x[:,2] -= mrc_org['zorg']       ;       x[:,2] += mrc_new['zorg']

        fnp_new = os.path.join(op['surface new out dir'], fn)
        assert not os.path.exists(fnp_new)
        print 'writing', fnp_new
        with open(fnp_new, 'w') as f:     GSAA.export_surf_ascii_simple(s, f)


if __name__ == '__main__':
    main()


