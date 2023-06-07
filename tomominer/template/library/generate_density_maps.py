#!/usr/bin/env python



'''
given a collection of pdb structures, generate density maps


~/ln/tomominer/tomominer/template/library/generate_density_maps.py

'''


import os, json, copy

import tomominer.structure.pdb.situs_pdb2vol as SP
import tomominer.image.vol.util as IVU
import tomominer.io.file as IF

def main():
    with open('generate_density_maps__op.json') as f:      op = json.load(f)
    op['out dir'] = os.path.abspath(op['out dir'])

    if not os.path.isdir(op['out dir']):    os.makedirs(op['out dir'])

    with open(op['pdb in']) as f:       ps = json.load(f)

    fs = {}    
    for pid in ps:
        cop = copy.deepcopy(op['situs'])
        cop['pdb_file'] = ps[pid]
        cop['pdb_id'] = pid
        r = SP.convert(cop)

        v = r['map']

        if 'out_map_size' in cop:
            v = IVU.resize_center(v=v, s=cop['out_map_size'], cval=0.0)

        fs[pid] = os.path.join(op['out dir'], pid + '.mrc')
        IF.put_mrc(v, fs[pid])

    with open(op['stat out'], 'w') as f:     json.dump(fs, f, indent=2)

if __name__ == '__main__':
    main()

