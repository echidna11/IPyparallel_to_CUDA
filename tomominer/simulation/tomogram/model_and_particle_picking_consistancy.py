#!/usr/bin/env python


# for each complex, count the number of instances in the model and the number of instances picked by particle picking
# usage example: model_and_particle_picking_consistancy.py models.json data_config.json

if __name__ == '__main__':

    import sys
    model_file = sys.argv[1]
    data_file = sys.argv[2]

    import json
    with open(model_file) as f:     m = json.load(f)
    with open(data_file) as f:      d = json.load(f)


    mc = {}         # counter for model
    for mi in m:
        for r in m[mi]['instances']:
            pid = r['pdb_id']
            if pid not in mc:       mc[pid] = 0
            mc[pid] += 1

    md = {}         # counter for data of picked particles
    for r in d:
        pid = r['pdb_id']
        if pid not in md:       md[pid] = 0
        md[pid] += 1


    for pid in (set(mc.keys()) | set(md.keys())):
        if pid not in mc:   mc[pid] = 0
        if pid not in md:   md[pid] = 0
        print '%s\t%d\t%d'%(pid, mc[pid], md[pid])



