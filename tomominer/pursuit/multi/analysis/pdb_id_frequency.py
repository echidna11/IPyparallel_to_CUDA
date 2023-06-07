#!/usr/bin/env python



# given a data_json file, count the number of items with particular pdb_id


import sys
import json

if __name__ == '__main__':

    pdb_label_file = sys.argv[1]
    data_json_file = sys.argv[2]

    with open(pdb_label_file) as f:     pl = json.load(f)

    pd = {}
    for r in pl:
        pd[r['cluster_label']] = r['pdb_id']


    with open(data_json_file) as f:     dj = json.load(f)

    cd = {}
    for r in dj:
        if not 'cluster_label' in r:       continue
        if r['cluster_label'] not in cd:       cd[r['cluster_label']] = []
        cd[r['cluster_label']].append(r)


    for l in cd:
        print pd[l], '\t', len(cd[l])


