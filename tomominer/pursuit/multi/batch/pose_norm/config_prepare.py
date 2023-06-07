#!/usr/bin/env python



'''
given pose normalization clustering result, first split the data according to manual grouping, 
then create folder and save data

~/ln/tomominer/tomominer/pursuit/multi/batch/pose_norm/config_prepare.py
'''


import os, json, uuid

def main():
    with open('config_prepare__op.json') as f:          op = json.load(f)

    op['out dir'] = os.path.abspath(op['out dir'])

    with open(op['clustering_kmeans__op file']) as f:     ckp_op = json.load(f)
    
    with open(os.path.join(os.path.dirname(op['clustering_kmeans__op file']), ckp_op['input data json file'])) as f:     dj = json.load(f)
    with open(os.path.join(os.path.dirname(op['clustering_kmeans__op file']), ckp_op['output data json file'])) as f:     lbl = json.load(f)
    lbl = {_['subtomogram']:_['cluster_label'] for _ in lbl}


    with open(op['cluster groups']) as f:     cg = json.load(f)


    st = {}
    for grp_i, grp in enumerate(cg):
        grp_id = grp['id']
        assert grp_id == grp_i

        lbl_sel = set(grp['clusters'])

        djf = []
        for d in dj:
            d['cluster_label'] = lbl[d['subtomogram']]
            if d['cluster_label'] not in lbl_sel:       continue
            djf.append(d)

        st[grp_id] = {}
        st[grp_id]['uuid'] = str(uuid.uuid4())
        st[grp_id]['group_dir'] = os.path.join(op['out dir'], 'group', '%04d'%(grp_id,))
        if not os.path.isdir(st[grp_id]['group_dir']):       os.makedirs(st[grp_id]['group_dir'])

        st[grp_id]['data'] = os.path.join(st[grp_id]['group_dir'], 'data.json')
        with open(st[grp_id]['data'], 'w') as f:     json.dump(djf, f, indent=2)

        st[grp_id]['repeat'] = {}
        for rep_i in range(op['repeat num']):
            st[grp_id]['repeat'][rep_i] = {}
            st[grp_id]['repeat'][rep_i]['uuid'] = str(uuid.uuid4())
            st[grp_id]['repeat'][rep_i]['root_dir'] = os.path.join(st[grp_id]['group_dir'], 'rep', '%04d'%(rep_i,))
            st[grp_id]['repeat'][rep_i]['out_dir'] = os.path.join(st[grp_id]['repeat'][rep_i]['root_dir'], 'out')


    with open(op['stat out'], 'w') as f:        json.dump(st, f, indent=2)


if __name__ == '__main__':
    main()



'''
similiar code

/home/rcf-47/mxu/ln/tomominer/tomominer/average/genetic_algorithm/config/from_multi_pursuit.py


related code
~/ln/frequent_structure/data/out/method/pose-norm-level-set/scripts/cluster_selection.py

'''

