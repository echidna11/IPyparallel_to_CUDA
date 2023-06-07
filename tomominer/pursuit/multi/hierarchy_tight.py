'''
given several sets, and hierarchical clustering, find tight clusters according to these sets.
see 
https://docs.google.com/document/d/1QzllJyZ26MDDDRPmpb8hXQu95bULwV9oZy9LQdOxyas/edit

'''


import numpy as N
import tomominer.pursuit.multi.hierarchy_ssnr_calculate_select as PMH


'''
parameters: hc: return of hierarchy_ssnr_calculate_select.hierarchical_clustering()
'''

def tight_clusters(hc, sets):
    print 'tight_clusters()'
    
    hi = hc['info']
    root = hc['root_id']

    cover_set_ids = {}
    tight_set_ids = {}

    for c in sets:
        s = set(sets[c])

        # recursively look down the most tight cluster that contains all nodes in s
        cid = root          # current cluster id
        while True:
            if (hi[cid]['left'] is not None) and (s <= set(hi[hi[cid]['left']]['nodes'])):
                cid = hi[cid]['left']
                continue

            if (hi[cid]['right'] is not None) and (s <= set(hi[hi[cid]['right']]['nodes'])):
                cid = hi[cid]['right']
                continue

            break

        cover_set_ids[c] = cid
        size_d = {_:len(hi[_]['nodes']) for _ in PMH.hierarchy_traverse(hi=hi, root=cid)}            # collect size of all subclusters starting from root cid

        tight_set_ids[c] = []
        size_d_min_similiar = N.min(N.abs(N.array([size_d[_] for _ in size_d]) - len(s)))
        for i in size_d:
            if N.abs(size_d[i] - len(s)) != size_d_min_similiar:        continue
            tight_set_ids[c].append(i)

    all_ids = [cover_set_ids[_] for _ in cover_set_ids]         # collection of all hierarchical ids, used for correcponding data_json clusters
    for c in tight_set_ids:     all_ids.extend(tight_set_ids[c])
    all_ids = set(all_ids)

    return {'cover_set_ids':cover_set_ids, 'tight_set_ids':tight_set_ids, 'all_ids':all_ids}


'''
given sizes of seed clusters (i.e. those specific clusters selected in last iteration) in seed_sizes, we can first find cover_set_ids using tight_clusters() with hc, we can then traverse each root in cover_set_ids and remove all subclusters that is smaller than the corresponding seed_size
'''
def collect_large_subclusters_by_seed_cluster_sizes(hc, cover_set_ids, seed_sizes, ratio=1.0):
    ratio = float(ratio)

    hi = hc['info']

    ids = {}
    for c in cover_set_ids:
        ids[c] = []
        for i in PMH.hierarchy_traverse(hi=hi, root=cover_set_ids[c]):
            if len(hi[i]['nodes']) < (seed_sizes[c] * ratio):  continue
            ids[c].append(i)

    all_ids = []
    for c in ids:       all_ids.extend(ids[c])

    return {'ids':ids, 'all_ids':all_ids}




'''
simply finding the corresponding node ids
parameters:     dj_all: data json that corresponds to node ids,         djs: particular data json sets
'''
def data_json_to_node_id(dj_all, djs):
    dj_all_dict = {}
    for i, d in enumerate(dj_all):      dj_all_dict[d['subtomogram']] = i

    sets = {}
    for c in djs:        sets[c] = [dj_all_dict[_['subtomogram']] for _ in djs[c]]
    
    return sets

