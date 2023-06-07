#!/usr/bin/env python

# select templates of certain resolution and list

import sys, json, pickle

def select_from_situs_pickle(op):
    with open(op["map_file"]) as f:
        m = pickle.load(f)
    
    #pdb_ids = set(op['pdb_ids'])
    pdb_ids = set([_ for _ in m.keys() if _ not in op["reject pdbs"]])

    maps_sel = {}
    for i, pdb_id in enumerate(pdb_ids):
        resolution = op['resolution']
        spacing = op['spacing']
        assert      pdb_id not in maps_sel
        assert      pdb_id == m[pdb_id][spacing][resolution]['pdb_id']
        assert      float(resolution) == float(m[pdb_id][spacing][resolution]['resolution'])
        assert      float(spacing) == float(m[pdb_id][spacing][resolution]['spacing'])
        maps_sel[pdb_id] = { 'id':i, 'map':m[pdb_id][spacing][resolution]['map'], 'pdb_id':m[pdb_id][spacing][resolution]['pdb_id'], 'resolution':float(m[pdb_id][spacing][resolution]['resolution']), 'spacing':float(m[pdb_id][spacing][resolution]['spacing'])  }
    
    return maps_sel


def select_from_bsoft_mat(op):

    pdb_ids = set(op['pdb_ids'])
    import scipy.io as SIO
    m = SIO.loadmat(op['map_file']);    m = m['maps']
    maps_sel = {}

    for i0 in range(m.shape[0]):
        for i1 in range(m.shape[1]):
            for i2 in range(m.shape[2]):
                mt = m[i0, i1, i2][0,0]

                info = {'map':mt[0], 'header':mt[1], 'pdb_id':str(mt[2][0]), 'resolution':float(mt[3][0,0]), 'spacing':float(mt[4][0,0])}

                if info['pdb_id'] not in pdb_ids:   continue
                if info['resolution'] != float(op['resolution']):       continue
                if info['spacing'] != float(op['spacing']):         continue

                assert      info['pdb_id'] not in maps_sel
                maps_sel[info['pdb_id']] = info

    return maps_sel


def normalize(ms, op):
    print 'normalize()', op
    for pi in ms:
        print pi
        v = ms[pi]['map']

        if op['mode'] == 'max':
            v = v / v.max()
        elif op['mode'] == 'mean':
            v = v / v.mean()
        elif op['mode'] == 'std':
            v = v / v.std()
        else:
            raise Exception('normalize mode')

        ms[pi]['map'] = v


if __name__ == '__main__':

    with open('template_select__config.json') as f:        op = json.load(f)
    if True:
        maps_sel = select_from_situs_pickle(op)
    else:
        maps_sel = select_from_bsoft_mat(op)

    print 'selected', len(maps_sel), 'maps'
    if 'intensity normalize' in op:        normalize(maps_sel, op['intensity normalize'])
    with open('template_select.pickle', 'wb') as f:   pickle.dump(maps_sel, f)
