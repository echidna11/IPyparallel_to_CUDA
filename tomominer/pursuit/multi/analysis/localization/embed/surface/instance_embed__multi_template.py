#!/usr/bin/env python

'''
embed instances in form of isosurfaces, generated using template_isosurface.py
~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/instance_embed__multi_template.py

'''

import json, pickle
import tomominer.geometry.surface.util as GSU
import tomominer.geometry.surface.amira.util as GSAA

def main():
    with open('instance_embed__multi_template__op.json') as f:     op = json.load(f)

    print 'loading', op['mesh file']
    with open(op['mesh file'], 'rb') as f:     s = pickle.load(f)

    print 'loading', op['subtomogram file']
    with open(op['subtomogram file']) as f:     dj = json.load(f)       # Note: we assume the dj is already filtered with respect to a particular tomogram_id

    dj = sorted(dj, key=lambda _:(-_['score']))         # order according score in descending order

    for d in dj:        d['mid_co'] = d['peak']['loc']

    if 'tomogram file' in op:
        # use mrc information so that the exported vrml can be directly used by amira in combination with the tomogram
        import tomominer.io.file as IF
        mrc = IF.read_mrc_header(op['tomogram file'])['MRC']
    else:
        mrc = None


    ss = [None] * len(dj)
    inds = {}
    for i in s:
        inds[i] = [_ for _ in xrange(len(dj)) if (dj[_]['template']['subtomogram'] == s[i]['subtomogram'])]
        if len(inds[i]) == 0:      continue

        print '\nembedding', len(inds[i]), 'instances of', s[i]['subtomogram']
        ss_t = GSU.instance_embedding_list(dj=[dj[_] for _ in inds[i]], ts=s[i]['fv'], tc=s[i]['center'], reverse=True)

        for inds_i, inds_t in enumerate(inds[i]):      ss[inds_t] = ss_t[inds_i]



    if 'surface vertice min dist' in op:
        print 'remove intersecting instances'
        non_intersect_flag = GSU.instance_find_nonintersect_instances(dj, ss, min_dist=op['surface vertice min dist'])

        for i in xrange(len(ss)):
            if non_intersect_flag[i] != 1:     ss[i] = None

   
    total_count = 0
    se = {}
    for i in s:
        if len(inds[i]) == 0:      continue

        se[i] = {}
        se[i]['subtomogram'] = s[i]['subtomogram']

        ss_t = [ss[_] for _ in inds[i] if (ss[_] is not None)]
        se[i]['fv'] = GSU.concat_list(ss_t)

        if mrc is not None:         se[i]['fv']['vertices'] = GSAA.vertice_location_transform(se[i]['fv']['vertices'], mrc)

        total_count += len(ss_t)

    print 'in total', total_count, 'instances embedded'

    print 'saving trangular mesh to', op['mesh file out']
    with open(op['mesh file out'], 'wb') as f:          pickle.dump(se, f, protocol=-1)
    print 'done'


if __name__ == '__main__':
    main()

