#!/usr/bin/env python




'''
embed instances in form of isosurface, generated using template_isosurface.py
assuming that there is only one template, and all subtomograms are aligned to that template. This is for template matching or averaging


~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/instance_embed__single_template.py

'''


import json, pickle

import tomominer.geometry.surface.util as GSU


def main():
    with open('instance_embed__single_template__op.json') as f:     op = json.load(f)

    with open(op['mesh file'], 'rb') as f:     s = pickle.load(f)
    assert  len(s) == 1

    with open(op['subtomogram file']) as f:     dj = json.load(f)       # Note: we assume the dj is already filtered with respect to a particular tomogram_id

    for d in dj:        d['mid_co'] = d['peak']['loc']
    
    se = {}
    
    print 'enbedding', len(dj), 'instances of', s[0]['subtomogram']
    se[0] = {}
    se[0]['subtomogram'] = s[0]['subtomogram']
    se[0]['fv'] = GSU.instance_embedding(dj=dj, ts=s[0]['fv'], tc=s[0]['center'], reverse=True)

    print 'saving trangular mesh to', op['mesh file out'],
    with open(op['mesh file out'], 'wb') as f:          pickle.dump(se, f, protocol=-1)
    print 'done'

if __name__ == '__main__':
    main()

