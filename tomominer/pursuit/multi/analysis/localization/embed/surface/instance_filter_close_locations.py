#!/usr/bin/env python




'''
if two predicted instances are closer than certain distance, remove the one with smaller score.
motivation: in the plotted embedding, occationally there are overlapping instances, dont know why. Maybe due to filtering in particle picking?


~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/instance_filter_close_locations.py
'''



import json
import numpy as N

from scipy.spatial.distance import cdist, pdist

def main():
    with open('instance_filter_close_locations__op.json') as f:     op = json.load(f)

    dist_cut = float(op['distance cut'])

    with open(op['data in']) as f:      dj = json.load(f)

    scores = N.array([_['score'] for _ in dj], dtype=N.float)
    x = N.array([_['peak']['loc'] for _ in dj], dtype=N.float)
    inds = N.array(xrange(len(dj)), dtype=N.int)

    keep = N.ones(len(dj)).astype(N.int)
    for i in xrange(len(dj)):
        if keep[i] != 1:     continue

        print '\r', 'progress', i, len(dj), '  ',

        x0 = x[i,:].reshape((1,3))
        d = cdist(x[keep == 1, :], x0)

        inds_t = inds[keep == 1][d.flatten() <= dist_cut]
        inds_t = inds_t[inds_t != i]

        if (len(inds_t) > 0) and (N.max(scores[inds_t]) >= scores[i]):      keep[i] = 0



    dj_new = [dj[_] for _ in xrange(len(dj)) if (keep[_] == 1)]
    with open(op['data out'], 'w') as f:        json.dump(dj_new, f, indent=2)

    print '\r', len(dj_new), 'out of', len(dj), 'entries are kept', '                 '

    x_new = N.array([_['peak']['loc'] for _ in dj_new], dtype=N.float)
    print 'minimum pariwise distance of kept instances', pdist(x_new).min()


if __name__ == '__main__':
    main()

