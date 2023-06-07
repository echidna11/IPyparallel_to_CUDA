#!/usr/bin/env python


'''
manually select a number of references, and find nearest subtomograms (peaks) from a list of subtomograms

~/ln/tomominer/tomominer/pick/reference/manual_select.py

'''


import os, json, copy, csv
import numpy as N

from scipy.spatial.distance import cdist


def main():
    with open('pick_reference_manual_select__op.json') as f:    op = json.load(f)

    with open(op['data in']) as f:      dj = json.load(f)

    if 'peaks' in op:
        peaks = op['peaks']
    else:
        peaks = load_peaks(path=op['peaks file'], default_tomogram_id=op['default tomogram id'], delimiter=op['peaks file delimiter'])

    print 'processing', len(peaks), 'manually selected peaks'


    res = []
    for r in peaks:

        c = N.reshape(N.array(r['loc']), (1,3))
        djt = [_ for _ in dj if _['tomogram_id'] == r['tomogram_id']]

        if len(djt) == 0:   continue

        x = [_['peak']['loc'] for _ in djt]
        x = N.array(x)

        dists = cdist(x, c).flatten()
        assert  len(dists) == len(djt)

        dists_i = N.argmin(dists)

        if dists[dists_i] > op['dist cutoff']:
            print 'ignored', r
            continue

        d = copy.deepcopy(djt[dists_i])
        d['peak']['reference'] = copy.deepcopy(r)
        d['peak']['reference']['distance_to_peak'] = dists[dists_i]

        res.append(d)

    print len(res), 'of', len(peaks), 'references selected'

    with open(op['data out'], 'w') as f:       json.dump(res, f, indent=2)


def load_peaks(path, default_tomogram_id, delimiter):
    with open(path) as f:
        r = csv.reader(f, delimiter=str(delimiter))
        ps = list(r)

    peaks = []
    for p in ps:
        if len(p) == 3:
            peaks.append({'loc':p, 'tomogram_id':default_tomogram_id})
        elif len(p) == 4:
            peaks.append({'loc':p[:3], 'tomogram_id':int(p[3])})

    return peaks


if __name__ == '__main__':
    main()


