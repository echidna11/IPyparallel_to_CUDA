#!/usr/bin/env python

'''

given pose normalization list L resulted from run_json.py, and a data json file J contains a list of subtomograms before pose normalization, filter L using J

'''


if __name__ == '__main__':

    import json
    with open('pose_normalize_segmentation_batch_filter_json__op.json') as f:   op = json.load(f)

    with open(op['data config file']) as f:   dj = json.load(f)
    with open(op['pose normalize file in']) as f:  pn = json.load(f)

    s = set([_['subtomogram'] for _ in dj])
    pn_t = [_ for _ in pn if _['subtomogram'] in s]

    with open(op['pose normalize file out'], 'w') as f:     json.dump(pn_t, f, indent=2)
