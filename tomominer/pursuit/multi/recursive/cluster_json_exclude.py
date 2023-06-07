#!/usr/bin/env python


# exclude all subtomograms from all clusters of pursuit result, so that we can use remaining subtomograms to pursuit again
# this script can be use in combination with cluster_json_export.py, for consistancy, we can create a folder with name 'clus-rest'


'''

usage

~/ln/tomominer/tomominer/pursuit/multi/recursive/cluster_json_exclude.py    >   data_config.json

'''


if __name__ == '__main__':

    import json
    with open('cluster_json_exclude__config.json') as f:     op = json.load(f)

    with open(op['data_json__file']) as f:        pp = json.load(f)


    import cPickle as pickle
    with open(op['cluster_average_select__file'], 'rb') as f:   cs = pickle.load(f)

    st = cs['selected_templates']
    ti = cs['tk_info']

    tif = {_:ti[st[_]['subtomogram']]   for _ in st}            # keep only selected cluster's information, indixed by cluster id

    del cs, st, ti

    # collect the set of subtomograms
    ss = set()
    for c in tif:           ss.update(_['subtomogram'] for _ in tif[c]['data_json'])

    import sys
    json.dump([_ for _ in pp if (_['subtomogram'] not in ss)], sys.stdout, indent=2)


