#!/usr/bin/env python


# in order for other scripts' fast loading, read cluster_average_select.pickle and keep only the info of the selected clusters of that pass.

import os
import json
import cPickle as pickle


if __name__ == '__main__':

    with open('cluster_average_select_filter__config.json') as f:       op = json.load(f)

    with open(op['cluster_average_select__file'], 'rb') as f:   cs = pickle.load(f)

    st = cs['selected_templates']
    ti = cs['tk_info']

    tif = {_:ti[st[_]['subtomogram']]   for _ in st}

    with open('cluster_average_select__filtered.pickle', 'wb') as f:        pickle.dump({'selected_templates':st, 'tk_info':tif}, f, protocol=-1)

