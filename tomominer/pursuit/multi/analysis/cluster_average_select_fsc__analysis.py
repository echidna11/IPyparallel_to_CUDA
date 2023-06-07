#!/usr/bin/env python


'''

analysis of cluster_average_select_fsc(), to see why some clusters are selected

'''

if __name__ == '__main__':

    import json
    with open('../classify_config.json') as f:      op = json.load(f)

    import cPickle as pickle
    with open('cluster_info.pickle', 'rb') as f:    ci = pickle.load(f)

    from tomominer.pursuit.multi.util import tomominer.cluster_average_select_fsc

    cluster_average_select_fsc(None, cluster_info=ci, op=op, debug=True)

