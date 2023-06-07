#!/usr/bin/env python

# split each cluster into halves


import json, copy, random
import cPickle as pickle




def split(ti):
    s = {}
    for c in ti:
        dj = copy.deepcopy(ti[c]['data_json'])
        random.shuffle(dj)         # randomly shuffle dj, so that the result is more objective

        s[c] = {}
        for i, d in enumerate(dj):
            l = i % 2
            if l not in s[c]:   s[c][l] = []
            s[c][l].append(d)
    return s


if __name__ == '__main__':

    with open('cluster_average_select__filtered.pickle', 'rb') as f:    ti = pickle.load(f)['tk_info']

    s = split(ti)

    with open('half_split__out.json', 'w') as f:       json.dump(s, f, indent=2)

