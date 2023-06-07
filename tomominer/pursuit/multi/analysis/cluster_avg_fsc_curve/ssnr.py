#!/usr/bin/env python



# extract ssnr based fsc from tomominer.pursuit results (filtered cluster_average_select.pickle)


import os
import json
import cPickle as pickle

if __name__ == '__main__':


    with open('cluster_average_select__filtered.pickle'), 'rb') as f:   cs = pickle.load(f)

    ti = cs['tk_info']

    fsc = {}

    for c in ti:        fsc[c] = [float(_) for _ in ti[c]['fsc']]

    with open('ssnr__out.json', 'w') as f:     json.dump(fsc, f, indent=-1)

