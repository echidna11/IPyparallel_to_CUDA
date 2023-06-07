#!/usr/bin/env python



'''

for batch calculating gold standard FSC

prepare configuration, including splitting data, and prepare a statistics file


~/ln/tomominer/tomominer/pursuit/multi/analysis/statistics/gold_standard_fsc/config_prepare.py

'''


import os, json, pickle, random, shutil


def main():
    with open('config_prepare__op.json') as f:       op = json.load(f)

    op['out_dir'] = os.path.abspath(op['out_dir'])
    if os.path.isdir(op['out_dir']):        shutil.rmtree(op['out_dir'])
    if not os.path.isdir(op['out_dir']):        os.makedirs(op['out_dir'])

    with open(op['collect_fsc_across_multiple_experiments__out'], 'rb') as f:       cfame = pickle.load(f)                   # usually       collect_fsc_across_multiple_experiments__out.pickle

    # prepare random split
    s = {}
    for fn in cfame:
        s[fn] = {}
        for c in cfame[fn]:

            dj = cfame[fn][c]['dj']
            random.shuffle(dj)
            
            st = {}
            for i, d in enumerate(dj):
                l = i % 2
                if l not in st:   st[l] = []
                st[l].append(d)

            print 'data size', len(dj), [len(st[_]) for _ in st]
            
            s[fn][c] = st

    conf_l = []         # list of configurations
    for fn in cfame:
        for c in cfame[fn]:
            cfame_t = cfame[fn][c]
            conf_t = {'file_stat':fn, 'cluster':c, 'count':len(conf_l), 'template':cfame_t['template']}
            conf_l.append(conf_t)


    for co in conf_l:
        co['root_dir'] = os.path.join(op['out_dir'], '%04d'%(co['count']))
        if not os.path.isdir(co['root_dir']):        os.makedirs(co['root_dir'])

        co['out_dir'] = os.path.join(co['root_dir'], 'out')
        if not os.path.isdir(co['out_dir']):        os.makedirs(co['out_dir'])



    for co in conf_l:
        co['split'] = os.path.join(co['root_dir'], 'split.pickle')
        with open(co['split'], 'wb') as f:       pickle.dump(s[co['file_stat']][co['cluster']], f, protocol=-1)


    with open(op['stat_out'], 'w') as f:        json.dump(conf_l, f, indent=2)


if __name__ == '__main__':
    main()



'''
related code

/home/rcf-47/mxu/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_fsc_curve/half_split.py
'''


