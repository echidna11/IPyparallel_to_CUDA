#!/usr/bin/env python

'''
collect path to all patterns generated at the final step of mpp,
output format should be readable by embedding and other analysis programs, specifically,
~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/data_prepare_from_fsc_collect.py

also collect all subtomogram alignments


~/ln/tomominer/tomominer/pursuit/multi/batch/pose_norm/analysis/embed/collect_data.py
'''


import json, copy, pickle

def main():
    with open('collect_data__op.json') as f:     op = json.load(f)
    with open(op['stat in']) as f:      st_in = json.load(f)
    with open(op['peak in']) as f:      p_in = json.load(f)
    with open(op['pattern in']) as f:      avg_in = json.load(f)
    avg_in = set([_['subtomogram'] for _ in avg_in])

    p_in = {_['subtomogram']:_ for _ in p_in}

    dj_tem = {}         # data for templates
    for st_in_t in st_in:
        with open(st_in_t['file_stat']) as f:      fs = json.load(f)
        ps = fs['passes']
        ps = {int(_):ps[_] for _ in ps}
        ps = ps[fs['pass_i_current']]

        with open(ps['cluster_average_select_file']) as f:      cas = pickle.load(f)

        st = cas['selected_templates']
        ti = cas['tk_info']

        for st_i in st:
            k = st[st_i]['subtomogram']
            dj_tem[k] = ti[k]['data_json']


    # collect all final alignments
    dj = []
    for st_in_t in st_in:
        with open(st_in_t['file_stat']) as f:      fs = json.load(f)
        ps = fs['passes']
        ps = {int(_):ps[_] for _ in ps}
        ps = ps[fs['pass_i_current']]


        with open(ps['data_json_file']) as f:   dj_t = json.load(f)

        for d in dj_t:
            d['peak'] = p_in[d['subtomogram']]['peak']
            if 'uuid' in p_in[d['subtomogram']]:            d['uuid'] = p_in[d['subtomogram']]['uuid']

            dj.append(d)

    print 'total subtomogram number', len(dj)


    # filter alignments to only those belong to a pattern
    dj_tem_set = []
    for k in dj_tem:        dj_tem_set.extend(_['subtomogram'] for _ in dj_tem[k])
    dj_tem_set = set(dj_tem_set)


    dj = [_ for _ in dj if _['subtomogram'] in dj_tem_set]
    with open(op['data out'], 'w') as f:        json.dump(dj, f, indent=2)
    print 'total number of subtomograms in a pattern', len(dj)



if __name__ == '__main__':
    main()



'''
example embedding
/auto/cmb-08/fa/mxu/proj/imaging/electron/frequent_structure/data/out/analysis/jensen/tocheva/130401-complex/classification/classification/pn/0000/analysis/plot/embed
'''


