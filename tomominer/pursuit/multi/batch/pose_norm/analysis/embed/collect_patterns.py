#!/usr/bin/env python

'''
collect path to all patterns generated at the final step of mpp,
output format should be readable by embedding and other analysis programs, specifically,


~/ln/tomominer/tomominer/pursuit/multi/batch/pose_norm/analysis/embed/collect_patterns.py
'''


import json
import pickle

def main():
    with open('collect_patterns__op.json') as f:     op = json.load(f)
    with open(op['stat in']) as f:      st_in = json.load(f)

    avgs = []
    dj = []
    for st_in_t in st_in:
        with open(st_in_t['file_stat']) as f:      fs = json.load(f)
        ps = fs['passes']
        ps = {int(_):ps[_] for _ in ps}
        ps = ps[fs['pass_i_current']]
        with open(ps['cluster_average_select_file']) as f:      cas = pickle.load(f)

        cascf = cas['selected_templates_common_frame']
        for c in cascf:
            avgs.append({'subtomogram':cascf[c]['subtomogram'], 'mask':cascf[c]['mask']})

        with open(ps['data_json_file']) as f:   dj_t = json.load(f)
        dj.extend(dj_t)

    for i, a in enumerate(avgs):        a['id'] = i

    with open(op['stat out'], 'w') as f:        json.dump(avgs, f, indent=2)



if __name__ == '__main__':
    main()



'''
example embedding
/auto/cmb-08/fa/mxu/proj/imaging/electron/frequent_structure/data/out/analysis/jensen/tocheva/130401-complex/classification/classification/pn/0000/analysis/plot/embed
'''


