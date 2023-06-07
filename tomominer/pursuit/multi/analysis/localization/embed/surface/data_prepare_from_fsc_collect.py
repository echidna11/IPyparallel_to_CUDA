#!/usr/bin/env python


'''


from output of
~/ln/tomominer/tomominer/pursuit/multi/analysis/statistics/collect_fsc_across_multiple_experiments.py

and a list of manually selected averages, export the averages and collected subtomogram and alignment info into data_json



~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/data_prepare_from_fsc_collect.py

'''




import json, csv
import cPickle as pickle


def main():
    with open('data_prepare_from_fsc_collect__op.json') as f:       op = json.load(f)


    with open(op['peak file']) as f:        pk = json.load(f)

    pk = {_['subtomogram']:_ for _ in pk}


    with open(op['collect file in'], 'rb') as f:        co = pickle.load(f)

    cod = {}
    for fn in co:
        for c in co[fn]:
            co_t = co[fn][c]
            cod[co_t['template']['subtomogram']] = co_t


    with open(op['selected template file']) as f:
        r = csv.reader(f)
        st = list(r)

    st = [_[op['selected template file column index']] for _ in st]

    dj = []
    for t in st:
        for d in cod[t]['dj']:
            d['template'] = cod[t]['template']      # we force replace the template information
            d['peak'] = pk[d['subtomogram']]['peak']
            dj.append(d)

    with open(op['data file out'], 'w') as f:       json.dump(dj, f, indent=2)
    with open(op['average file out'], 'w') as f:    json.dump([{'subtomogram':_, 'id':id_} for id_, _ in enumerate(st)], f, indent=2)

    print 'selected', len(st), 'averages', 'containing', len(dj), 'entries'


if __name__ == '__main__':
    main()


