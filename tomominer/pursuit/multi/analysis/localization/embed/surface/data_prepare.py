#!/usr/bin/env python


"""
Prepare a combined subtomogram alignment file (merged with peak location information) for embedding and plotting


~/ln/tomominer/tomominer/pursuit/multi/analysis/localization/embed/surface/data_prepare.py

"""


import json

def main():
    with open('data_prepare__op.json') as f:        op = json.load(f)
    
    with open(op['peak file']) as f:        pk = json.load(f)

    pk = {_['subtomogram']:_ for _ in pk}

    dj_all = []
    for dft in op['data files in']:
        with open(dft['data']) as f:    djt = json.load(f)

        dj = []
        for d in djt:
            if 'template' in dft:                d['template'] = dft['template']                # we force replacing the template
            d['peak'] = pk[d['subtomogram']]['peak']
            dj.append(d)

        assert  len(dj) > 0

        dj_all.extend(dj)

    with open(op['data file out'], 'w') as f:        json.dump(dj_all, f, indent=2)
    print len(dj_all), 'entries collected'

    
    avg = set([_['template']['subtomogram'] for _ in dj_all])
    avg = list(avg)
    avg = [{'subtomogram':_} for _ in avg]

    with open(op['average file out'], 'w') as f:        json.dump(avg, f, indent=2)
    print len(avg), 'averages collected'


if __name__ == "__main__":
    main()



