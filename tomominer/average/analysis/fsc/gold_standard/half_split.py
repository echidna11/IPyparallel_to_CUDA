#!/usr/bin/env python


'''
given a set of subtomograms, randomly split into two halves
Then can use programs under following folder to calculate FSC
~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_fsc_curve/half_avg_independent.py


~/ln/tomominer/tomominer/average/analysis/fsc/gold_standard/half_split.py
'''



import json, copy, random



def main():
    with open('half_split__op.json') as f:      op = json.load(f)

    with open(op['data in']) as f:      dj = json.load(f)

    random.shuffle(dj)         # randomly shuffle dj, so that the result is more objective

    s = {}
    for i, d in enumerate(dj):
        l = i % 2
        if l not in s:   s[l] = []
        s[l].append(d)


    with open(op['data out'], 'w') as f:     json.dump({0:s}, f, indent=2)



if __name__ == '__main__':
    main()
    




'''
related code

~/ln/tomominer/tomominer/pursuit/multi/analysis/cluster_avg_fsc_curve/half_split.py
'''

