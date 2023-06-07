#!/usr/bin/env python



'''
compare two groups of scores, in paired and other modes

~/ln/tomominer/tomominer/template/search/analysis/align_score_compare.py
'''



import os, json, uuid

import numpy as N
import scipy.stats as SS


def main():
    with open('align_score_compare__op.json') as f:     op = json.load(f)

    with open(op['data0']) as f:        dj0 = json.load(f)
    dj0 = {_['subtomogram']:_['score'] for _ in dj0}

    with open(op['data1']) as f:        dj1 = json.load(f)
    dj1 = {_['subtomogram']:_['score'] for _ in dj1}

    s0 = []
    s1 = []
    for k in set(dj0.keys()) & set(dj1.keys()):
        s0.append(dj0[k])
        s1.append(dj1[k])

    s0 = N.array(s0)
    s1 = N.array(s1)

    print 'common entity number', len(s0)

    print 'data0', 'mean', s0.mean(), 'median', N.median(s0), 'std', s0.std()
    print 'data1', 'mean', s1.mean(), 'median', N.median(s1), 'std', s1.std()

    print 'ttest_rel', SS.ttest_rel(s0, s1)
    print 'ranksum', SS.ranksums(s0, s1)
    print 'Mann-Whitney rank test', SS.mannwhitneyu(s0, s1)
    print 'wilcoxon', SS.wilcoxon(s0, s1)

    print 'pearson correlation', SS.pearsonr(s0, s1)


    # save data matrix as txt file to plot scatter plots
    print 'for further inspection, saving values to', op['values file out']
    N.savetxt(op['values file out'], N.vstack((s0,s1)).T, delimiter='\t')

    # save data as GSEA format for GSEA analysis


if __name__ == '__main__':
    main()




'''
# R commands for inspection

x = read.table('/tmp/value_out.txt')

plot(x[,1], x[,2], main="Alignment scores", xlab="Template", ylab="Pattern", pch=19)

'''

