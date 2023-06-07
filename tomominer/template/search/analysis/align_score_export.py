#!/usr/bin/env python



'''
export alignment score of a group of subtomograms into a matrix format, so that can be analysed using R or matlab


~/ln/tomominer/tomominer/template/search/analysis/align_score_export.py

'''



import os, json, pickle

def main():
    with open('align_score_export__op.json') as f:     op = json.load(f)

    fns = {}
    djs = {}
    for i, fn in enumerate(op['data']):
        fns[i] = fn

        with open(fn) as f:        dj = json.load(f)
        dj = {_['subtomogram']:_['score'] for _ in dj}
        djs[i] = dj

    with open(op['out file'], 'wb') as f:   pickle.dump({'djs':djs, 'fns':fns}, f, protocol=-1)


if __name__ == '__main__':
    main()




'''
# Example R commands for inspection, executed using rPython installed on u64 virtual machine


sshfs hpc-cmb:/staging/fa/mxu /home/mxu/tmp -o reconnect,idmap=user         # first we need to mount a shared folder

library('rPython')
python.exec('def load(fn):\n\timport pickle\n\twith open(fn, \'rb\') as f:\n\t\to = pickle.load(f)\n\treturn o')
python.exec('o = load(\'/home/mxu/tmp/align_score_export__out.pickle\')')
python.exec('fns = o[\'fns\']')
python.exec('djs = o[\'djs\']')


python.exec('t0 = [djs[0][_] for _ in djs[0]]')
python.exec('t1 = [djs[0][_] for _ in djs[0] if _ in djs[1]]')
python.exec('t2 = [djs[0][_] for _ in djs[0] if _ in djs[2]]')


t0 = python.get('t0')
t1 = python.get('t1')
t2 = python.get('t2')

boxplot(list(t0,t1,t2), names=c('Ribosome template', 'Pattern 1', 'Pattern 2'), ylab='Alignment scores')

'''



'''
# example rpy2 commands for inspection
'''



