#!/usr/bin/env python




'''
batch submit jobs

ignore those jobs that has been marked as done?

~/ln/tomominer/tomominer/template/search/standard_scanning/batch/run_batch.py

'''


import os, sys, json, subprocess, uuid
import cPickle as pickle


def qsub(job_file, op):

    proot_file = os.path.abspath(op['proot'])
    assert os.path.isfile(proot_file)

    exec_file = os.path.abspath(op['exec_file'])
    assert os.path.isfile(exec_file)

    job_file = os.path.abspath(job_file)
    assert os.path.isfile(job_file)

    job_dir = os.path.dirname(job_file)
    idt = str(uuid.uuid4())
    args = ['qsub', 
            '-l', op['resource'], 
            '-o', os.path.join(job_dir, 'qsub-out-'+idt), 
            '-e', os.path.join(job_dir, 'qsub-err-'+idt),
            '-q', op['queue'],
            '-V',
            '-v', 'exec_file=%s,job_file=%s'%(exec_file, job_file),
            proot_file
            ]

    # print args
    subprocess.call(args=args, stdout=sys.stdout, stderr=sys.stderr)
     

def main():
    with open('run_batch__op.json') as f:       op = json.load(f)
    
    with open(op['config_file']) as f:      cop = json.load(f)
    with open(cop['config_stat_out_file']) as f:      tfl = json.load(f)['job_file_list']

    submit_count = 0
    for t in tfl:
        with open(t['job_file']) as f:     j = pickle.load(f)
        if os.path.isfile(j['stat_out']):   continue            # ignore those completed jobs
        qsub(t['job_file'], op=op['qsub'])
        submit_count += 1

        if ('submit_count__max' in op) and (submit_count >= op['submit_count__max']):    break


if __name__ == '__main__':
    main()

