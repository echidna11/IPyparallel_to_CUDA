#!/usr/bin/env python




'''

read the output of
~/ln/tomominer/tomominer/average/genetic_algorithm/config/from_multi_pursuit.py

then create multiple processes to process each configurations


~/ln/tomominer/tomominer/average/genetic_algorithm/config/from_multi_pursuit__batch_run.py

'''


import os, sys, json, subprocess, time


def main():
    with open('from_multi_pursuit__batch_run__op.json') as f:       op = json.load(f)

    with open(op['config stat file']) as f:       conf_stat = json.load(f)

    ps = []
    for c in conf_stat['clusters']:
        for s in conf_stat['clusters'][c]:
            p = {}
            p['stat'] = s
            p['args'] = [op['average program']]
            p['stdout'] = open(os.path.join(s['out_dir'], 'log.out'), 'a')
            p['stderr'] = open(os.path.join(s['out_dir'], 'log.err'), 'a')
            p['proc'] = subprocess.Popen(args=p['args'], stdout=p['stdout'], stderr=p['stderr'], cwd=s['op_dir'])
            
            ps.append(p)


    while True:
        for pst in ps:                  pst['poll'] = pst['proc'].poll()

        active_proc = [_ for _ in ps if (_['poll'] is None)]

        print '\r%5d active processes'%(len(active_proc),),
        sys.stdout.flush()

        if len(active_proc) == 0:    break

        time.sleep(5)



    for pst in ps:
        pst['stdout'].close()
        pst['stderr'].close()



if __name__ == '__main__':
    main()

