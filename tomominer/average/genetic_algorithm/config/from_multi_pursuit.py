#!/usr/bin/env python


'''

generate multiple configuration files from multi pursuit results generated using
~/ln/tomominer/tomominer/pursuit/multi/recursive/config_prepare.py


IMPORTANT:  currently we just hard code different options



~/ln/tomominer/tomominer/average/genetic_algorithm/config/from_multi_pursuit.py

'''


import os, json, copy

def main():

    stat = {}
    stat['op file'] = os.path.abspath('average_genetic_algorithm_config_from_multi_pursuit__op.json')

    with open(stat['op file']) as f:     op = json.load(f)
    op['out dir'] = os.path.abspath(op['out dir'])

    with open(op['average op file']) as f:   avg_op = json.load(f)

    with open(op['pursuit recursive config prepare stat file']) as f:     prs = json.load(f)

    stat['clusters'] = {}
    for c in prs['clusters']:           stat['clusters'][c] = prepare_one_cluster(out_dir=os.path.join(op['out dir'], '%03d'%(int(c), )), data_file=prs['clusters'][c]['tem__data_config_file'], avg_op=avg_op, repeat_num=op['repeat num'])

    with open(op['stat file out'], 'w') as f:       json.dump(stat, f, indent=2)

    print sum([len(stat['clusters'][_]) for _ in stat['clusters']]), 'configurations prepared'



def prepare_one_cluster(out_dir, data_file, avg_op, repeat_num):
    
    stat = []

    for aligned in [False, True]:
        for smooth in [False, True]:
            for repeat_i in range(repeat_num):

                op_dir_t = os.path.join(out_dir, ('align' if aligned else 'noalign'), ('smooth' if smooth else 'nosmooth'), '%04d'%(repeat_i,))
                if not os.path.isdir(op_dir_t):        os.makedirs(op_dir_t)

                avg_op_t = copy.deepcopy(avg_op)
                avg_op_t['data_file'] = data_file
                avg_op_t['out_dir'] = os.path.join(op_dir_t, 'out')
                if not os.path.isdir(avg_op_t['out_dir']):        os.makedirs(avg_op_t['out_dir'])

                avg_op_t['random_initial_alignment'] = not aligned
                
                assert 'smooth' in avg_op_t['averaging']
                if not smooth:
                    del avg_op_t['averaging']['smooth']
                    assert 'smooth' not in avg_op_t['averaging']

                avg_op_t_file = os.path.join(op_dir_t, 'average__op.json')
                if not os.path.isfile(avg_op_t_file):
                    print 'generating', avg_op_t_file
                    with open(avg_op_t_file, 'w') as f:     json.dump(avg_op_t, f, indent=2)
                else:
                    print 'ignoring', avg_op_t_file

                stat.append({'aligned':aligned, 'smooth':smooth, 'repeat_i':repeat_i, 'op_dir':op_dir_t, 'op_file':avg_op_t_file, 'out_dir':avg_op_t['out_dir']})


    return stat


if __name__ == '__main__':
    main()

