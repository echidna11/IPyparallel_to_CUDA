#!/usr/bin/env python



'''
given any parameter file written in json, we generate different configurations by varying some of the parameters defined in option

~/ln/tomominer/tomominer/common/param/vary_param_generator.py

'''

import os, copy, json
from Queue import Queue
from tomominer.common.prog_util import multi_for
from tomominer.common.dict_util import hi_get, hi_set, hi_del

def generate():
    
    stat = {}

    with open('common_param_vary_param_generator__op.json') as f:       op = json.load(f)           # op['vary'] are those parameters that we can choose from multiple values,  op['disable'] are those options that we can chose to disable (delete)

    if 'verbose' not in op:       op['verbose'] = False


    with open(op['option file in']) as f:      op_org = json.load(f)

    op_vary = vary_config_indexing(op['vary'])
    stat['op_vary'] = op_vary


    stat['config'] = {}

    param_count = 0
    for pdi in multi_for(map(xrange, [2]*len(op['disable']))):
        for pvi in multi_for(map(xrange, op_vary['param_nums'])):

            op_t = get_one_config(op=op_org, vary_key=op_vary['key'], vary_list=op_vary['list'], pvi=pvi, disable_key=op['disable'], pdi=pdi)
            if op_t is None:        continue

            if op['verbose']:       print pdi, pvi, op_t

            op_dir = os.path.join(op['out dir'], '%04d'%(param_count,))
            if not os.path.isdir(op_dir):       os.makedirs(op_dir)
            op_file_t = os.path.join(op_dir, op['option file name out'])
            with open(op_file_t, 'w') as f:       json.dump(op_t, f, indent=2)

            stat['config'][param_count] = {'param_count':param_count, 'op_file':op_file_t, 'pdi':pdi, 'pvi':pvi}

            param_count += 1


    with open(op['stat file out'], 'w') as f:       json.dump(stat, f, indent=2)


'''
assign each option a index, and collect corresoponding indics, also calculate number of choices for each option
'''
def vary_config_indexing(op):
    op = copy.deepcopy(op)

    param_key = {}
    vary_list = {}

    i = 0
    q = Queue()

    assert type(op) is dict
    for k in op:            q.put([k])

    while not q.empty():
        ks = q.get()

        o = hi_get(op, ks)

        if type(o) is dict:
            for k in o:
                ks_t = copy.deepcopy(ks)
                ks_t.append(k)
                q.put(ks_t)
            continue

        if type(o) is list:
            param_key[i] = ks
            vary_list[i] = copy.deepcopy(o)
            i += 1
            continue

        raise Exception('Dont know how to process this type')

    param_nums = [None] * (max([_ for _ in vary_list]) + 1)
    for i in vary_list:        param_nums[i] = len(vary_list[i])

    return  {'key':param_key, 'list':vary_list, 'param_nums':param_nums}
    


def get_one_config(op, vary_key, vary_list, pvi, disable_key, pdi):

    op = copy.deepcopy(op)
    
    for i in range(len(pdi)):
        if pdi[i] == 0:
            try:
                hi_del(h=op, ks=disable_key[i])
            except KeyError:
                pass

    set_value_count = 0         # count the number of values set, if 0, this means no field can be set value, in that case we ignore this configuration
    for i in range(len(pvi)):
        try:
            hi_set(h=op, ks=vary_key[i], v=vary_list[i][pvi[i]])
            set_value_count += 1
        except KeyError:
            # if some key does not exist, this means it is deleted, in such case we do not wang to set any variations inside the corresponding subfields
            pass

    return (op if (set_value_count > 0) else None)






if __name__ == '__main__':
    generate()




'''

# test file and code

vim common_param_vary_param_generator__op.json

{
    "option file in": "./op.json",
    "option file name out": "op.json",
    "stat file out": "common_param_vary_param_generator__stat.json",

    "out dir": "./pv",

    "verbose": true,

    "vary": {
        
        "a": {
            "ab":[2,3]
        },

        "b": {
            "d": {
                "k": ["la", "rm"]
            }
        }

    },

    "disable": [
        ["a"],
        ["a", "ab"],
        ["b", "d", "g"]
    ]
}



vim op.json

{
    "a": {
        "ab": 1,
        "c": "d"
    },

    "b": {
        "a": 0.5,
        "d": {
            "k": "la",

            "g": 1.3
        }
    }
}



'''




