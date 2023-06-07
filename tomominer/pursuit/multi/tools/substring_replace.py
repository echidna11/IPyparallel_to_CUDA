#!/usr/bin/env python



'''
recursively replace substring with complex objects of all json and pickle files
'''


import os, json
import cPickle as pickle

import tomominer.common.obj_util as COU
import tomominer.common.string as CS



def replace(op):

    for root, dirs, files in os.walk(op['dir']):
        for f_t in files:
            f_full = os.path.join(root, f_t)

            if f_t.endswith('.pickle'):
                print f_full
                with open(f_full, 'rb') as f:        o = pickle.load(f)
                o = COU.recursive_elementary_apply(o=o, f=CS.substring_replace, fop={'startswith':op['startswith'], 'match':op['match'], 'replace':op['replace']}, verbose=op['verbose'])
                with open(f_full, 'wb') as f:        pickle.dump(o, f, protocol=-1)

            if f_t.endswith('.json'):
                print f_full
                with open(f_full) as f:     o = json.load(f)
                o = COU.recursive_elementary_apply(o=o, f=CS.substring_replace, fop={'startswith':op['startswith'], 'match':op['match'], 'replace':op['replace']}, verbose=op['verbose'])
                with open(f_full, 'w') as f:        json.dump(o, f, indent=2)




if __name__ == '__main__':
    with open('pursuit_multi_tools_substring_replace__op.json') as f:       op = json.load(f)
    replace(op)

