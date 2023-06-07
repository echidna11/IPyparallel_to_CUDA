#!/usr/bin/env python


'''
merge two lists according to key inside particular field

~/ln/tomominer/tomominer/json_/list_/merge_key.py

'''


import json
import tomominer.common.dict_util as CDU

if __name__ == '__main__':
    with open('merge_key__op.json') as f:     op = json.load(f)

    with open(op['extracted list']) as f:    dj_extracted = json.load(f)
    dj_extracted = {CDU.hi_get(_, op['key field']):_ for _ in dj_extracted}
    print "Extracted list: ", len(dj_extracted)

    with open(op["pass list"]) as f:    dj_pass = json.load(f)
    dj_pass = {CDU.hi_get(_, op['key field']):_ for _ in dj_pass}
    print "Pass list: ", len(dj_pass)

    with open(op['original list']) as f:    dj_original = json.load(f)
    dj_original = {CDU.hi_get(_, op['key field']):_ for _ in dj_original}
    print "Original list: ", len(dj_original)

    dj_final = []
    for k in dj_pass:
        if k not in dj_extracted:   continue
        if k not in dj_original:    continue

        r_p = dj_pass[k]
        r_o = dj_original[k]

        for rk in r_o:
            if rk in r_p:
                continue
            r_p[rk] = r_o[rk]

        dj_final.append(r_p)

    print "Final list: ", len(dj_final)
    with open(op['output list'], 'w') as f:     json.dump(dj_final, f, indent=2)
