#!/usr/bin/env python

'''
given one field name (can be subfield names), and a json list file 
remove those redundant records

~/ln/tomominer/tomominer/json/list/remove_redundant.py

'''



import json
from tomominer.json.list.util import get_field_val



def main():

    with open('json_list_remove_redundant__op.json') as f:   op = json.load(f)

    with open(op['file in']) as f:     rec_json = json.load(f)

    val_set = set() 

    rec_new = []

    for ri, rec in enumerate(rec_json):
        
        val = get_field_val(rec, op['fields'])
        if val in val_set:      continue
        val_set.add(val)

        rec_new.append(rec)

    with open(op['file out'], 'w') as f:       json.dump(rec_new, f, indent=2)

    print "New length of data:" ,len(rec_new)
    print "Unique Values for ", op['fields'], ": ", len(val_set)
    print "Original length of data: ", len(rec_json)


if __name__ == '__main__':
    main()


