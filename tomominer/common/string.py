
'''
utility for string manipulation
'''

import copy
def substring_replace(o, op=None):
    if (type(o) is not str) and (type(o) is not unicode):   return     None
    o = copy.copy(o)

    if op['startswith']:
        if not o.startswith(op['match']):   return
        o = op['replace'] + o[len(op['match']):]
    else:
        if op['match'] not in o:    return
        o = o.replace(op['match'], op['replace'])

    return o



def print_string(o):
    if (type(o) is not str) and (type(o) is not unicode):   return     None
    print o
    return o

def print_matched_string(o, op):
    if (type(o) is not str) and (type(o) is not unicode):   return     None
    if op['match'] not in o:    return  None
    print o
    return o
