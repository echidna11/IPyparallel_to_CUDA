'''
utility funtions for dictionary operation

~/ln/tomominer/tomominer/common/dict_util.py
'''



'''
get a field inside a dictionary hierachy h, accoridng to a list of keywords
'''

def hi_get(h, ks):
    for k in ks:
        h = h[k]
    return h

'''
get a field to v inside a hierachy h, accoridng to a list of keywords
'''

def hi_set(h, ks, v):

    hp = None
    for k in ks:
        hp = h
        h = h[k]

    hp[k] = v


'''
delete one key from hierarchy dictionary
'''

def hi_del(h, ks):
    hp = None
    for k in ks:
        hp = h
        h = h[k]

    del hp[k]



