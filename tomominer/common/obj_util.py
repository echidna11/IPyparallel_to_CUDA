'''
recursively apply an action to elementary objects inside a complex object.
currently supported object types: general object, dict, list

parameters: o: object,  f: action,      fop: options for the action

the main purpose is to recursively replace path names in complex objects, so that we can duplicate project into another folder
'''

import numpy as N
from collections import defaultdict as collections_defaultdict

def recursive_elementary_apply(o, f, fop=None, verbose=False):

    if o is None:   return None

    if  hasattr(o, '__call__'):     return None            # o is a function

    if fop is not None:
        o1 = f(o=o, op=fop)
    else:
        o1 = f(o=o)

    if o1 is not None:
        if verbose:     print o, o1
        return o1       # in such case, o is some eligible elementary objects

    # for object of type dict
    if (type(o) is dict) or (type(o) is collections_defaultdict):
        o1 = {}
        for k in o:
            
            if fop is not None:            k1 = f(o=k, op=fop)        # we also apply the action to the key of a dict       !!!
            else:       k1 = f(o=k)

            r = o[k]
            r1 = recursive_elementary_apply(o=r, f=f, fop=fop, verbose=verbose)
            if r1 is None:      r1 =r

            if k1 is not None:
                assert  k1 not in o1
            else:
                assert k not in o1
                k1 = k

            o1[k1] = r1

            del k, k1, r, r1
        return o1

    # for object of type list
    if type(o) is list:
        o1 = [None] * len(o)
        for i, t in enumerate(o):
            t1 = recursive_elementary_apply(o=t, f=f, fop=fop, verbose=verbose)
            if t1 is None:      t1 = t
            o1[i] = t1
            del i, t, t1

        return o1

    if type(o) is set:
        o1 = set()
        for t in o:
            t1 = recursive_elementary_apply(o=t, f=f, fop=fop, verbose=verbose)
            if t1 is None:  t1 = t
            assert  t1 not in o1
            o1.add(t1)
        return o1


    if not isinstance(o, (type, str, unicode, int, float, bool, N.ndarray)):
        # for generic object (not str, int, ... builtin_function_or_method, instancemethod, function etc), in such case we DO NOT make a copy!!!
        for k in dir(o):
            if k.startswith('_'):   continue        # we do not process internal attributes
            a = getattr(o, k)
            a1 = recursive_elementary_apply(a, f, fop, verbose=verbose)
            if a1 is not None:      setattr(o, k, a1)           # if a1 is not some method or non-replacable stuff, we set the corresponding attribute
            del k, a

    return o


def print_elementary_obj(o):

    print o
    return o



'''
# test code

import copy
import tomominer.common.obj_util as CO

a = CO.Object()
a.a = {'a':'apple', 'b':'task', 'c':'apask', '1':CO.Object(), 2:'land'}
a.a['1'].b = {'a':'apple', 'b':'task', 'c':'apask', 'apple_key':'apnana', 'task':'apple', 'task_key':['apply', 'dont ask me apple']}
a.apple = copy.deepcopy(a)

getattr(a, 'a')

import tomominer.common.string as CS
CO.recursive_elementary_apply(o=a, f=CS.print_string)


a1 = copy.deepcopy(a)
a2 = CO.recursive_elementary_apply(o=a1, f=CS.substring_replace, fop={'startswith':False, 'match':'ap', 'replace':'replace_bf_'}, verbose=True)

CO.recursive_elementary_apply(o=a, f=CS.print_string)
CO.recursive_elementary_apply(o=a1, f=CS.print_string)
CO.recursive_elementary_apply(o=a2, f=CS.print_string)

a3 = copy.deepcopy(a)
a4 = CO.recursive_elementary_apply(o=a3, f=CS.substring_replace, fop={'startswith':False, 'match':'ask', 'replace':'rep_swer_'})
CO.recursive_elementary_apply(o=a4, f=CS.print_string)

dir(a4)
dir(a4.apple)
a4.a
a4.a['1']
a4.a['1'].b

a4.apple.a
a4.apple.a['1']
a4.apple.a['1'].b


'''


