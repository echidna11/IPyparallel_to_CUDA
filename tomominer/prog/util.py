'''
utility functions for python programming

~/ln/tomominer/tomominer/prog/util.py
'''



'''
arbitrary nested loop,
see http://code.activestate.com/recipes/502194-a-generator-for-an-arbitrary-number-of-for-loops/

usage example:
for i in multi_for(map(xrange, [2, 2, 1, 2])):      print i

'''
def multi_for(iterables):
    if not iterables:
        yield ()
    else:
        for item in iterables[0]:
            for rest_tuple in multi_for(iterables[1:]):
                yield (item,) + rest_tuple


'''
# test code

from tomominer.common.prog_util import multi_for
for i in multi_for([[1,2], ['a', 'b'], 'def']):     print i

'''


'''
similiar to multi_for()
'''
def mrange(minvec, maxvec=None):
    if maxvec is None:
        maxvec = minvec
        minvec = [0] * len(maxvec)
    vec = list(minvec)
    unitpos = len(vec) - 1
    maxunit = maxvec[unitpos]
    _tuple = tuple
    while 1:
        if vec[unitpos] == maxunit:
            i = unitpos
            while vec[i] == maxvec[i]:
                vec[i] = minvec[i]
                i -= 1
                if i == -1:
                    return
                vec[i] += 1
        yield _tuple(vec)
        vec[unitpos] += 1

