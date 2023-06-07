import numpy as N

# simple linear regression to get slope
def one_var(x, y):
    A = N.array([x, N.ones(len(x))])
    w = N.linalg.lstsq(A.T, y)
    return w[0]

