
import numpy as np

# given a set of feautres in x (a matrix with columns as features and rows as data points), find the set of pareto optimal data points by minimizing each column
def pareto_min__set_sel(x):

    inds = set(range(len(x)))

    while True:

        non_optimal_removed = False
        for i in inds:
            for j in inds:
                if i==j:    continue
                
                if all(x[i] <= x[j]) and any(x[i] < x[j]): 
                    inds.remove(j) 
                    non_optimal_removed = True
                    break

            if non_optimal_removed:     break        
            
        if not non_optimal_removed:     break
    
    return inds        
                           

# given a set of vectors (supposed to be pareto optimal), perform linear regression, then calculate combined coefficient score
def pareto_min__regress(x):

    n_row = x.shape[0]

    xt = np.hstack( (np.ones( (n_row, 1) ), x) )

    n_col = xt.shape[1]

    ls = {}
    ls['x'], ls['residuals'], ls['rank'], ls['s'] = np.linalg.lstsq( xt[:,:(n_col-1)], -xt[:,(n_col-1)]  )       # perform least squares fitting, c_0 + c_1* x_0 + c_2* x_1 +...c_{k+1)* x_k = 0, fix the coefficient c_k for the last feature to be 1

    # construct coefficients for all
    co = np.append(ls['x'], 1.0)
    assert(all(co[1:] >= 0))        # make sure that all coefficients (except c_0) are non-negative, so that smaller score correspond with better solution

    # regression scores
    score = xt.dot(co)


    return {'coef':co, 'score':score, 'ls':ls}

# both selection and regression
def pareto_min__set_sel_regress(x):
    i = pareto_min__set_sel(x)
    i = list(i)

    re = pareto_min__regress(x[i,:])

    re['i'] = i
    re['x'] = x


    return re

