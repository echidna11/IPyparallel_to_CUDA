

import numpy as N


# given proportions, generate a randomized vector of ids that corresponds to each proportion
def proportion_sample(size, prop):

    sample_size_s = proportion_sample_size(total_size=size, prop=prop)

    ids = N.zeros(  size, dtype=N.int  )
    ids_len = 0

    for i in prop:
        sample_len = sample_size_s[i]

        ids[ ids_len : (ids_len + sample_len) ] = i
        ids_len += sample_len

        if ids_len >= len(ids):     break


    ids = N.random.permutation(ids)

    return ids

# generate sample sizes that equal to proportion, and total sample size roughly equal to total_size
def proportion_sample_size(total_size, prop):
    prop_sum = sum(prop[_] for _ in prop);          assert(prop_sum > 0)
    prop = {_:(prop[_] / prop_sum) for _ in prop}

    sample_size = {}
    for i in prop:
        sample_size[i] = N.int(N.round(total_size * prop[i]))
  
    return sample_size



def mad(data, axis=None):
    return N.median(   N.absolute(data - N.median(data, axis)), axis   )


