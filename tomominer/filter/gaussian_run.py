#!/usr/bin/env python



'''

command line gaussian filtering

~/ln/tomominer/tomominer/filtering/gaussian_run.py

'''


if __name__ == '__main__':

    import json
    with open('gaussian_run__op.json') as f:   op = json.load(f)

    import os
    assert(not os.path.exists(op['file out']))  # just a little protection of overwriting existing file

    import tomominer.io.file as IOF

    v_i = IOF.read_mrc_vol(op['file in'])


    sigma = float(op['sigma'])
    if 'voxel size' in op:      sigma /= float(op['voxel size'])

    if False:
        import tomominer.filter.gaussian as FG
        v_o = FG.smooth(v_i, sigma)
    else:
        import scipy.ndimage as SN
        v_o = SN.gaussian_filter(input=v_i, sigma=sigma)

    import numpy as N
    v_o = N.array(v_o, order='F')
    
    IOF.write_mrc_by_chunk(v=v_o, path=op['file out'])




