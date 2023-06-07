
# command line alignment of subtomogram 2 against subtomogram 1

import os
import sys

if __name__ == '__main__':
    
    sys.path.append( os.path.dirname( os.path.dirname(__file__) ))
    
    import tomo

    v1f = sys.argv[1]
    m1f = sys.argv[2]
    v2f = sys.argv[3]
    m2f = sys.argv[4]

    out_path = sys.argv[5]

    import tomominer.io.file as fio 
    v1 = fio.get_mrc(v1f)
    m1 = fio.get_mrc(m1f)
    v2 = fio.get_mrc(v2f)
    m2 = fio.get_mrc(m2f)


    import util
    re = util.align_vols(v1, m1, v2, m2, L=36)
    
    v2r = tomo.rotate_vol_pad_mean_py(v2, re['ang'], re['loc'])


    fio.put_mrc(v2r, out_path)


