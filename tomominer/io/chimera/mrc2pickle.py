#!/usr/bin/env python


if __name__=='__main__':

    import os
    import sys

    mrc_file = sys.argv[1]
    pickle_file_root = sys.argv[2]


    from mrc_format_ext import MRC_Data_ext
    im = MRC_Data_ext(mrc_file, 'mrc')
    print 'volume header', im.matrix_size,

    v = im.read_matrix(ijk_origin=[0,0,0], ijk_size=im.matrix_size, ijk_step=[1,1,1], progress=None)

    import numpy as N
    v = N.transpose(v, axes=[2,1,0])            # for some strange reason, we must perform this transpose so that the size of v is consistent with im.matrix_size
    assert      N.all(      N.array(im.matrix_size) == N.array(v.shape)     )

    print 'volume read', 'type', v.dtype

    import pickle
    with open(pickle_file_root+'-header.pickle', 'wb') as f:      pickle.dump(im, f, protocol=-1)

    N.save(pickle_file_root+'-vol.npy', v)

'''
#Usage example

~/ln/tomominer/tomominer/input_output/chimera/mrc2pickle.sh /home/rcf-47/mxu/ln/frequent_structure/data/out/analysis/beck/u2os/classification/tomograms/original/9.rec /tmp/t

'''



'''
# tests for transpose data to make sure that the loaded volume is in correct order


'''



