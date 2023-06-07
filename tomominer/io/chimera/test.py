

# commands for testing chimera routines


# result: the chimera python libraries loads a number of dynamic libraries. So we need to use chimera's own python2.7 to correctly load these dynamical libraries

import sys
sys.path.append('/auto/cmb-panasas2/mxu/proj/imaging/electron/util/chimera/1.8/share')
sys.path.append('/auto/cmb-panasas2/mxu/proj/imaging/electron/util/chimera/1.8/lib')


mrc_file = '/home/rcf-47/mxu/ln/frequent_structure/data/out/analysis/beck/u2os/classification/tomograms/original/9.rec'

import VolumeData.mrc.mrc_format as VMMF
header = VMMF.MRC_Data(mrc_file, 'mrc')


import VolumeData.mrc.mrc_grid as VMMG
im = VMMG.MRC_Grid(mrc_file)

v = im.read_matrix(ijk_origin=[0,0,0], ijk_size=[100,100,100], ijk_step=[1,1,1], progress=None)


v = im.read_matrix(ijk_origin=[0,0,0], ijk_size=im.mrc_data.matrix_size, ijk_step=[1,1,1], progress=None)




# test through which permutation can convert the read_matrix into correct order, in comparison to matlab version

import numpy as N
import tomominer.io.file as IF
vm = IF.read_mrc__matlab(       mrc_file,   sub_region=N.array([[0,400], [0,400], [0,400]])     )
vc_full = IF.read_mrc__chimera(mrc_file)
N.transpose(vc_full, axes=[2,1,0]).shape

N.where(vm==vm.max())

w = N.array(        N.where(vc_full==vm[375,120,377])       )
w[:, N.logical_and(w[0,:]==377, w[1,:]==120)]

# the result of above shows that after transpose, there is still a small shift of 5 voxels, between vc_full and vm, dont know why



vc = N.transpose(vc_full, axes=[1,2,0])[:100,:100,:100]
vm - vc




