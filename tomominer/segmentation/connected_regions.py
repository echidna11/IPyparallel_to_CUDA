
import numpy as N
import tomominer.core as core


def connected_regions(msk):
    msk = N.array(msk, dtype=N.int32, order='F')
    return core.connected_regions(msk)

