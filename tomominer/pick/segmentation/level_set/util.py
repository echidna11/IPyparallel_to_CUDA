# utility functions for segmentation based particle picking

import tomominer.segmentation.active_contour.chan_vese.segment as SACS

def segment(v, op):
    phi = SACS.segment_with_postprocessing(v, op)

    return phi

