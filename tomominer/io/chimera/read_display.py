#!/usr/bin/env python


# use chimera to load a mrc file, then display the volume slice by slice to see if the map is correctly loaded in a qualitative way. For visual check can compare imod display result.

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Qt4Agg')
 
    import sys
    mrc_file = sys.argv[1]

    print 'loading map...'
    import tomominer.io.file as IF
    v = IF.read_mrc__chimera(mrc_file)

    if len(sys.argv) > 2:
        transpose_mode = int(sys.argv[2])
    else:
        transpose_mode = 0

    import numpy as N
    if transpose_mode == 0:
        pass
    elif transpose_mode == 1:
        v = N.transpose(v, axes=[2,0,1])
    elif transpose_mode == 2:
        v = N.transpose(v, axes=[1,2,0])
    else:
        raise Exception('transpose_mode')

    #v = N.rot90(v, k=2)     # rotate so that the displayed image has same orientation as imod's display


    import tomominer.plot.marker_display as PM
    PM.display(v)


