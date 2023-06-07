#!/usr/bin/env python



'''
concat multiple surfaces in amira format


~/ln/tomominer/tomominer/geometry/surface/amira/concat_surfaces.py

'''

import sys

import tomominer.geometry.surface.amira.util as GSAU
import tomominer.geometry.surface.util as GSU

if __name__ == '__main__':

    s = None
    for i in range(len(sys.argv)):
        if i < 1:   continue

        with open(sys.argv[i]) as f:    st = GSAU.surf_parse(f)

        s = GSU.concat(s, st)


    GSAU.export_surf_ascii_simple(s, sys.stdout)

