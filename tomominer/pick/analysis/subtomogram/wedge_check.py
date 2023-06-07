#!/usr/bin/env python


'''
given a data json file, average the subtomograms, and perform fft transform, then save magnitute so that we can use it to check the correctness of wedge mask

pick/analysis/subtomogram/wedge_check
~/ln/tomominer/tomominer/pick/analysis/subtomogram/wedge_check.py

related code
/home/rcf-47/mxu/ln/tomominer/tomominer/image/vol/wedge/tilt_angle_inference.py

'''


import os, sys, json, random
import numpy as N
from numpy.fft import fftshift, fftn

import tomominer.io.file as IF
import tomominer.statistics.util as SU

if __name__ == '__main__':

    with open('pick_analysis_subtomogram_wedge_check__op.json') as f:       op = json.load(f)
    
    with open(op['data json file']) as f:       dj = json.load(f)

    # collect subtomograms for different wedge masks
    dc = {}
    for d in dj:
        if d['mask'] not in dc:     dc[d['mask']] = []

        dc[d['mask']].append(d['subtomogram'])


    if not os.path.isdir(op['out dir']):        os.makedirs(op['out dir'])

    # make averages
    avg = {}
    for mk in dc:

        dt = dc[mk]
        if ('sample num' in op) and (len(dt) > op['sample num']):      dt = random.sample(dt, op['sample num'])

        vs = None
        for i, s in enumerate(dt):
            print '\r', '%0.4f'%(float(i) / len(dt),),        ;       sys.stdout.flush()

            v = IF.read_mrc_vol(s)
            if vs is None:
                vs = v.astype(N.float)
            else:
                vs += v

        avg[mk] = vs / len(dt)




    i = 0
    for mk in avg:
        vf = N.log(     N.abs(fftshift(fftn(avg[mk])))      )
        IF.put_mrc(vf, os.path.join(op['out dir'], '%03d-vol-f.mrc'%(i,)))

        if os.path.isfile(mk):
            m = IF.read_mrc_vol(mk)
            IF.put_mrc(m, os.path.join(op['out dir'], '%03d-mask.mrc'%(i,)))

        i += 1


