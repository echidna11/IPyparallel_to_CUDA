#!/usr/bin/env python


'''
a simple command line program that uses fastaniso to perform anistropic diffusion filtering
'''


import os, json

import tomominer.io.file as IF


if __name__ == '__main__':

    with open('filtering_anistropic_diffusion_fastaniso_cli__op.json') as f:       op = json.load(f)

    v = IF.read_mrc_vol(op['file in'])

    import fastaniso
    v = fastaniso.anisodiff3(v, niter=op['iter'], kappa=op['kappa'], gamma=op['gamma'], step=tuple(op['step']), option=op['option'])

    IF.put_mrc(v, op['file out'], overwrite=True)


