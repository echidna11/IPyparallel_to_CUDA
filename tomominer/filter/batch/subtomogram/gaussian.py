#!/usr/bin/env python


'''
given a list of subtomograms, filter each one, and generate a new list

~/ln/tomominer/tomominer/filtering/batch/subtomogram/gaussian.py

'''

import os, json, uuid

import tomominer.io.file as IF
import tomominer.filter.gaussian as FG


def main():

    with open('filtering_batch_subtomogram_gaussian__op.json') as f:        op = json.load(f)

    with open(op['data in']) as f:      dj = json.load(f)

    op['out dir'] = os.path.abspath(op['out dir'])

    assert not os.path.isdir(op['out dir'])
    os.makedirs(op['out dir'])

    dj_new = []
    for i, d in enumerate(dj):
        v = IF.read_mrc_vol(d['subtomogram'])

        v = FG.smooth(v, sigma=op['sigma'])

        out_file = os.path.join(op['out dir'], 'st--%05d--%s.mrc'%(i, str(uuid.uuid4())))
        IF.put_mrc(v, out_file)

        dj_new.append({'subtomogram':out_file, 'mask':d['mask'], 'original':{'subtomogram':d['subtomogram']}})


    with open(op['data out'], 'w') as f:        json.dump(dj_new, f, indent=2)



if __name__ == "__main__":
    main()


