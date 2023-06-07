#!/usr/bin/env python



'''
given a list of references, prepare mirrors, and put into a new list (together with original)

~/ln/tomominer/tomominer/pick/reference/prepare_mirrors.py

'''


import os, json, shutil

import tomominer.io.file as IF
import tomominer.filter.gaussian as FG


def main():
    with open('pick_reference_prepare_mirrors__op.json') as f:    op = json.load(f)

    with open(op['reference list in']) as f:    dj = json.load(f)

    op['out dir'] = os.path.abspath(op['out dir'])

    if os.path.isdir(op['out dir']):        shutil.rmtree(op['out dir'])
    os.makedirs(op['out dir'])

    dj_new = []

    for i, d in enumerate(dj):
        v = IF.read_mrc_vol(d['subtomogram'])

        if 'gaussian sigma' in op:      v = FG.smooth(v=v, sigma=op['gaussian sigma'])

        mrc_out = os.path.join(op['out dir'], '%03d.mrc'%(i,))
        assert  not os.path.isfile(mrc_out)
        IF.put_mrc(v, mrc_out, overwrite=True)

        dt = {'subtomogram':mrc_out, 'mask':d['mask'], 'original':d, 'mirror_dim':-1}
        dj_new.append(dt)

        for dim_i in range(3):
            if dim_i == 0:
                vt = v[::-1, :, :]
            elif dim_i == 1:
                vt = v[:, ::-1, :]
            elif dim_i == 2:
                vt = v[:, :, ::-1]
            else:
                assert False

            mrc_out = os.path.join(op['out dir'], '%03d-%d.mrc'%(i, dim_i))
            assert  not os.path.isfile(mrc_out)

            IF.put_mrc(vt, mrc_out, overwrite=True)

            dt = {'subtomogram':mrc_out, 'mask':d['mask'], 'original':d, 'mirror_dim':dim_i}
            dj_new.append(dt)

    with open(op['mirror list out'], 'w') as f:       json.dump(dj_new, f, indent=2)


if __name__ == '__main__':
    main()


