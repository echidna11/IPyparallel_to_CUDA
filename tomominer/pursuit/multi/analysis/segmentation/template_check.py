#!/usr/bin/env python


# combined the level set phi and the template for visual check if the segmentation is successful. Store the output to a temporory folder: /tmp/seg

if __name__ == '__main__':

    import json
    with open('pursuit-op.json') as f: op = json.load(f)
    op = op['segmentation']


    avg_prefix = 'clus_vol_avg_'
    phi_prefix = 'clus_vol_seg_phi_'
    out_dir_root = '/tmp/seg'


    import os
    import shutil
    import numpy as N
    import tomominer.io.file as IF

    for root, dirs, fs in os.walk("."):
        for f in fs:
            if not(f.startswith(phi_prefix) and f.endswith('.mrc')):  continue

            out_dir = os.path.join(out_dir_root, root)
            if not os.path.isdir(out_dir):  os.makedirs(out_dir)
            
            phi_f = os.path.join(root, f)

            # make a copy of template file
            avg_f = os.path.join(root, f.replace(phi_prefix, avg_prefix))
            shutil.copyfile(avg_f, os.path.join(out_dir_root, avg_f))

            phi = IF.read_mrc(phi_f)['value']
            v = IF.read_mrc(avg_f)['value']

            ind = phi > (phi.max() * op['phi_propotion_cutoff'])
            v[N.logical_not(ind)] = float('nan')    # v[ind].mean()

            IF.put_mrc(v, os.path.join(out_dir_root, avg_f.replace('.mrc', '-seg.mrc')), overwrite=True)

            if False:
                # make a copy of level set file
                shutil.copyfile(phi_f, os.path.join(out_dir_root, phi_f))
            else:
                IF.put_mrc(phi, os.path.join(out_dir_root, avg_f.replace('.mrc', '-phi.mrc')), overwrite=True)

            print f


