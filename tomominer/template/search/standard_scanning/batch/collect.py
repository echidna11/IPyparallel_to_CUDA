#!/usr/bin/env python

'''
collect and summarize results


~/ln/tomominer/tomominer/template/search/standard_scanning/batch/collect.py

'''


import os, sys, json
import cPickle as pickle
import numpy as N

import tomominer.io.file as IF
import tomominer.image.vol.util as IVU

def main():
    with open('config_prepare__op.json') as f:      cop = json.load(f)
    with open(cop['config_stat_out_file']) as f:     cop_out = json.load(f)


    mrc = IF.read_mrc_header(cop['tomogram'])['MRC']
    siz = [mrc['nx'], mrc['ny'], mrc['nz']]


    c_max = N.zeros(siz) - N.inf
    phi_max = N.zeros(siz)
    theta_max = N.zeros(siz)
    psi_max = N.zeros(siz)

    for jl_i, jl in enumerate(cop_out['job_file_list']):
        with open(jl['job_file']) as f:     j = pickle.load(f)
        with open(j['stat_out']) as f:      r = json.load(f)

        se = j['job']['partition']['base']

        c = N.load(r['c'])
        c_full = N.zeros(siz) - N.inf
        IVU.paste_to_whole_map__se(c_full, c, se)
        
        ind = c_full > c_max
        c_max[ind] = c_full[ind]
        del c, c_full

        phi_t = N.load(r['phi'])
        phi_full = N.zeros(siz)
        IVU.paste_to_whole_map__se(phi_full, phi_t, se)
        phi_max[ind] = phi_full[ind]
        del phi_t, phi_full

        theta_t = N.load(r['theta'])
        theta_full = N.zeros(siz)
        IVU.paste_to_whole_map__se(theta_full, theta_t, se)
        theta_max[ind] = theta_full[ind]
        del theta_t, theta_full

        psi_t = N.load(r['psi'])
        psi_full = N.zeros(siz)
        IVU.paste_to_whole_map__se(psi_full, psi_t, se)
        psi_max[ind] = psi_full[ind]
        del psi_t, psi_full

        print '\r %d / %d        '%(jl_i, len(cop_out['job_file_list'])),
        sys.stdout.flush()



    st = {}
    out_dir = os.path.join(cop['out_dir'], 'collect')
    st['out_dir'] = out_dir
    if not os.path.isdir(out_dir):      os.makedirs(out_dir)

    st['c_mrc'] = os.path.join(out_dir, 'c.mrc')
    IF.put_mrc(c_max, st['c_mrc'])

   
    st['c'] = os.path.join(out_dir, 'c.npy')
    N.save(st['c'], c_max)

    st['phi'] = os.path.join(out_dir, 'phi.npy')
    N.save(st['phi'], phi_max)

    st['theta'] = os.path.join(out_dir, 'theta.npy')
    N.save(st['theta'], theta_max)

    st['psi'] = os.path.join(out_dir, 'psi.npy')
    N.save(st['psi'], psi_max)

    
    with open(cop['collect_stat_out_file'], 'w') as f:       json.dump(st, f, indent=2)


if __name__ == "__main__":
    main()

