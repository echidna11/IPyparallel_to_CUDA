#!/usr/bin/env python


# given classification iterations, list Fourier Shell Correlation / resolution information for each cluster



if __name__ == '__main__':

    import os
    import pickle

    import tomominer.statistics.vol as SV

    pass_i = 0
    while True:
        cpf  = os.path.join(os.getcwd(), 'pass_%03d' % (pass_i), 'cluster_ssnr_fsc.pickle')
        if not os.path.exists(cpf) :       break

        with open(cpf) as f:   cp = pickle.load(f)

        print pass_i,
        for c in cp:
            print '\t %d %d %.1f'%(c, SV.resolution_ind(cp[c]['fsc']), cp[c]['fsc'].sum()),

        print

        pass_i += 1

