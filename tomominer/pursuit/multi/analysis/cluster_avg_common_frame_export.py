#!/usr/bin/env python


'''
export the cluster average in common frame
'''


import os, sys, json, shutil
import cPickle as pickle

def main():
    
    with open("export_cluster_averages_in_common_frame.json") as f:
        op = json.load(f)

    out_dir = os.path.abspath(op['out_dir'])

    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)

    os.makedirs(out_dir)

    fs_file = sys.argv[1]

    with open(fs_file) as f:        fs = json.load(f)

    fsp = fs['passes']
    fsp = {int(_):fsp[_] for _ in fsp}

    for pass_i in fsp:
        fspt = fsp[pass_i]
        with open(fspt['cluster_average_select_file'], 'rb') as f:   stcf = pickle.load(f)['selected_templates_common_frame']
        for c in stcf:
            shutil.copyfile(stcf[c]['subtomogram'], os.path.join(out_dir, '%03d-%03d-vol-avg.mrc'%(pass_i, c)))
            shutil.copyfile(stcf[c]['mask'], os.path.join(out_dir, '%03d-%03d-msk-avg.mrc'%(pass_i, c)))


    print 'exported to', out_dir

if __name__ == '__main__':
    main()

