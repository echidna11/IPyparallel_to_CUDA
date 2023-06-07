#!/usr/bin/env python


import os, sys, json, shutil

def main():
    stat_file = sys.argv[1]
    dest_dir = sys.argv[2]

    assert      not os.path.isdir(dest_dir)
    os.makedirs(dest_dir)

    with open(stat_file) as f:      stat = json.load(f)

    for rep_i, stat_ in stat.iteritems():
        rep_i = int(rep_i)
        src_dir = os.path.join(stat_['pass_dir'], 'common_frame')

        for c, fsc_stat in stat_['fsc_stat'].iteritems():
            c = int(c)

            src_msk_file = os.path.join(src_dir, 'clus_mask_avg_%03d.mrc'%(c,))
            dest_msk_file = os.path.join(dest_dir, '%03d-%03d-msk.mrc'%(rep_i,c))           ;       assert  not os.path.exists(dest_msk_file)
            
            src_vol_file = os.path.join(src_dir, 'clus_vol_avg_%03d.mrc'%(c,))
            dest_vol_file = os.path.join(dest_dir, '%03d-%03d-vol.mrc'%(rep_i,c))           ;       assert  not os.path.exists(dest_vol_file)

            shutil.copyfile(src_msk_file, dest_msk_file)
            shutil.copyfile(src_vol_file, dest_vol_file)




if __name__ == "__main__":
    
    main()

