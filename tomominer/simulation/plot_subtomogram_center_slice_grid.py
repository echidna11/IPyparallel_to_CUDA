#!/usr/bin/env python


'''
Plot center slice of a number of mrc files, 
'''

import os, json, random
import numpy as N
import tomominer.io.file as IF
import tomominer.filter.gaussian as FG
import imageio

def main():
    with open('plot_subtomogram_center_slice_grid_op.json') as f:     op = json.load(f)
    with open(op['input_file']) as f:   info = json.load(f) #info file containing infoformation about all the extracted subtomograms
    for plot_pdb in op["plot_pdb"]:
        if plot_pdb == 1:
            out_dir = os.path.join(op["out_dir"], "True_Labeled")
        elif plot_pdb == 2:
            out_dir = os.path.join(op["out_dir"], "Noise")
        else:
            out_dir = os.path.join(op["out_dir"], "Mix")
        if not os.path.isdir(out_dir):    os.makedirs(out_dir)
        # Read mrc size
        mrc_file = info[0]['subtomogram']
        mrc_size = IF.read_mrc_vol(mrc_file).astype(N.float).shape
        # Calculate grid size based on mrc size and number of subtomograms we want to show on grid
        grid_width = mrc_size[0]*op['width'] + (op['width']-1)*op['gap']
        grid_height = mrc_size[1]*op['height'] + (op['height']-1)*op['gap']
        # get list of tomograms ids
        tom_ids = set([info[_]['tomogram_id'] for _ in range(len(info))])
        data = {}
        ext = ""
        for i in tom_ids:
            print "Tomogram id: ", i
            # initialize empty grid for original and smooth verson of grid
            grid_v = N.zeros((grid_height,grid_width)).astype(N.float)
            grid_vs = N.zeros((grid_height,grid_width)).astype(N.float)
            # Read all subtomograms of particular tomogram id
            if plot_pdb==1:
                data[i] = [info[_]['subtomogram'] for _ in range(len(info)) if info[_]['tomogram_id']==i and 'pdb_id' in info[_]]
                ext = "_true_label.png"
            elif plot_pdb==2:
                data[i] = [info[_]['subtomogram'] for _ in range(len(info)) if info[_]['tomogram_id']==i and 'pdb_id' not in info[_]]
                ext = "_noise.png"
            else:
                data[i] = [info[_]['subtomogram'] for _ in range(len(info)) if info[_]['tomogram_id']==i]
                ext = ".png"
           # Select random sample out of it
            if len(data[i]) >= op['width']*op['height']:
                data[i] = random.sample(data[i], op['width']*op['height'])
            k = 0
            l = 0
            # Place each subtomogram on grid
            for j in range(len(data[i])):
                v = IF.read_mrc_vol(data[i][j]).astype(N.float)
                vs = FG.dog_smooth(v, 7)
                c = (N.array(v.shape) / 2.0).astype(N.int)
                grid_v[k:(k+v.shape[0]), l:(l+v.shape[1])] = v[:,:,c[2]]
                grid_vs[k:(k+vs.shape[0]), l:(l+vs.shape[1])] = vs[:,:,c[2]]
                #f_out = os.path.join(op['out_dir'], str(k)+"-"+str(l)+".png")
                #SM.imsave(f_out, v[:,:,c[2]])
                #save_png(N.squeeze(v[:,:,c[2]]), f_out)
                l += mrc_size[1] + op['gap']
                if l-op['gap'] == grid_v.shape[1]:
                    k += mrc_size[0] + op['gap']
                    l = 0
                if k-op['gap'] == grid_v.shape[0]:  break
            f_out = os.path.join(out_dir, op['output_grid_file']+ i + ext)
            fs_out = os.path.join(out_dir, op['output_grid_file']+ i + "_smooth" + ext)
            #save_png(N.squeeze(grid_v), f_out)
            imageio.imwrite(f_out, grid_v)
            imageio.imwrite(fs_out, grid_vs)

if __name__ == '__main__':
    main()


'''
related code
~/ln/tomominer/tomominer/image/vol/plot/mrc_center_slice.py
'''

