import os, json, copy, gc
import numpy as N
#import matplotlib.pyplot as plt
import tifffile as tiff
import tomominer.io.file as IF

def plot_single_channel(x, c, title, f_name):
    plt.imshow(N.rot90(x[c, :, :, x.shape[3]/2]), cmap="gray")
    plt.title(title)
    plt.xlabel("X axis")
    plt.ylabel("Y axis")
    plt.xticks([])
    plt.yticks([])
    plt.savefig(f_name)
    plt.close()

def plot_all_channels(x, channels, outdir):
    if not os.path.isdir(outdir):   os.makedirs(outdir)
    for i in sorted(channels.keys()):
        f_name = os.path.join(outdir, i + ".png")
        plot_single_channel(x, int(i), channels[i], f_name)

def extract_coordinates(im_dna_i):
    start = []
    end = []
    for i in range(im_dna_i.shape[0]):
        if im_dna_i[i, :, :].sum()!=0:
            start.append(i)
            break
    for i in range(start[0]+1, im_dna_i.shape[0]):
        if im_dna_i[i, :, :].sum()==0:
            end.append(i-1)
            break
    for i in range(im_dna_i.shape[1]):
        if im_dna_i[:, i, :].sum()!=0:
            start.append(i)
            break
    for i in range(start[1]+1, im_dna_i.shape[1]):
        if im_dna_i[:, i, :].sum()==0:
            end.append(i-1)
            break
    for i in range(im_dna_i.shape[2]):
        if im_dna_i[:, :, i].sum()!=0:
            start.append(i)
            break
    for i in range(start[2]+1, im_dna_i.shape[2]):
        if im_dna_i[:, :, i].sum()==0:
            end.append(i-1)
            break
    return start, end

def extract_cells(im, outdir, tmp_file, file_i):
    if not os.path.isdir(outdir):   os.makedirs(outdir)
    im_dna = copy.deepcopy(im[6, :, :, :]) # Segmented DNA
    values = N.trim_zeros(sorted(list(set(im_dna.flatten()))))
    file_name = outdir[outdir.find("/")+1:]
    print file_name, "have", len(values), "segmented cells"
    if len(values) == 0:
      print file_name, "doesn't have any segmented cells"
      return 0
    
    max_value = max(values) + 1
    for cell in range(len(values)):
        print "Saving mrcs for cell:", cell, "in", file_name
        f_out = [os.path.join(outdir, str(cell)+"_seg_dna.mrc"), os.path.join(outdir, str(cell)+"_nuc_contour.mrc"), os.path.join(outdir, str(cell)+"_seg_EGFP.mrc")]
        flag = 0
        for f_o in range(len(f_out)): # skip processing only if all three files are available
            if not os.path.isfile(f_out[f_o]):
                flag = 1
                break
        if flag == 1:
            # Segmented Nucleus Mask
            im_dna_i = copy.deepcopy(im_dna)
            # place 0 and 1 in array
            N.place(im_dna_i, im_dna_i!=values[cell], [0])
            N.place(im_dna_i, im_dna_i==values[cell], [1])
            assert sorted(list(set(im_dna_i.flatten()))) == [0,1]
            # extract coordinates to crop the cell from full field matrix
            start, end = extract_coordinates(im_dna_i)
            # crop the matrix
            im_dna_i = im_dna_i[start[0]-1:end[0]+2, start[1]-1:end[1]+2, start[2]-1:end[2]+2]
            # save the mrc for segmented dna
            IF.put_mrc(im_dna_i, f_out[0])
            # change the pixel values of mrc using alterheader command
            tmp_file_out = "alterheader_tmp_output_" + file_i[:file_i.find(".")] + ".txt"
            tmp = "alterheader " + f_out[0] + " < " + tmp_file + " > " + tmp_file_out
            tmp = os.system(tmp)
            mask = N.ma.masked_where(im_dna_i!=0, im_dna_i)
            mask = mask.mask
            del im_dna_i
            gc.collect()

            # Nucleus contour
            im_contour = copy.deepcopy(im[8, :, :, :])
            # Crop and mask the contour with segmented dna array
            im_contour = im_contour[start[0]-1:end[0]+2, start[1]-1:end[1]+2, start[2]-1:end[2]+2]
            im_contour[~mask]=0
            N.place(im_contour, im_contour!=0, [1])
            # save the nucleus contour
            IF.put_mrc(im_contour, f_out[1])
            # change the pixel values of mrc using alterheader command
            tmp = "alterheader " + f_out[1] + " < " + tmp_file + " > " + tmp_file_out
            tmp = os.system(tmp)
            del im_contour 
            gc.collect()

            # Segmented EGFP
            im_egfp = copy.deepcopy(im[4, :, :, :])
            # Crop and mask the contour with segmented dna array
            im_egfp = im_egfp[start[0]-1:end[0]+2, start[1]-1:end[1]+2, start[2]-1:end[2]+2]
            im_egfp[~mask]=0
            N.place(im_egfp, im_egfp!=0, [1])
            # save the nucleus contour
            IF.put_mrc(im_egfp, f_out[2])
            # change the pixel values of mrc using alterheader command
            tmp = "alterheader " + f_out[2] + " < " + tmp_file + " > " + tmp_file_out
            tmp = os.system(tmp)
            del start, end, tmp_file_out, tmp, mask, im_egfp
            gc.collect()

        if flag == 0: # risk cover. In case they were saved but header not updated
            tmp_file_out = "alterheader_tmp_output_" + file_i[:file_i.find(".")] + ".txt"
            tmp = "alterheader " + f_out[0] + " < " + tmp_file + " > " + tmp_file_out
            tmp = os.system(tmp)
            tmp = "alterheader " + f_out[1] + " < " + tmp_file + " > " + tmp_file_out
            tmp = os.system(tmp)
            tmp = "alterheader " + f_out[2] + " < " + tmp_file + " > " + tmp_file_out
            tmp = os.system(tmp)
            del tmp_file_out, tmp, f_out
            gc.collect()
    
    count = len(values)
    del values, im_dna, file_name, max_value
    gc.collect()
    return count

def pool_function(file_i, op):
    os.chdir(op["working_dir"])
    im = tiff.TiffFile(os.path.join(op["input_dir"], file_i)) # read ome.tif file
    assert im.is_ome
    imm = im.ome_metadata # get metadata
    imm = imm['Image']['Pixels']
    # number of pixels in each direction
    sizeX = imm['SizeX']
    sizeY = imm['SizeY']
    sizeZ = imm['SizeZ']
    # real pixel distance in each direction
    psizeX = imm['PhysicalSizeX']
    psizeY = imm['PhysicalSizeY']
    psizeZ = imm['PhysicalSizeZ']
    # unit of real distance
    unitX = imm['PhysicalSizeXUnit']
    unitY = imm['PhysicalSizeYUnit']
    unitZ = imm['PhysicalSizeZUnit']
    # Dimension order
    dimOrder = imm['DimensionOrder']
    channels = {_['ID'][_['ID'].find(":")+3:]:_['Name'] for _ in imm['Channel']}
    im = im.asarray() # Read actual image data
    # the dimension order is ZCYX, make it CXYZ
    im = N.swapaxes(im, 0,1)
    im = N.swapaxes(im, 1,3)
    # plot all channels just to see
    #outdir = os.path.join(op["plot_dir"], file_i[:file_i.find(".")])
    #plot_all_channels(im, channels, outdir)

    outdir = os.path.join(op["mrc_dir"], file_i[:file_i.find(".")])
    tmp_txt = "del\n" + str(N.round(psizeX*10000, decimals=2)) + " " + str(N.round(psizeY*10000, decimals=2)) + " " + str(N.round(psizeZ*10000, decimals=2)) + "\ndone"
    tmp_file = "alterheader_input_" + file_i[:file_i.find(".")] + ".txt"
    tmp_f = open(tmp_file, "w")
    tmp_f.write(tmp_txt)
    tmp_f.close()
    while not os.path.isfile(tmp_file):
        flag = 0
    count = extract_cells(im, outdir, tmp_file, file_i)
    del im, imm, sizeX, sizeY, sizeZ, psizeX, psizeY, psizeZ, unitX, unitY, unitZ, dimOrder, channels, outdir, tmp_txt, tmp_file, tmp_f, flag
    gc.collect()
    return count

def main():
    return

if __name__=="__main__":
    main()

