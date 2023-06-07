

'''
functions to handle tiff images

~/ln/tomominer/tomominer/image/tiff.py
'''



'''
from 
http://peterthomasweir.blogspot.com/2015/01/using-pylibtiff-to-read-tiff-files-in.html

'''


import numpy as np
from libtiff import TIFFfile
from libtiff import TIFF

def read(fileName):
    """
    Script to import tif file from imageJ,
    usage: zstack =  tiff.read(inFileName)
    PTW 2015/01/29
    """
    tiff = TIFFfile(fileName)
    
    subfile_type = tiff.get_subfile_types()
    assert len(subfile_type) == 1
    subfile_type = list(subfile_type)[0]

    samples, sample_names = tiff.get_samples(subfile_type=subfile_type)
    
    outList = []
    for sample in samples:
        outList.append(np.copy(sample)[...,np.newaxis])
        
    out = np.concatenate(outList,axis=-1)
    out = np.rollaxis(out,0,3)
    out = np.flipud(out)
    
    tiff.close()
    
    return out

def write_stack_color(fileName, array, compression=None):
    """
    Script to export tif file to imageJ,
    usage: tiff.write_stack(outFileName, array, compression=None)
    PTW 2015/04/22
    """
    #NEED: to document and decide about dealing with array with fewer dimensions than 4...
    outTiff = TIFF.open(fileName, mode='w')
    a = np.flipud(array)
    a = np.rollaxis(a, 3, 0)
    
    for zInd in range(a.shape[3]):
        outTiff.write_image(a[:,:,:,zInd],compression=compression, write_rgb=True)
        
    outTiff.close()
    
    return None






def write_stack_grey(fileName, array, compression=None):
    outTiff = TIFF.open(fileName, mode='w')

    a = array    
    for zInd in range(a.shape[2]):
        outTiff.write_image(a[:,:,zInd],compression=compression, write_rgb=True)
        
    outTiff.close()
    
    return None




