#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
compress current code

example usage
~/ln/tomominer/tomominer/common/zip_code.py ./ ~/ln/tomominer/tomominer
'''

import os
import sys
import time
from file_util import zipdir

if __name__ == '__main__':

    dest_dir = sys.argv[1]
    src_dir = os.path.realpath(sys.argv[2])
    
    zip_file = os.path.join(dest_dir, 'code-bak-%d.zip'%(int(time.time())))
    
    zipdir(zipFilePath=zip_file, dirPath=src_dir)
    
    
