#!/usr/bin/env python

# collect all subtomogram files and corresponding wedge files, create a json file for classification

import os
import sys
import json

import random as rd



def process_dir(cur_dir, wedge_mask_file = None):
    json_data = []

    if wedge_mask_file is None:
        files = [f for f in os.listdir(cur_dir) if (os.path.isfile(os.path.join(cur_dir, f)) and f.startswith('wedge-mask') and f.endswith('.mrc'))]
	if len(files) > 1:
	    raise RuntimeError('too many wedge-mask files in a single dir');

	if len(files) == 1:
	    wedge_mask_file = os.path.join(cur_dir, files[0])

    files = [f for f in os.listdir(cur_dir) if (os.path.isdir(os.path.join(cur_dir, f)) and (not f.startswith('_')))]       # those dirs starts with '_' will be ignored
    for f in files:
        json_data.extend( process_dir(os.path.join(cur_dir, f), wedge_mask_file = wedge_mask_file) )


    if wedge_mask_file is None:
	return json_data

    files = [f for f in os.listdir(cur_dir) if (os.path.isfile(os.path.join(cur_dir, f)) and (not f.startswith('wedge-mask')) and (not f.startswith('_')) and f.endswith('.mrc'))]

    for f in files:
        json_data.append({'subtomogram':os.path.join(cur_dir, f), 'mask':wedge_mask_file, 'angle':[rd.random()*6.28, rd.random()*6.28,rd.random()*6.28], 'loc':[0.0, 0.0, 0.0]})    # generate random initial orientations to reduce bias introduced by same oriented missing wedge

    return json_data

if __name__ == '__main__':

    json_data = process_dir(os.getenv('PWD'))       # os.getcwd() will return physical location, but os.getenv('PWD') will use symbolic link
    json.dump(json_data, sys.stdout, indent=2)
    sys.stdout.flush()


