

# recursive search and export information needed for inspection or downstream analysis


import os
import sys
import pickle

import export as ex

if __name__ == '__main__':
    input_dir = sys.argv[1]
    output_dir_root = os.getcwd()
    
    for root, dirs, files in os.walk(input_dir):
        dump_file = os.path.join(root, 'dump.pickle')
        if not os.path.exists(dump_file):  continue

        with open(dump_file, 'rb') as f: dump = pickle.load(f)
        ex.export_dump(d=dump, out_dir=(output_dir_root + root) )

