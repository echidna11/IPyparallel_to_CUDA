#!/usr/bin/env python



'''
identify the lines that contain comments

~/ln/tomominer/tomominer/prog/code_export/comment_identify_remove.py
'''

import os, sys, astor
from tomominer.prog.code_export.export_from_profile import export_code


def main(file_in, file_out, header=None):
    print 'loading', file_in, 'and exporting', file_out
    ast = astor.parsefile(file_in)
    src = astor.to_source(ast)
    export_code(src=src, file_out=file_out, executable=True)


if __name__ == '__main__':
    file_in = sys.argv[1]
    file_out = sys.argv[2]

    if len(sys.argv) >= 4:
        # read header
        with open(sys.argv[3]) as f:
            header = f.read()
    else:
        header = None

    main(file_in, file_out, header=header)

