

'''
convert image format

~/ln/tomominer/tomominer/image/format_convert.py
'''

import os

def convert_a2ping(in_file, out_file, resolution=72):

    assert not os.path.isfile(out_file)

    args = ['a2ping', '--Resolution=%d'%(resolution,), in_file, out_file]

    import subprocess
    subprocess.call(args)

    assert os.path.isfile(out_file)



