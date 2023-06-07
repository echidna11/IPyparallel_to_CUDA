#!/usr/bin/env python


'''

manually extract one subtomomogram


~/ln/tomominer/tomominer/pick/subtomogram/manual_extract.py

'''


import json
import tomominer.io.file as IF


def main():

    with open('pick_subtomogram_manual_extract__op.json') as f:     op = json.load(f)

    v = IF.read_mrc_vol(op['tomogram'])

    vc = CV.cut_from_whole_map(whole_map=v, c=N.array(op['loc']), siz=op['size'])

    IF.put_mrc(op['subtomogram file out'])


if __name__ == '__main__':
    main()


