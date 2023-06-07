#!/usr/bin/env python


'''
display header information


~/ln/tomominer/tomominer/image/vol/format/mrc/disp_header_info.py
'''

def main():
    import sys
    import tomominer.io.file as IF
    im  = IF.read_mrc(sys.argv[1], read_data=False)
    mrc = im['header']['MRC']
    print 'xorg yorg zorg', mrc['xorg'], mrc['yorg'], mrc['zorg']
    print 'nx ny nz', mrc['nx'], mrc['ny'], mrc['nz']
    print 'xlen ylen zlen', mrc['xlen'], mrc['ylen'], mrc['zlen']

if __name__ == '__main__':
    main()

