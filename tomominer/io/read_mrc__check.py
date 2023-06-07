#!/usr/bin/env python



# compare the matlab and python version of read_mrc to see if the volumeric data loaded are consistent


import numpy as N
import tomominer.io.file as IF

# check the consistancy of read_mrc() and the tom_mrcread() in tom_toolbox
def read_mrc__check(mrc_file, sub_region=None, show_progress=False, matlab_session=None):

    print 'read_mrc__check()'

    if sub_region is not None:          sub_region = N.array(sub_region)

    try:
        vt = IF.read_mrc__matlab(mrc_file, sub_region=sub_region, show_progress=show_progress, session=matlab_session)
    except:
        print 'read error'
        return

    i = IF.read_mrc(mrc_file, show_progress=show_progress)
    v = i['value']

    if sub_region is not None:     v = v[sub_region[0,0]:sub_region[0,1], sub_region[1,0]:sub_region[1,1], sub_region[2,0]:sub_region[2,1]]

    print 'size', v.shape, 'min', v.min(), 'max', v.max(), 'mean', v.mean(), 'std', v.std(), 'std32', v.astype(N.float32).std(), 'std64', v.astype(N.float64).std()
    print 'size', vt.shape, 'min', vt.min(), 'max', vt.max(), 'mean', vt.mean(), 'std', vt.std(), 'std32', vt.astype(N.float32).std(), 'std64', vt.astype(N.float64).std()


    dif = abs(v - vt).max()     # this number should be 0, if the volumeric data can be correctly read
    print 'difference to the matlab loaded version', dif, '----- if this value is 0, then the two versions of mrc reader is consistent'


    del v
    del vt

    return dif



if __name__ == '__main__':

    matlab_tom_2008_code_path = '/home/rcf-47/mxu/proj/imaging/electron/matlab_tom/2008/TOM_Release_2008'

    import pymatlab
    session = pymatlab.session_factory(options='-nodesktop -nodisplay')
    session.run( 'addpath(  genpath(\'%s\') )'%(matlab_tom_2008_code_path,) )

    import os

    for root, dirs, files in os.walk('.'):
        for file_t in files:
            if (not file_t.lower().endswith('.mrc')) and (not file_t.lower().endswith('.rec')):      continue
            
            file_path = os.path.join(root, file_t)
            file_path = os.path.abspath(file_path)

            print 'checking', file_path
            read_mrc__check(file_path, sub_region=[[0,200], [0,200], [0,200]], matlab_session=session)

