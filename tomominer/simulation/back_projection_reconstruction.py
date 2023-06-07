#!/usr/bin/env python


'''
back projection reconstruction.

~/ln/tomominer/tomominer/simulation/back_projection_reconstruction.py
'''


import os
import random
from StringIO import StringIO
import pymatlab
import tempfile
import numpy as N

import tomominer.io.file as IF


# given a density map dm, call matlab routine to perform backprojection reconstruction
def do_reconstruction(session, dm, op, verbose=False):
    if verbose:     print 'do_reconstruction()', 'SNR', op['model']['SNR'], 'wedge_angle', op['model']['missing_wedge_angle']

    session.run('clear all;')


    session.putvalue('m', dm)

    ppc = param_prepare(op)
    session.run(ppc)
    #print ppc

    session.run('mbp = GenerateSimulationMap.backprojection_reconstruction(reconstruction_param, m, reconstruction_param.model.SNR);')

    mbp = session.getvalue('mbp')

    if ('inverse_intensity' in op) and op['inverse_intensity']:        mbp = -mbp

    return mbp


# when reconstructing large amount of subtomograms, the do_reconstruction() is very unstable, in such case, we try to use following function to reconstruct a bounch of subtomograms at a time, and the values are transferred through temporary files instead of using putvalue and get value (indicated by through_memory=False)
def do_reconstruction__batch(session, dms, op, through_memory=False, verbose=False):
    if verbose:     print 'do_reconstruction()', 'SNR', op['model']['SNR'], 'wedge_angle', op['model']['missing_wedge_angle']

    session.run('clear all;')

    session.run(param_prepare(op))

    re = {}
    for i in dms:
        if through_memory:
            session.putvalue('m', dms[i])
        else:
            f, fn = tempfile.mkstemp();     os.close(f)
            IF.put_mrc(dms[i], fn, overwrite=True)
            session.putvalue('fn', fn)
            session.run('im = tom_mrcread(fn);')
            session.run('m = im.Value;')
            os.remove(fn)

        session.run('mbp = GenerateSimulationMap.backprojection_reconstruction(reconstruction_param, m, reconstruction_param.model.SNR);')

        if ('inverse_intensity' in op) and op['inverse_intensity']:        session.run('mbp = -mbp')

        if through_memory:
            re[i] = session.getvalue('mbp')
        else:
            assert  not os.path.exists(fn)
            session.run('tom_mrcwrite(mbp, \'name\', fn);')
            re[i] = IF.read_mrc(fn)['value']
            os.remove(fn)

    return re



def param_prepare(op):
    out = StringIO()

    if 'random_seed' in op:
        print >>out, 'RandStream.setGlobalStream(RandStream(\'mt19937ar\',\'seed\', %d));'%(op['random_seed'])

    print >>out, 'reconstruction_param = struct();'

    print >>out, 'reconstruction_param.model = struct();'

    print >>out, 'reconstruction_param.model.missing_wedge_angle = %d;'%(op['model']['missing_wedge_angle'])

    print >>out, 'reconstruction_param.model.ctf = struct();'
    print >>out, 'reconstruction_param.model.ctf.pix_size = %f;'%(op['ctf']['pix_size'])      # pixel size in nm unit
    print >>out, 'reconstruction_param.model.ctf.Dz = %f;'%(op['ctf']['Dz'])
    print >>out, 'reconstruction_param.model.ctf.voltage = %f;'%(op['ctf']['voltage'])
    print >>out, 'reconstruction_param.model.ctf.Cs = %f;'%(op['ctf']['Cs'])
    if 'sigma' in op['ctf']:         print >>out, 'reconstruction_param.model.ctf.sigma = %f;'%(op['ctf']['sigma'])

    print >>out, 'reconstruction_param.model.do_reconstruction = true;'

    print >>out, 'reconstruction_param.model.backprojection_reconstruction = true;'
    print >>out, 'reconstruction_param.model.with_missing_wedge = true;'
    print >>out, 'reconstruction_param.model.titlt_angle_step = %d;'%(op['model']['titlt_angle_step'])
    print >>out, 'reconstruction_param.model.band_pass_filter = false;'
    print >>out, 'reconstruction_param.model.ctf.do_correction = false;'

    print >>out, 'reconstruction_param.model.SNR = %f;'%(op['model']['SNR'])

    print >>out, 'addpath(\'%s\');'%(op['matlab_code_path'])
    print >>out, 'addpath( genpath(\'%s\') );'%(op['matlab_tom_2008_code_path'])
    print >>out, 'addpath( genpath(\'%s\') );'%(op['matlab_tom_2005_code_path'])


    #print out.getvalue()        # print out commands for debugging purpose, so that you can paste it into matlab to see errors
    cmd = out.getvalue()

    out.close()

    return cmd



# just a test given all parameters fixed
def bp_test(v):


    print 'WARNING: if current version of matlab cannot be started, change to version 2011a instead'


    session = pymatlab.session_factory(options='-nodesktop -nodisplay')


    op = {'matlab_code_path':'/home/rcf-47/mxu/ln/frequent_structure/code', 'matlab_tom_2008_code_path':'/home/rcf-47/mxu/proj/imaging/electron/matlab_tom/2008/TOM_Release_2008', 'matlab_tom_2005_code_path':'/home/rcf-47/mxu/proj/imaging/electron/util/tom_2005', 'model':{'missing_wedge_angle':30, 'titlt_angle_step':1, 'SNR':N.nan}, 'ctf':{'pix_size':1.0, 'Dz':-15.0, 'voltage':300, 'Cs':2.2, 'sigma':0.4}}

    vb = do_reconstruction(session, v, op, verbose=False)
    return vb



if __name__ == '__main__':

    class Object:
        pass

    self = Object()
    self.matlab_session = pymatlab.session_factory(options='-nodesktop -nodisplay')

    self.matlab_session.run( 'addpath(\'%s\')'%(op['matlab_code_path']) )
    self.matlab_session.run( 'addpath(\'%s\')'%(op['tom_toolbox_2008_path']) )
    self.matlab_session.run( 'addpath(\'%s\')'%(op['tom_toolbox_2005_path']) )



