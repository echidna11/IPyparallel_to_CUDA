'''
same as back_projection_reconstruction.py, but use octave instad of matlab


~/ln/tomominer/tomominer/simulation/back_projection_reconstruction__octave.py

'''


import os
import random
from StringIO import StringIO
import tempfile
import numpy as N
import oct2py

import tomominer.io.file as IF


# given a density map dm, call matlab routine to perform backprojection reconstruction
def do_reconstruction(session, dm, op, verbose=False):
    if verbose:     print 'do_reconstruction()', 'SNR', op['model']['SNR'], 'wedge_angle', op['model']['missing_wedge_angle']

    #session.eval('clear all;')


    session.push('m', dm)

    ppc = param_prepare(op)
    session.eval(ppc)
    #print ppc

    session.eval('mbp = backprojection_reconstruction(reconstruction_param, m, reconstruction_param.model.SNR);')

    mbp = session.pull('mbp')

    if ('inverse_intensity' in op) and op['inverse_intensity']:        mbp = -mbp

    return mbp


# when reconstructing large amount of subtomograms, the do_reconstruction() is very unstable, in such case, we try to use following function to reconstruct a bounch of subtomograms at a time, and the values are transferred through temporary files instead of using putvalue and get value (indicated by through_memory=False)
def do_reconstruction__batch(session, dms, op, through_memory=False, verbose=False):
    if verbose:     print 'do_reconstruction()', 'SNR', op['model']['SNR'], 'wedge_angle', op['model']['missing_wedge_angle']

    session.eval('clear all;')

    session.eval(param_prepare(op))

    re = {}
    for i in dms:
        if through_memory:
            session.push('m', dms[i])
        else:
            f, fn = tempfile.mkstemp();     os.close(f)
            IF.put_mrc(dms[i], fn, overwrite=True)
            session.push('fn', fn)
            session.eval('im = tom_mrcread(fn);')
            session.eval('m = im.Value;')
            os.remove(fn)

        session.eval('mbp = backprojection_reconstruction(reconstruction_param, m, reconstruction_param.model.SNR);')

        if ('inverse_intensity' in op) and op['inverse_intensity']:        session.eval('mbp = -mbp')

        if through_memory:
            re[i] = session.pull('mbp')
        else:
            assert  not os.path.exists(fn)
            session.eval('tom_mrcwrite(mbp, \'name\', fn);')
            re[i] = IF.read_mrc(fn)['value']
            os.remove(fn)

    return re



def param_prepare(op):
    out = StringIO()

    if 'random_seed' in op:
        print >>out, 'rand(\'seed\', %d);'%(op['random_seed'])

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


    #print out.pull()        # print out commands for debugging purpose, so that you can paste it into matlab to see errors
    cmd = out.getvalue()

    out.close()

    return cmd



# just a test given all parameters fixed
def bp_test(v):


    import oct2py
    session = oct2py.Oct2Py()

    op = {'matlab_code_path':'/home/rcf-47/mxu/ln/tomominer/tomominer/simulation', 'matlab_tom_2008_code_path':'/home/rcf-47/mxu/proj/imaging/electron/matlab_tom/2008/TOM_Release_2008--octave', 'matlab_tom_2005_code_path':'/home/rcf-47/mxu/proj/imaging/electron/util/tom_2005', 'model':{'missing_wedge_angle':30, 'titlt_angle_step':1, 'SNR':0.1,'SNR bak':N.nan}, 'ctf':{'pix_size':1.0, 'Dz':-15.0, 'voltage':300, 'Cs':2.2, 'sigma':0.4}, 'random_seed':1234}

    vb = do_reconstruction(session, v, op, verbose=False)
    return vb



if __name__ == '__main__':

    import tomominer.model.util as MU
    vb = bp_test(MU.generate_toy_model(dim_siz=128))
    IF.put_mrc(vb, '/tmp/v.mrc')

