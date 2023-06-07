
'''
functions to calculate spherical harmonics features
'''
import numpy as N

# SH expansion of v from center, by calling the matlab routine SphericalHarmonicsUtil.adp_em_sh_expansion()
def sh_expansion_matlab(v, op):
    if 'session' not in op:
        import pymatlab
        op['session'] = pymatlab.session_factory(options='-nodesktop -nodisplay')
        op['session'].run( 'addpath(\'%s\')'%(op['expansion code path'],))
        op['session'].run( 'addpath(\'%s\')'%(op['DOTM code path'],))


    session = op['session']
    session.run('clear all')
    session.run('max_l = %d'%(op['max_l']))
    session.putvalue('v', v.astype(N.float))

    session.run('save /tmp/t')

    session.run('[coef_real, coef_imag, rads] = SphericalHarmonicsUtil.adp_em_expansion(v, max_l, [0, 0.5])')

    rads = session.getvalue('rads')

    coef_real = session.getvalue('coef_real')
    coef_imag = session.getvalue('coef_imag')

    return {'rads':rads, 'coef_real':coef_real, 'coef_imag':coef_imag}


'''
# test code

op = {'max_l':3, 'expansion code path':'/home/rcf-47/mxu/ln/frequent_structure/code', 'DOTM code path': '/home/rcf-47/mxu/ln/electron/util/spherical_harmonics/DOTM/src'}


import tomominer.model.util as MU
v = MU.generate_toy_model()
import tomominer.geometry.rotate as GR
import tomominer.geometry.ang_loc as GA
vr = GR.rotate(v, rm=GA.random_rotation_matrix())

import tomominer.image.feature.spherical_harmonics as IFS

s = IFS.sh_expansion_matlab(vr, op)

f = IFS.rotation_invariant_feature(s['coef_real'])
print f['f']

'''


# c is the SH obtained from the real part of the volume
def rotation_invariant_feature(c):
    rad = N.real(c[:,0].flatten())
    l_s = N.real(c[:,1].flatten()).astype(N.int)
    m_s = N.real(c[:,2].flatten()).astype(N.int)
    c_s = c[:,3].flatten()

    rad_set = list(set(rad))        ;       rad_set = sorted(rad_set)
    l_set = list(set(l_s))      ;       l_set = sorted(l_set)

    rif = []
    rif_r = []
    rif_l = []
    for r in rad_set:
        for l in l_set:
            rif.append(N.linalg.norm(c_s[(rad == r) & (l_s == l)]))
            rif_r.append(r)
            rif_l.append(l)

    return  {'f':rif, 'r':rif_r, 'l':rif_l}

