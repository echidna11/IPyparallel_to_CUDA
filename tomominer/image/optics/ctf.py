

import numpy as N


'''

2D and 3D ctf function, written according to 


~/ln/electron/matlab_tom/2005/Analysis/tom_create_ctf.m

~/ln/tomominer/tomominer/image/optics/ctf.py

'''


def create(Dz, size, pix_size=0.72, voltage=300.0, Cs=2.0, sigma=None, display_info=False):

    Cs *= (1e-3)
    voltage *= 1000.0
    pix_size *= (1e-9)

    Dz *= (1e-6)
    Dzn = Dz * 1000000.0          # for display
    Csn =  Cs * 1000.0            # for display
    voltagen = voltage/1000.0       # for display
    voltagest = voltage * (1.0 + voltage / 1022000.0)   # for relativistic calc
    lambda_t = N.sqrt(150.4 / voltagest) * (1e-10)


    Ny = 1.0 / (2.0 * pix_size)
    nyqvist = 2.0 * pix_size * (1e9)

    if display_info:        print 'CTF is calculated for: Defocus', Dzn, '\mum Voltage = ', voltagen, 'kV, Nyqvist = ', nyqvist, 'nm'


    if len(size) == 2:
        g = N.mgrid[0:size[0], 0:size[1]]
    elif len(size) == 3:
        g = N.mgrid[0:size[0], 0:size[1], 0:size[2]]
    else:
        raise Exception('2D or 3D array only')

    g = g.astype(N.float)

    for dim_i in range(len(g)):
        g[dim_i] *= Ny / (size[dim_i] / 2.0)
        g[dim_i] -= Ny

    r = N.sqrt((g ** 2).sum(axis=0))

    vol = N.sin((N.pi / 2.0) * ( Cs* (lambda_t**3.0) * (r**4.0) - 2.0*Dz*lambda_t*(r**2.0) ))
    amplitude = N.cos( (N.pi / 2.0) * (Cs * (lambda_t**3.0) * (r**4.0) - 2.0*Dz*lambda_t*(r**2.0)))

    if sigma is not None:
        # gaussian smoothing....
        vol *=  N.exp(-(r/(sigma*Ny))**2.0)
        amplitude *= N.exp(-(r/(sigma*Ny))**2.0)

    assert N.all(N.isfinite(vol))
    
    return {'ctf':vol, 'amplitude':amplitude, 'g':g}



'''

# code for verification with matlab version

import numpy as N
import numpy.random as NR
import tomominer.image.optics.ctf as IOC

# parameters for example, see ~/ln/frequent_structure/data/out/method/imputation-classification/tomogram-subtomogram/0001/back_projection_reconstruction__config.json

Dz = NR.normal(-15.0, 2)
pix_size = NR.uniform(0.1, 2)
voltage = NR.uniform(250, 350)
Cs = NR.uniform(1.2, 3.2)

#sigma = NR.uniform(0.2, 1.5)
sigma = None

#size = NR.randint(30, 60, size=3)
size = NR.randint(30, 60) * N.ones(3)

c = IOC.create(Dz=Dz, size=size, pix_size=pix_size, voltage=voltage, Cs=Cs, sigma=sigma)


tom_toolbox_2005_path = '/home/rcf-47/mxu/proj/imaging/electron/util/tom_2005'
import pymatlab
session = pymatlab.session_factory(options='-nodesktop -nodisplay')
session.run( 'addpath(\'%s\')'%(tom_toolbox_2005_path,) )


session.run('clear all')
session.run('Dz = %f'%(Dz,))
session.run('pix_size = %f'%(pix_size,))
session.run('voltage = %f'%(voltage,))
session.run('Cs = %f'%(Cs,))

if sigma is not None:   session.run('sigma = %f'%(sigma))

session.putvalue('size', size)

session.run('vol = zeros(size(1), size(2), size(3))')


if sigma is not None:               
    session.run('[ctf, amplitude] = tom_create_ctf(Dz, vol, pix_size, voltage, Cs, sigma)')

else:
    session.run('[ctf, amplitude] = tom_create_ctf(Dz, vol, pix_size, voltage, Cs)')

ctf_t = session.getvalue('ctf')
ctf_t = N.swapaxes(ctf_t, 0, 1)         # this is the fault of matlab's mashgrid()
amplitude_t = session.getvalue('amplitude')
amplitude_t = N.swapaxes(amplitude_t, 0, 1)     # this is the fault of matlab's mashgrid()

(c['ctf'] - ctf_t).max() / ctf_t.max()
rel_c = N.abs((c['ctf'] - ctf_t) / ctf_t)         # see how much relative differences are there
rel_c.max()
i = N.argmax(rel_c)
ctf_t.flatten()[i] / ctf_t.max()                # see where the max difference happens

(c['amplitude'] - amplitude_t).max() / ctf_t.max()
rel_a = N.abs((c['amplitude'] - amplitude_t) / amplitude_t)
rel_a.max()
i = N.argmax(rel_a)
amplitude_t.flatten()[i] / amplitude_t.max()        # see where the max difference happens




# the above tests show some small differences between the two versions

'''





