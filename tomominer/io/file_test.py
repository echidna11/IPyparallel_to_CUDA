import os
import sys
import pickle
import shutil
import time
import struct
import uuid

import array
import numpy as N



def get_mrc(path, retry_interval=1.0, max_retry=5):
    """
    Use filesystem to fetch data.
    """
    path = os.path.realpath(str(path))

    import tomominer.core as tomo

    v = None

    retry = 0
    while retry < max_retry:

        try :
            v = None
            v = tomo.read_mrc(path)
            break

        except:
            retry += 1
            time.sleep(retry_interval)      # if load failed, retry a number of times, this is because sometimes some file system is not stable

    if v is None:    raise IOError('cannot load ' + path)

    return v


'''
translated according to tom_toolbox's mrcread()

WARNING: this function cannot correctly read etomo generated full tomogram mrc file, in order to correctly read it! 
There is one way to work around this peoblem: first use tom toolbox's tom_mrcread to load it then save using tom_mrcwrite as classical style!!!
'''

def read_mrc(path, read_data=True, show_progress=False):

    '''
    %Structure of MRC-data files:
    %MRC Header has a length of 1024 bytes
    % SIZE  DATA    NAME    DESCRIPTION
    %   4   int     NX      number of Columns    (fastest changing in map)
    %   4   int     NY      number of Rows
    %   4   int     NZ      number of Sections   (slowest changing in map)
    %   4   int     MODE    Types of pixel in image
    %                       0 = Image     unsigned bytes
    %                       1 = Image     signed short integer (16 bits)
    %                       2 = Image     float
    %                       3 = Complex   short*2
    %                       4 = Complex   float*2     
    %	4   int     NXSTART Number of first COLUMN  in map (Default = 0)
    %   4   int     NYSTART Number of first ROW     in map      "
    %   4   int     NZSTART Number of first SECTION in map      "
    %   4   int     MX      Number of intervals along X
    %   4   int     MY      Number of intervals along Y
    %   4   int     MZ      Number of intervals along Z
    %   4   float   Xlen    Cell Dimensions (Angstroms)
    %   4   float   Ylen                 "
    %   4   float   Zlen                 "
    %   4   float   ALPHA   Cell Angles (Degrees)
    %   4   float   BETA                 "
    %   4   float   GAMMA                "
    %   4   int     MAPC    Which axis corresponds to Columns  (1,2,3 for X,Y,Z)
    %   4   int     MAPR    Which axis corresponds to Rows     (1,2,3 for X,Y,Z)
    %   4   int     MAPS    Which axis corresponds to Sections (1,2,3 for X,Y,Z)
    %   4   float   AMIN    Minimum density value
    %   4   float   AMAX    Maximum density value
    %   4   float   AMEAN   Mean    density value    (Average)
    %   2   short   ISPG    Space group number       (0 for images)
    %   2   short   NSYMBT  Number of bytes used for storing symmetry operators
    %   4   int     NEXT    Number of bytes in extended header
    %   2   short   CREATID Creator ID
    %   30    -     EXTRA   Not used. All set to zero by default
    %   2   short   NINT    Number of integer per section
    %   2   short   NREAL   Number of reals per section
    %   28    -     EXTRA2  Not used. All set to zero by default
    %   2   short   IDTYPE  0=mono, 1=tilt, 2=tilts, 3=lina, 4=lins
    %   2   short   LENS    
    %   2   short   ND1   
    %   2   short   ND2
    %   2   short   VD1 
    %   2   short   VD2
    %   24  float   TILTANGLES
    %   4   float   XORIGIN X origin
    %   4   float   YORIGIN Y origin
    %   4   float   ZORIGIN Z origin
    %   4   char    CMAP    Contains "MAP "
    %   4   char    STAMP   
    %   4   float   RMS 
    %   4   int     NLABL   Number of labels being used
    %   800 char    10 labels of 80 character
    %
    %Extended Header (FEI format and IMOD format)
    %The extended header contains the information about a maximum of 1024 images. 
    %Each section is 128 bytes long. The extended header is thus 1024 * 128 bytes 
    %(always the same length, regardless of how many images are present
    %   4   float   a_tilt  Alpha tilt (deg)
    %   4   float   b_tilt  Beta tilt (deg)
    %   4   float   x_stage  Stage x position (Unit=m. But if value>1, unit=???m)
    %   4   float   y_stage  Stage y position (Unit=m. But if value>1, unit=???m)
    %   4   float   z_stage  Stage z position (Unit=m. But if value>1, unit=???m)
    %   4   float   x_shift  Image shift x (Unit=m. But if value>1, unit=???m)
    %   4   float   y_shift  Image shift y (Unit=m. But if value>1, unit=???m)
    %   4   float   z_shift  Image shift z (Unit=m. But if value>1, unit=???m)
    %   4   float   defocus  Defocus Unit=m. But if value>1, unit=???m)
    %   4   float   exp_time Exposure time (s)
    %   4   float   mean_int Mean value of image
    %   4   float   tilt_axis   Tilt axis (deg)
    %   4   float   pixel_size  Pixel size of image (m)
    %   4   float   magnification   Magnification used
    %   4   float   remainder   Not used (filling up to 128 bytes)   
    %


    for python unpacking of binary data, see http://docs.python.org/2/library/struct.html#struct.calcsize

    '''

    path = os.path.realpath(path)

    with open(path, 'rb') as f:

        mrc = {}
        mrc['nx'] = int(struct.unpack('i', f.read(4))[0])       # MRC.nx = fread(fid,[1],'int');        %integer: 4 bytes
        mrc['ny'] = int(struct.unpack('i', f.read(4))[0])       # MRC.ny = fread(fid,[1],'int');        %integer: 4 bytes
        mrc['nz'] = int(struct.unpack('i', f.read(4))[0])       # MRC.nz = fread(fid,[1],'int');        %integer: 4 bytes
        mrc['mode'] = struct.unpack('i', f.read(4))[0]       # MRC.mode = fread(fid,[1],'int');      %integer: 4 bytes
        mrc['nxstart'] = struct.unpack('i', f.read(4))[0]       # MRC.nxstart= fread(fid,[1],'int');    %integer: 4 bytes
        mrc['nystart'] = struct.unpack('i', f.read(4))[0]       # MRC.nystart= fread(fid,[1],'int');    %integer: 4 bytes
        mrc['nzstart'] = struct.unpack('i', f.read(4))[0]       # MRC.nzstart= fread(fid,[1],'int');    %integer: 4 bytes
        mrc['mx'] = struct.unpack('i', f.read(4))[0]       # MRC.mx= fread(fid,[1],'int');         %integer: 4 bytes
        mrc['my'] = struct.unpack('i', f.read(4))[0]       # MRC.my= fread(fid,[1],'int');         %integer: 4 bytes
        mrc['mz'] = struct.unpack('i', f.read(4))[0]       # MRC.mz= fread(fid,[1],'int');         %integer: 4 bytes
        mrc['xlen'] = struct.unpack('f', f.read(4))[0]       # MRC.xlen= fread(fid,[1],'float');     %float: 4 bytes
        mrc['ylen'] = struct.unpack('f', f.read(4))[0]       # MRC.ylen= fread(fid,[1],'float');     %float: 4 bytes
        mrc['zlen'] = struct.unpack('f', f.read(4))[0]       # MRC.zlen= fread(fid,[1],'float');     %float: 4 bytes
        mrc['alpha'] = struct.unpack('f', f.read(4))[0]       # MRC.alpha= fread(fid,[1],'float');    %float: 4 bytes
        mrc['beta'] = struct.unpack('f', f.read(4))[0]       # MRC.beta= fread(fid,[1],'float');     %float: 4 bytes
        mrc['gamma'] = struct.unpack('f', f.read(4))[0]       # MRC.gamma= fread(fid,[1],'float');    %float: 4 byte
        mrc['mapc'] = struct.unpack('i', f.read(4))[0]       # MRC.mapc= fread(fid,[1],'long');       %integer: 4 bytes
        mrc['mapr'] = struct.unpack('i', f.read(4))[0]       # MRC.mapr= fread(fid,[1],'long');       %integer: 4 bytes
        mrc['maps'] = struct.unpack('i', f.read(4))[0]       # MRC.maps= fread(fid,[1],'long');       %integer: 4 bytes
        mrc['amin'] = struct.unpack('f', f.read(4))[0]       # MRC.amin= fread(fid,[1],'float');     %float: 4 bytes
        mrc['amax'] = struct.unpack('f', f.read(4))[0]       # MRC.amax= fread(fid,[1],'float');     %float: 4 bytes
        mrc['amean'] = struct.unpack('f', f.read(4))[0]       # MRC.amean= fread(fid,[1],'float');    %float: 4 bytes
        mrc['ispg'] = struct.unpack('h', f.read(2))[0]       # MRC.ispg= fread(fid,[1],'short');     %integer: 2 bytes
        mrc['nsymbt'] = struct.unpack('h', f.read(2))[0]       # MRC.nsymbt = fread(fid,[1],'short');  %integer: 2 bytes
        mrc['next'] = struct.unpack('i', f.read(4))[0]       # MRC.next = fread(fid,[1],'int');      %integer: 4 bytes
        mrc['creatid'] = struct.unpack('h', f.read(2))[0]       # MRC.creatid = fread(fid,[1],'short'); %integer: 2 bytes
        mrc['unused1'] = struct.unpack('c'*30, f.read(30))[0]       # MRC.unused1 = fread(fid,[30]);        %not used: 30 bytes
        mrc['nint'] = struct.unpack('h', f.read(2))[0]       # MRC.nint = fread(fid,[1],'short');    %integer: 2 bytes
        mrc['nreal'] = struct.unpack('h', f.read(2))[0]       # MRC.nreal = fread(fid,[1],'short');   %integer: 2 bytes
        mrc['unused2'] = struct.unpack('c'*28, f.read(28))[0]       # MRC.unused2 = fread(fid,[28]);        %not used: 28 bytes
        mrc['idtype'] = struct.unpack('h', f.read(2))[0]       # MRC.idtype= fread(fid,[1],'short');   %integer: 2 bytes
        mrc['lens'] = struct.unpack('h', f.read(2))[0]       # MRC.lens=fread(fid,[1],'short');      %integer: 2 bytes
        mrc['nd1'] = struct.unpack('h', f.read(2))[0]       # MRC.nd1=fread(fid,[1],'short');       %integer: 2 bytes
        mrc['nd2'] = struct.unpack('h', f.read(2))[0]       # MRC.nd2 = fread(fid,[1],'short');     %integer: 2 bytes
        mrc['vd1'] = struct.unpack('h', f.read(2))[0]       # MRC.vd1 = fread(fid,[1],'short');     %integer: 2 bytes
        mrc['vd2'] = struct.unpack('h', f.read(2))[0]       # MRC.vd2 = fread(fid,[1],'short');     %integer: 2 bytes

        mrc['tiltangles'] = struct.unpack('f'*6, f.read(4*6))      # for i=1:6;  MRC.tiltangles(i)=fread(fid,[1],'float');%float: 4 bytes,   24 bytes in total
            
        mrc['xorg'] = struct.unpack('f', f.read(4))[0]       # MRC.xorg = fread(fid,[1],'float');    %float: 4 bytes
        mrc['yorg'] = struct.unpack('f', f.read(4))[0]       # MRC.yorg = fread(fid,[1],'float');    %float: 4 bytes
        mrc['zorg'] = struct.unpack('f', f.read(4))[0]       # MRC.zorg = fread(fid,[1],'float');    %float: 4 bytes
        mrc['cmap'] = struct.unpack('c'*4, f.read(4))       # MRC.cmap = fread(fid,[4],'char');     %Character: 4 bytes
        mrc['stamp'] = struct.unpack('c'*4, f.read(4))       # MRC.stamp = fread(fid,[4],'char');    %Character: 4 bytes
        mrc['rms'] = struct.unpack('f', f.read(4))[0]       # MRC.rms=fread(fid,[1],'float');       %float: 4 bytes
        mrc['nlabl'] = struct.unpack('i', f.read(4))[0]       # MRC.nlabl = fread(fid,[1],'int');     %integer: 4 bytes
        mrc['labl'] = struct.unpack('c'*800, f.read(800))       # MRC.labl = fread(fid,[800],'char');   %Character: 800 bytes

        size = [mrc['nx'], mrc['ny'], mrc['nz']]
        n_voxel = N.prod(size)

        extended = {}
        extended['magnification'] = [0]
        extended['exp_time'] = [0]
        extended['pixelsize'] = [0]
        extended['defocus'] = [0]
        extended['a_tilt'] = [0] * mrc['nz']
        extended['tiltaxis'] = [0]

        if mrc['next'] != 0:        # extended header
            nbh = mrc['next'] / 128     # 128 = length of FEI extended header
            if nbh == 1024:         # FEI extended header
                for lauf in range(nbh):
                    extended['a_tilt'][lauf] = struct.unpack('f', f.read(4))[0]    # Extended.a_tilt(lauf)= fread(fid,[1],'float');        %float: 4 bytes
                    extended['b_tilt'][lauf] = struct.unpack('f', f.read(4))[0]    # Extended.b_tilt(lauf)= fread(fid,[1],'float');        %float: 4 bytes
                    extended['x_stage'][lauf] = struct.unpack('f', f.read(4))[0]    # Extended.x_stage(lauf)= fread(fid,[1],'float');       %float: 4 bytes
                    extended['y_stage'][lauf] = struct.unpack('f', f.read(4))[0]    # Extended.y_stage(lauf)=fread(fid,[1],'float');        %float: 4 bytes
                    extended['z_stage'][lauf] = struct.unpack('f', f.read(4))[0]    # Extended.z_stage(lauf)=fread(fid,[1],'float');        %float: 4 bytes
                    extended['x_shift'][lauf] = struct.unpack('f', f.read(4))[0]    # Extended.x_shift(lauf)=fread(fid,[1],'float');        %float: 4 bytes
                    extended['y_shift'][lauf] = struct.unpack('f', f.read(4))[0]    # Extended.y_shift(lauf)=fread(fid,[1],'float');        %float: 4 bytes
                    extended['defocus'][lauf] = struct.unpack('f', f.read(4))[0]    # Extended.defocus(lauf)=fread(fid,[1],'float');        %float: 4 bytes
                    extended['exp_time'][lauf] = struct.unpack('f', f.read(4))[0]    # Extended.exp_time(lauf)=fread(fid,[1],'float');       %float: 4 bytes
                    extended['mean_int'][lauf] = struct.unpack('f', f.read(4))[0]    # Extended.mean_int(lauf)=fread(fid,[1],'float');       %float: 4 bytes
                    extended['tiltaxis'][lauf] = struct.unpack('f', f.read(4))[0]    # Extended.tiltaxis(lauf)=fread(fid,[1],'float');       %float: 4 bytes
                    extended['tiltaxis'][lauf] = struct.unpack('f', f.read(4))[0]    # Extended.tiltaxis(lauf)=fread(fid,[1],'float');       %float: 4 bytes
                    extended['pixelsize'][lauf] = struct.unpack('f', f.read(4))[0]    # Extended.pixelsize(lauf)=fread(fid,[1],'float');      %float: 4 bytes
                    extended['magnification'][lauf] = struct.unpack('f', f.read(4))[0]    # Extended.magnification(lauf)=fread(fid,[1],'float');  %float: 4 bytes
                    f.seek(offset=128-52, whence=1)         # fseek(fid,128-52,0);
                else:
                    # IMOD extended Header
                    f.seek(offset=MRC.next, whence=1)       # fseek(fid,MRC.next,'cof');%go to end end of extended Header

        if read_data:
            slice_voxel_num = mrc['nx'] * mrc['ny']
            v = None
            for i in range(mrc['nz']):
                if show_progress:
                    print ('\r', i, '   '),
                    sys.stdout.flush()

                # Note, when reading large array, numpy.fromfile should be much faster than struct.unpack
                if mrc['mode'] == 0:
                    if v is None:       v = N.zeros(size, dtype=N.int8)            # important, in the new specification, format could be either uint8 or int8, see http://bio3d.colorado.edu/imod/doc/mrc_format.txt
                    #data_read = struct.unpack('b'*slice_voxel_num, f.read(slice_voxel_num))         # Data_read(:,:,i) = fread(fid,[MRC.nx,MRC.ny],'int8');
                    data_read = N.fromfile(f, dtype=N.int8, count=slice_voxel_num)                  # Data_read(:,:,i) = fread(fid,[MRC.nx,MRC.ny],'int8');
                elif mrc['mode'] == 1:
                    if v is None:       v = N.zeros(size, dtype=N.int16)
                    #data_read = struct.unpack('<'+'h'*slice_voxel_num, f.read(slice_voxel_num*2))         # Data_read(:,:,i) = fread(fid,[MRC.nx,MRC.ny],'int16');
                    data_read = N.fromfile(f, dtype=N.int16, count=slice_voxel_num)                     # Data_read(:,:,i) = fread(fid,[MRC.nx,MRC.ny],'int16');
                    #data_read = array.array('i');                       data_read.fromfile(f, slice_voxel_num);                     data_read = N.array(data_read, dtype=N.int16)

                elif mrc['mode'] == 2:
                    if v is None:       v = N.zeros(size, dtype=N.float32)
                    #data_read = struct.unpack('f'*slice_voxel_num, f.read(slice_voxel_num*4))         # Data_read(:,:,i) = fread(fid,[MRC.nx,MRC.ny],'float');
                    data_read = N.fromfile(f, dtype=N.float32, count=slice_voxel_num)                   # Data_read(:,:,i) = fread(fid,[MRC.nx,MRC.ny],'float');
                else:
                    raise Exception('Sorry, i cannot read this as an MRC-File !!!')
                    data_read = None

                if data_read.size != slice_voxel_num:
                    import pdb;     pdb.set_trace()

                v[:,:,i] = N.reshape(data_read, (mrc['nx'], mrc['ny']), order='F')



        else:

            v = None


        h = {}
        h['Voltage'] = None
        h['Cs'] = None
        h['Aperture'] = None
        h['Magnification'] = extended['magnification'][0]
        h['Postmagnification'] = None
        h['Exposuretime'] = extended['exp_time'][0]
        h['Objectpixelsize'] =  extended['pixelsize'][0] * 1e9
        h['Microscope'] = None
        h['Pixelsize'] = None
        h['CCDArea'] = None
        h['Defocus'] = extended['defocus'][0]
        h['Astigmatism'] = None
        h['AstigmatismAngle'] = None
        h['FocusIncrement'] = None
        h['CountsPerElectron'] = None
        h['Intensity'] = None
        h['EnergySlitwidth'] = None
        h['EnergyOffset'] = None
        h['Tiltangle'] = extended['a_tilt'][:mrc['nz']]
        h['Tiltaxis'] = extended['tiltaxis'][0]
        h['Username'] = None
        h['Date'] = None
        h['Size'] = [mrc['nx'], mrc['ny'], mrc['nz']]  
        h['Comment'] = None
        h['Parameter'] = None
        h['Fillup'] = None
        h['Filename'] = path
        h['Marker_X'] = None
        h['Marker_Y'] = None
        h['MRC'] = mrc


    return {'header':h, 'value':v}


def read_mrc_vol(path, show_progress=False):
    return read_mrc(path=path, show_progress=show_progress)['value']


def read_mrc_header(path):
    return read_mrc(path=path, read_data=False)['header']



'''


import tomominer.io.file as IF
mrc_file = '/home/rcf-47/mxu/ln/frequent_structure/data/out/analysis/beck/u2os/classification/tomograms/original/1.rec'
i = IF.read_mrc(mrc_file, read_data=False)
print i['header']['MRC']




#mrc_file = '/home/rcf-47/mxu/ln/frequent_structure/data/out/analysis/beck/extract/negative_stain/TS_CT/original/subtomograms/part_1.mrc'
#mrc_file = '/home/rcf-47/mxu/ln/frequent_structure/data/out/analysis/beck/u2os/classification/tomograms/original/1.rec'
mrc_file = '/tmp/1-full.mrc'
import tomominer.io.file as IF
IF.read_mrc__check(mrc_file)



'''


def read_mrc__matlab(mrc_file, session=None, matlab_tom_2008_code_path = '/home/rcf-47/mxu/proj/imaging/electron/matlab_tom/2008/TOM_Release_2008', show_progress=False, slice_chunk=10, sub_region=None):
    mrc_file = os.path.realpath(mrc_file)

    if session is None:
        import pymatlab
        session = pymatlab.session_factory(options='-nodesktop -nodisplay')
        session.run( 'addpath(  genpath(\'%s\') )'%(matlab_tom_2008_code_path,) )

    session.run('clear all;')
    session.run('im = tom_mrcread(\'%s\');'%(mrc_file,))
    session.run('v = im.Value;')
    session.run('siz = int64(size(v));')
    siz = session.getvalue('siz')

    if sub_region is not None:
        sub_region = N.copy(sub_region)
        for dim_i in range(len(siz)):            sub_region[dim_i,1] = min(sub_region[dim_i,1], siz[dim_i])         # constrain size so that the region does not exceed the size
        session.run(    'vt = v(%d:%d, %d:%d, %d:%d);'%(sub_region[0,0]+1, sub_region[0,1], sub_region[1,0]+1, sub_region[1,1], sub_region[2,0]+1, sub_region[2,1])   )
        return session.getvalue('vt')


    if True:
        if False:
            # problem: session.getvalue() cannot handle large arrays
            vt = session.getvalue('v')
        else:
            # read by slice chunks
            vt = None
            slices = list(      range(siz[2])    )
            while slices:
                slices_t = slices[:slice_chunk]

                if show_progress:   print ('\r', slices_t, '            '),            ;       sys.stdout.flush()

                session.putvalue('slices_t', slices_t)
                session.run( 'vs = v(:,:,(slices_t+1))' )       # profiling show that this is the most time consuming step

                vs = session.getvalue('vs')
                   
                if vs.ndim < 3:     vs = N.expand_dims(vs, 2)

                if vt is None:      vt = N.zeros(siz, dtype=vs.dtype)

                vt[:,:, slices_t] = vs

                slices = slices[slice_chunk:]
    else:

        session.run('save(\'/tmp/t.mat\', \'v\', \'-v7.3\');')

        if True:
            import scipy.io as SI
            vt = SI.loadmat('/tmp/t.mat')['v']
        else:
            import h5py
            f = h5py.File('/tmp/t.mat', 'r')
            import pdb;     pdb.set_trace()

    return vt


def read_mrc__octave(mrc_file, matlab_tom_2008_code_path = '/home/rcf-47/mxu/proj/imaging/electron/matlab_tom/2008/TOM_Release_2008'):
    from oct2py import octave
    octave.addpath(matlab_tom_2008_code_path)

    im = octave.tom_mrcread(mrc_file)

    import pdb;     pdb.set_trace()

    return
 

def read_mrc__chimera(mrc_file):
    mrc_file = os.path.realpath(mrc_file)

    #raise       Exception('test.py shows that after transpose, there is still a small displacement, need to correct')

    import uuid
    out_file_root = os.path.join('/tmp', str(uuid.uuid1()))
    mrc2pickle = os.path.join(os.path.dirname(__file__), 'chimera', 'mrc2pickle.sh')

    import subprocess
    subprocess.call([mrc2pickle, mrc_file, out_file_root])

    out_file_vol = out_file_root+'-vol.npy'
    v = N.load(out_file_vol)
    os.remove(out_file_vol)

    out_file_header = out_file_root+'-header.pickle'
    os.remove(out_file_header)

    return v

    '''
    with open(out_file_header, 'rb') as f:      header = pickle.load(f)

    return {'header':header, 'value':v}
    '''


    


def put_mrc(mrc, path, overwrite=True):
    path = os.path.realpath(path)

    if mrc.dtype != N.float:        mrc = N.array(mrc, order='F', dtype=N.float)
    if not mrc.flags['F_CONTIGUOUS']:      mrc = N.array(mrc, order='F', dtype=N.float)

    path = str(path)
    if (overwrite == False) and os.path.isfile(path):        # mxu: to prevent trying to overwrite an existing file owned by another user
        return

    import tomominer.core as core
    core.write_mrc(mrc, path)

'''
translated from tom_mrcwrite(), write classic style only
'''
def write_mrc(mrc, path):
    raise Exception('to be implemented')
    return


'''
when we have a large volume, tomo.write_mrc does not work properly. In such case, we first write a series of tif images, and put these images using tif2mrc
then the tif images are put together using newstack program in imod package
'''
def write_mrc_by_chunk(v, path, n_chunk=5):
   from scipy.misc import imsave
     
   batch_id = 'mrc-chunk--' + str(uuid.uuid4())
   inds = list(range(v.shape[2]))

   file_names = []
   while inds:
       inds_t = inds[:n_chunk]

       f_t = os.path.join('/tmp', '%s--%05d.mrc'%(batch_id,inds_t[0]))
       file_names.append(f_t)
       put_mrc(path=f_t, mrc=v[:,:,inds_t])
       inds = inds[n_chunk:]


   cmd = ['newstack']
   cmd.extend(file_names)
   cmd.append(path)

   import subprocess
   subprocess.call(cmd)

   for f_t in file_names:      os.remove(f_t)

'''
# test code

from tomominer.model.util import generate_toy_model
v = generate_toy_model()

import tomominer.io.file as IF
IF.write_mrc_by_chunk(v, '/tmp/out.mrc', n_chunk=1)

IF.put_mrc(v, '/tmp/out-true.mrc')

'''


def write_mrc__matlab(v, path, session=None, matlab_tom_2008_code_path = '/home/rcf-47/mxu/proj/imaging/electron/matlab_tom/2008/TOM_Release_2008'):
    path = os.path.realpath(path)

    if session is None:
        import pymatlab
        session = pymatlab.session_factory(options='-nodesktop -nodisplay')
        session.run( 'addpath(  genpath(\'%s\') )'%(matlab_tom_2008_code_path,) )

    session.run('clear all;')
    session.putvalue('v', v)
    session.putvalue('fn', str(path))
    session.run('tom_mrcwrite(v, \'name\', fn);')





def np_load(fn):
    fn = os.path.realpath(str(fn))

    with open(fn, 'rb') as f:
        o = N.load(f)

    return o    

# mxu: save a numpy file only when it does not exist
def np_save(fn, arr):
    fn = os.path.realpath(str(fn))

    if not os.path.isfile(fn):        # mxu: to prevent trying to overwrite an existing file owned by another user
        with open(fn, 'wb') as f:
            N.save(f, arr)





