# had to install armadillo specific version
# blas and lapack was installed on cluster, I just used libraries (mentioned in include_dirs, lib_dirs list below)
# for gccm default version on UCLA cluster was 4.4.x, but for -std=c++11 to work we need 4.7.x version
# So, first I checked available modules using "module av gcc"
# Then loaded required version using "module load gcc/4.7.x"

import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from numpy import get_include

def get_core_sources(root_dir='tomominer/core'):
    src_list = []
    for root, dirs, files in os.walk(root_dir):
        for ft in files:
            ft_n, ft_e = os.path.splitext(ft)
            ft_e = ft_e.lower()
            if ft_e == '.cpp':   src_list.append(os.path.join(root, ft))
            if ft_e == '.pyx':   src_list.append(os.path.join(root, ft))
    return src_list

cpp_core = Extension('tomominer.core.core',
                    #sources = get_core_sources()
                    sources = [ 'tomominer/core/cython/wrap_core.cpp',
                                'tomominer/core/cython/core.pyx',
                                'tomominer/core/src/affine_transform.cpp',
                                'tomominer/core/src/align.cpp',
                                'tomominer/core/src/arma_extend.cpp',
                                'tomominer/core/src/dilate.cpp',
                                'tomominer/core/src/fft.cpp',
                                'tomominer/core/src/geometry.cpp',
                                'tomominer/core/src/interpolation.cpp',
                                'tomominer/core/src/io.cpp',
                                'tomominer/core/src/legendre.cpp',
                                'tomominer/core/src/rotate.cpp',
                                'tomominer/core/src/sht.cpp',
                                'tomominer/core/src/wigner.cpp',
                                'tomominer/core/src/segmentation/watershed/watershed_segmentation.cpp'
                              ],
                    libraries           = ['m', 'fftw3', 'armadillo', ],
                    # include_dirs        = [get_include(), '/usr/include', 'tomominer/core/src', '/u/local/apps/fftw3/3.3.8-gcc/include/', 'ext_libraries/armadillo/armadillo-destination-dir/usr/include/'],
                    include_dirs        = [get_include(), '/usr/include', 'tomominer/core/src'],
                    # library_dirs        = ['/u/local/apps/fftw3/3.3.8-gcc/lib', 'ext_libraries/armadillo/armadillo-destination-dir/usr/lib64/'],
                    extra_compile_args  = ['-std=c++11'],
                    language='c++',
)


def get_packages(root_dir='tomominer', exclude_dir_roots=['tomominer/core/src', 'tomominer/core/cython']):
    pkg = []
    for root, dirs, files in os.walk(root_dir):
        exclude = False
        for d in exclude_dir_roots:
            if root.startswith(d):  exclude = True
        if exclude:     continue
        pkg.append(root.replace('/', '.'))
    return pkg

setup(  name        = 'tomominer',
        version     = '1.0.0',
        author      = 'Alber Lab (USC)',
        description = 'Subtomogram Analysis and Mining Software',
        license     = 'GPLv3',
        url         = '',
        platforms   = ['x86_64'],
        ext_modules = [cpp_core],

        packages    = get_packages(),
        cmdclass   = {'build_ext': build_ext},
     )
