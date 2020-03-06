import setuptools
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
# Fortran extension
#from numpy.distutils.core import setup, Extension

setup(
    name = 'pymvg',
    version = '0.0.1',
    author="Example Author",
    author_email="esepulveda@protonmail.com",
    description="Multivariate Gaussianisation package",
    long_description='''''',
    long_description_content_type="text/markdown",
    packages=['pymvg'],
    ext_modules = [
        Extension(
            name='ppmt_interface',
            sources= [
                'src/ppmt_fortran/ppmt.pyf',
                'src/ppmt_fortran/acorni.f90',
                'src/ppmt_fortran/backtrmod.f90',
                'src/ppmt_fortran/eig.f90',
                'src/ppmt_fortran/normaldist.f90',
                'src/ppmt_fortran/quicksort.f90',
                'src/ppmt_fortran/solve.f90',
                'src/ppmt_fortran/sortem.f90',
                'src/ppmt_fortran/ppmt_module.f90',
                'src/ppmt_fortran/ppmt_ext.f90',
                ],
            libraries = ['lapack'],
        ),
    ],
)
