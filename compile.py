from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import numpy

ext = Extension("*", ["src/*.pyx"],
    include_dirs = [numpy.get_include()])

setup(name='PAHA',
      version='0.1',
      description='Python for Atomic Halo\'s Analysis',
      author='Fernando Becerra',
      author_email='fbecerra@cfa.harvard.edu',
      url='https://www.cfa.harvad.edu/~fbecerra',
      ext_modules=cythonize(ext),
      cmdclass = {'build_ext': build_ext})
