from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import sys,os,time
from numpy import argmax

L = os.listdir('.')
L = [f for f in L if f.endswith('.pyx')]
L.sort()
T = [os.path.getmtime(f)  for f in L]
##print L
##print T

last_modified = L[argmax(T)]
last_modified_no_ext = os.path.splitext(last_modified)[0]
 
print last_modified

ext_modules = [Extension(last_modified_no_ext, [last_modified])]

import numpy
setup(
  name = 'my app',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules,
  include_dirs = [numpy.get_include(),],
)

