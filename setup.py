import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

#os.environ["CC"] = "/Users/Kevin/bin/gcc/bin/gcc"
#os.environ["CXX"] = "/Users/Kevin/bin/gcc/bin/g++"

eca = ['-arch', 'x86_64']
ela = eca

setup(
	  cmdclass = {'build_ext':build_ext},
	  include_dirs = [np.get_include()],
	  ext_modules = [Extension("fmesh",["fmesh.pyx"],
							   libraries=["m"],
							   extra_objects = ["fmesh_wrapper.o","mesh.o"],
							   extra_compile_args=eca,
							   extra_link_args=ela,
							   language="c++")]
	  )
