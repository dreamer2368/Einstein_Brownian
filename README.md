# Einstein_Brownian
Adjoint_particle approach for Brownian motion

mesh.f90: base fortran90 code
fmesh_wrapper.f90: fortran90 wrapper that generate c-object
fmesh.pyx: cython code that uses functions from c-object (c-object generated from fmesh_wrapper.f90)
setup.py: a python file that links all the codes above. You need to execute 'python setup.py build_ext --inplace' first before compiling your code

Once you prepared all 4 kinds of file described above, you need to:
1. compile .f90 files
2. compile setup.py file

These are executed by commands
gfortran -c mesh.f90 fmesh_wrapper.f90;
python setup.py build_ext --inplace;

And this command is written in file
build.sh

THUS ONCE YOU PREPARED ALL 4 KINDS OF FILE DESCRIBED ABOVE, JUST TYPE THIS BEFORE YOU ACTUALLY COMPILE IT:
sh build.sh;
