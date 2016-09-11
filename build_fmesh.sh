gfortran -c mesh.f90 fmesh_wrapper.f90;
python setup.py build_ext --inplace;
