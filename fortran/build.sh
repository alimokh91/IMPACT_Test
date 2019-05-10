#!/bin/bash

export HDF5_DIR=/home/lukasd/src/build_whale/opt_scientific/libs/petsc
gfortran -g -c mr_protocol.f90 -I${HDF5_DIR}/include
gfortran -g -c mr_io.f90 -I${HDF5_DIR}/include
gfortran -g -c mr_io_test.f90 -I${HDF5_DIR}/include
gfortran *.o -o mr_io_test -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
