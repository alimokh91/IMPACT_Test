#!/bin/bash

export HDF5_DIR=/home/lukasd/src/build_whale/opt_scientific/libs/petsc
gfortran -g -c mr_protocol.f90 -I${HDF5_DIR}/include
gfortran -g -c mr_io.f90 -I${HDF5_DIR}/include
gfortran -g -c *.f90 -I${HDF5_DIR}/include
gfortran mr_io_test.o mr_io.o mr_protocol.o -o mr_io_test -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
gfortran mr_io_test_space_time.o mr_io.o mr_protocol.o -o mr_io_test_space_time -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran

gfortran mr_io_test_reader.o mr_io.o mr_protocol.o -o mr_io_test_reader -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
gfortran mr_io_test_reader_space_time.o mr_io.o mr_protocol.o -o mr_io_test_reader_space_time -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
