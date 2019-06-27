#!/bin/bash

set -euxo pipefail

export HDF5_DIR=/home/lukasd/src/build_whale/opt_scientific/libs/petsc
export MPI_DIR=/usr/include/mpich
mpifort -g -c -fbacktrace mr_protocol.f90 -I${HDF5_DIR}/include
mpifort -g -c -cpp -fbacktrace mr_io.f90  -I${MPI_DIR} -I${HDF5_DIR}/include
mpifort -g -c -cpp -fbacktrace mr_io_parallel.f90  -I${MPI_DIR} -I${HDF5_DIR}/include
mpifort -g -c -cpp -fbacktrace  *.f90 -I${MPI_DIR}  -I${HDF5_DIR}/include
mpifort mr_io_example_spatial.o mr_io.o mr_io_parallel.o mr_protocol.o -o mr_io_example_spatial -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_example_space_time.o mr_io.o mr_io_parallel.o  mr_protocol.o -o mr_io_example_space_time -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran

mpifort mr_io_test_reader.o mr_io.o mr_io_parallel.o  mr_protocol.o -o mr_io_test_reader -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_reader_space_time.o mr_io.o mr_io_parallel.o  mr_protocol.o -o mr_io_test_reader_space_time -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran

mpifort mr_io_test_parallel_reader.o mr_io.o mr_io_parallel.o  mr_protocol.o -o mr_io_test_parallel_reader -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran

mpifort mr_io_test_reader_writer.o mr_io.o mr_io_parallel.o  mr_protocol.o -o mr_io_test_reader_writer -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_reader_writer_space_time.o mr_io.o mr_io_parallel.o  mr_protocol.o -o mr_io_test_reader_writer_space_time -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran

mpifort mr_io_test_parallel_reader_writer.o mr_io.o mr_io_parallel.o  mr_protocol.o -o mr_io_test_parallel_reader_writer -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
