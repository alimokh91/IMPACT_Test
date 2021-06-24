#!/bin/bash

set -euxo pipefail

export IMPACT_DIR=/home/lukasd/src/IMPACT/prog
export HDF5_DIR=/home/lukasd/src/build_whale/opt_scientific/libs/petsc
export MPI_DIR=/usr/include/mpich

mpifort -g -c -fbacktrace mr_io_protocol.f90 -I${HDF5_DIR}/include
mpicxx  -g -c -cpp mr_io_locking_utils.cxx  -I${MPI_DIR} -I${HDF5_DIR}/include
mpicc   -g -c -cpp mr_io_locking_utils.c  -I${MPI_DIR} -I${HDF5_DIR}/include
mpifort -g -c -cpp -fbacktrace mr_io_locking_utils.f90  -I${MPI_DIR} -I${HDF5_DIR}/include
mpifort -g -c -cpp -fbacktrace mr_io.f90  -I${MPI_DIR} -I${HDF5_DIR}/include
mpifort -g -c -cpp -fbacktrace mr_io_parallel.f90  -I${MPI_DIR} -I${HDF5_DIR}/include
mpifort -g -c -cpp -fbacktrace mr_io_parallel_spacetime.f90  -I${MPI_DIR} -I${HDF5_DIR}/include
mpifort -g -c -cpp -fbacktrace  *.f90 -I${IMPACT_DIR}  -I${MPI_DIR}  -I${HDF5_DIR}/include

cd test
./build.sh
cd ..
