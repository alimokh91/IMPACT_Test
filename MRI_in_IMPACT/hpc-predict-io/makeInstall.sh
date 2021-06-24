#!/bin/bash

rm -r build install \
    && mkdir build install \
    && cd build \
    && cmake -DHDF5_ROOT=$HDF5_DIR \
      -DCMAKE_INSTALL_PREFIX=../install \
      -DIMPACT_DIR=/project/c12/pcorso/hpc-predict/IMPACT/prog/ \
      ../ \
    && make VERBOSE=1 install 

    #&& cmake 
    #-DCMAKE_Fortran_COMPILER=mpifort \
    #  -DCMAKE_C_COMPILER=mpicc \
    #  -DCMAKE_CXXCOMPILER=mpicxx \

export HPC_PREDICT_IO_DIR=/project/c12/pcorso/hpc-predict/hpc-predict-io/install/
