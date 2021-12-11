#!/bin/sh
if [[ "$-" == *i* ]]; then
        bind '"\e[A":history-search-backward'
        bind '"\e[B":history-search-forward'
fi
# load necessary modules
if [[ ${HOST} == daint* ]]
then
  module load daint-gpu
  if [ $PE_ENV = "CRAY" ] ; then
    module swap PrgEnv-cray PrgEnv-gnu
  fi
fi

export PROJECT_DIR=$(pwd)

#module load vim
#module purge
module load cray-hdf5-parallel
module load CMake
module load cray-libsci
#module load cray-libsci_acc
#module load craype-accel-nvidia60
#module load cray-tpsl
#module load cray-petsc
module load Boost
module load GSL
#module load mpiP

# export compiler wrapper paths
export CC=cc
export CXX=CC
export FC=ftn
export F90=ftn
export F77=ftn

# make dynamic linking the default
export CRAYPE_LINK_TYPE=dynamic
export CRAY_ADD_RPATH=yes

export HPC_PREDICT_IO_DIR=$PROJECT_DIR/MRI_in_IMPACT/hpc-predict-io/install

# aliases
alias mpicc="cc" 
alias mpicxx="CC" 
alias mpifort="ftn" 
alias mpif90="ftn" 
alias mpif77="ftn"
