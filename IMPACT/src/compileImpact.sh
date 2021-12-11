# Compile IMPACT
cd /project/c12/pcorso/hpc-predict/IMPACT/src \
    && HPC_PREDICT_IO_DIR=/project/c12/pcorso/hpc-predict/hpc-predict-io/install make

# Compile the Fortran part of hpc-predict-io (with IMPACT tests)
#RUN cd /src/hpc-predict/hpc-predict-io \
#    && rm -rf build install \
#    && mkdir build install \
#    && cd build \
#    && IMPACT_DIR=/src/hpc-predict/IMPACT/prog/ \
#      cmake -DCMAKE_Fortran_COMPILER=mpifort \
#      -DCMAKE_C_COMPILER=mpicc \
#      -DCMAKE_CXX_COMPILER=mpicxx \
#      -DHDF5_ROOT=/usr/ \
#      -DCMAKE_INSTALL_PREFIX=../install \
#      ../ \
#    && make VERBOSE=1 install
