#!/bin/bash

rm -r build install \
    && mkdir build install \
    && cd build \
    && cmake -DHDF5_ROOT=/opt/cray/pe/hdf5-parallel/1.12.0.0/GNU/8.2 \
      -DCMAKE_INSTALL_PREFIX=../install \
      -DIMPACT_DIR=$PROJECT_DIR/MRI_in_IMPACT/IMPACTwoEnKF/prog \
      ../ \
    && make VERBOSE=1 install
