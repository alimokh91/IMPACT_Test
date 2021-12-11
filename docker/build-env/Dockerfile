FROM nvidia/cuda:10.0-cudnn7-devel-ubuntu18.04

ARG MPICH_VERSION=3.1.4

ENV MPICH_VERSION ${MPICH_VERSION}

# Update all packages and install the development toolset (newer version of GNU gcc and tools)
# Dependencies for encfs: libfuse, libssl, gettext
RUN apt-get update  && apt-get install --yes --no-install-recommends \
       wget git curl less tree vim unzip \
       gfortran libblas-dev liblapack-dev zlib1g-dev \
       python-dev python3-dev strace python3-venv \
       libfuse-dev libssl-dev gettext fuse encfs \
    && rm -rf /var/lib/apt/lists/*

# Install MPI
RUN cd /tmp \
    && echo "Using $(which gfortran)" \
    && gfortran --version \
    && wget -q http://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz \
    && tar xf mpich-${MPICH_VERSION}.tar.gz \
    && cd mpich-${MPICH_VERSION} \
    && ./configure --enable-g=debug --enable-fast=O0 --prefix=/usr \
    && make -j3 \
    && make install \
    && ldconfig \
    && cd .. \
    && rm -rf mpich-${MPICH_VERSION} \
    && rm -f mpich-${MPICH_VERSION}.tar.gz

# Install HDF5
RUN cd /tmp \
    && wget -q https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_10_6/source/hdf5-1.10.6.tar.gz \
    && tar xf hdf5-1.10.6.tar.gz \
    && cd hdf5-1.10.6 \
    && FC=mpifort CC=mpicc ./configure --enable-build-mode=debug --enable-symbols=yes --enable-profiling=yes --enable-parallel --enable-fortran --prefix=/usr \
    && make -j3 \
    && make install \
    && ldconfig \
    && cd .. \
    && rm -rf hdf5-1.10.6 \
    && rm -f hdf5-1.10.6.tar.gz

# Install CMake
RUN cd /tmp \
    && wget -q https://cmake.org/files/v3.12/cmake-3.12.3.tar.gz \
    && tar -xf cmake-3.12.3.tar.gz \
    && cd cmake-3.12.3 \
    && ./configure --prefix=/usr \
    && make -j3 \
    && make install \
    && cd .. \
    && rm -rf cmake-3.12.3* 

# Install bundle tooling for creating small Docker images
RUN cd root && \
    wget -q https://github.com/haampie/libtree/releases/download/v1.0.3/libtree_x86_64.tar.gz && \
    tar -xzf libtree_x86_64.tar.gz && \
    rm libtree_x86_64.tar.gz
