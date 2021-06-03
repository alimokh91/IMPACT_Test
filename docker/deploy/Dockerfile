# Multistage build: here we import the current source code
# into build environment image, build the project, bundle it
# and then extract it into a small image that just contains
# the binaries we need to run

ARG BUILD_ENV

# FROM $BUILD_ENV as builder

FROM $BUILD_ENV
# FROM hpc-predict:io-build-env

# Copy hpc-predict source code in
COPY . /src/hpc-predict/hpc-predict-io
 
# Compile the Fortran part of hpc-predict-io (without IMPACT tests) and setup python venv
RUN bash -c "cd /src/hpc-predict/hpc-predict-io \
    && mkdir build install \
    && cd build \
    && cmake -DCMAKE_Fortran_COMPILER=mpifort \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_CXX_COMPILER=mpicxx \
      -DHDF5_ROOT=/usr/ \
      -DCMAKE_INSTALL_PREFIX=../install \
      ../ \
    && make VERBOSE=1 install \
    && cd /src/hpc-predict/hpc-predict-io/python \
    && python3 -m venv venv \
    && source venv/bin/activate \
    && pip install --upgrade pip \
    && pip install -r requirements.txt \
    && patch venv/lib/python3.6/site-packages/black/cache.py requirements_black_cache.patch" # https://github.com/psf/black/pull/1224/files

# FIXME: Bundling remains TBD - question is especially how to apply this to python venvs/whether they should be built inside final image only 
#      \
#    && /root/libtree/libtree --chrpath --strip \
#      -d /root/hpc-predict.bundle \
#      /src/hpc-predict       
#      \
#    && rm -rf /src/hpc-predict


# FROM ubuntu:18.04
# 
# # This is the only thing necessary really from nvidia/cuda's ubuntu18.04 runtime image
# ENV NVIDIA_VISIBLE_DEVICES all
# ENV NVIDIA_DRIVER_CAPABILITIES compute,utility
# ENV NVIDIA_REQUIRE_CUDA "cuda>=10.2 brand=tesla,driver>=384,driver<385 brand=tesla,driver>=396,driver<397 brand=tesla,driver>=410,driver<411 brand=tesla,driver>=418,driver<419"
# 
# COPY --from=builder /root/hpc-predict.bundle /root/hpc-predict.bundle
# 
# # Make it easy to call our binaries.
# ENV PATH="/root/hpc-predict.bundle/usr/bin:$PATH"
# 
# RUN echo "/root/hpc-predict.bundle/usr/lib/" > /etc/ld.so.conf.d/hpc-predict.conf && ldconfig
# 
# WORKDIR /root/hpc-predict.bundle/usr/bin
