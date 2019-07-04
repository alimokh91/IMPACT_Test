# Fortran-Python MRI data communication layer for HPC-PREDICT

The functionality to communicate MRI objects across Python and Fortran is contained in the `python` and `fortran` directories. Integration apps accessing various data sources are kept in `mri_datasource` directory. 

## Communication library

The objects being communicated are currently a spatial voxel-based scalar field (`SpatialMRI`), a spatiotemporal voxel-based vectorial field as well as coordinates along each of the x-/y-/z-/t-axis (`SpaceTimeMRI`) and a variation of the latter with two spatiotemporal voxel-based vectorial fields for velocity mean and covariance (`HPCPredictMRI`) as specified in `fortran/mr_protocol.f90`

The HDF5 group name is used for the message label, every array is in a separate dataset under this group and the memory layout is according to Fortran (column major).
That is any necessary conversions (spatial index inversion to have forward spatial indexes in each programming language) are done in Python. For vector-valued fields on a space-time MRI we use the index ordering (x,y,z,t,i) in Python and (i,t,x,y,z) in Fortran (i is the component-index of the pointwise vector).

# Building, installing and usage from IMPACT

To build the Fortran library and install it, run 

```
export FC=mpifort # (export FC=ftn for daint)
export HDF5_ROOT=/path/where/hdf5/can/be/found
mkdir build
cd build
FC=mpifort cmake -DCMAKE_INSTALL_PREFIX=/path/where/you/want/to/install/hpc-predict-io ../
make install
```
For usage in IMPACT set the HPC_PREDICT_IO_DIR environment variable to /path/to/hpc-predict-io/install before running make in the src directory.

# Tests

To run Python-to-Fortran and round trip tests for hpc-predict-io execute

```
source python/venv/bin/activate
./test.sh
```

