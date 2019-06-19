# Fortran-Python MRI data communication layer for HPC-PREDICT

Objects being communicated is specified in `fortran/mr_protocol.f90`

The HDF5 group name is used for the message label, every array is in a separate dataset under this group and the memory layout is according to Fortran (column major).
That is any necessary conversions (spatial index inversion to have forward spatial indexes in each programming language) are done in Python. For vector-valued fields on a space-time MRI we use the index ordering (x,y,z,t,i) in Python and (i,t,x,y,z) in Fortran (i is the component-index of the pointwise vector).

To run Python-to-Fortran and round trip tests execute

```
source python/venv/bin/activate
./test.sh
```
