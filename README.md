# Fortran-Python MRI data communication layer for HPC-PREDICT

Reference object(s) being communicated is specified in `fortran/mr_protocol.f90`

An (unfinished) test can be run as follows

```
source python/venv/bin/activate
python python/mr_io_test_writer.py 
gdb fortran/mr_io_test_reader
```

Note that it's segfaulting in the Fortran reader code currently and not reconstructing the matrix of doubles correctly.
