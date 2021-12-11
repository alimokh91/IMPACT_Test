#ifdef __GFORTRAN__
! Spatial hyperslab macro enum
! Compute spatial hyperslab
#define mr_io_parallel_spatial_hyperslab_enum_compute   ( 0)
! Read specified hyperslab
#define mr_io_parallel_spatial_hyperslab_enum_specified ( 1)
! Inconsistent arguments (shape/offset/feature_array)
#define mr_io_parallel_spatial_hyperslab_enum_error     (-1)
#endif
