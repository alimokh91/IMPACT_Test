program mr_io_test_reader
  
    use mr_io
    use mr_protocol

    !include 'mpif.h'

    implicit none

    INTEGER err
    character(len=23) :: path = "mr_io_test_parallel.h5"
    type(SpatialMRI) :: mri_dest    

    call MPI_Init(err)
    
    call mr_io_read_parallel_spatial(MPI_COMM_WORLD, MPI_INFO_NULL, path, mri_dest)

    print *, SpatialMRI_group_name
    print *, mri_dest%voxel_feature_dims    
    print *, mri_dest%voxel_feature

    call MPI_Finalize(err)
        
end program mr_io_test_reader
