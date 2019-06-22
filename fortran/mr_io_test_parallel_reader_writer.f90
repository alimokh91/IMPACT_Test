program mr_io_test_parallel_reader_writer
  
    use mr_io
    use mr_protocol

    !include 'mpif.h'

    implicit none

    INTEGER err
    character(len=30) :: in_path = "mr_io_test_parallel_in.h5"
    character(len=30) :: out_path = "mr_io_test_parallel_out.h5"
    type(DistSpatialMRI) :: mri_dest    

    call MPI_Init(err)
    
    call mr_io_read_parallel_spatial(MPI_COMM_WORLD, MPI_INFO_NULL, in_path, mri_dest)
    print *, "Array properties:"
    print *, mri_dest%voxel_feature%dims
    print *, mri_dest%voxel_feature%offset
    print *, shape(mri_dest%voxel_feature%array)
    call mr_io_write_parallel_spatial(MPI_COMM_WORLD, MPI_INFO_NULL, out_path, mri_dest)

    call MPI_Finalize(err)
        
end program mr_io_test_parallel_reader_writer
