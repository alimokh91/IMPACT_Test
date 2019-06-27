program mr_io_test_parallel_reader_writer_space_time
  
    use mr_io_parallel_spacetime

    !include 'mpif.h'

    implicit none

    INTEGER err
    character(len=40) :: in_path = "mr_io_test_parallel_space_time_in.h5"
    character(len=40) :: out_path = "mr_io_test_parallel_space_time_out.h5"
    type(DistSpacetimeMRI) :: mri_dest

    integer, dimension(5) :: voxel_feature_shape

    call MPI_Init(err)
    
    call mr_io_read_parallel_spacetime(MPI_COMM_WORLD, MPI_INFO_NULL, in_path, mri_dest)

    print *, SpatialMRI_group_name

    voxel_feature_shape = shape(mri_dest%voxel_feature%array)

    print *, mri_dest%voxel_feature%dims
    print *, mri_dest%voxel_feature%offset
    print *, voxel_feature_shape(3:5)

    print *, mri_dest%voxel_feature%time_dim
    print *, mri_dest%voxel_feature%time_offset
    print *, voxel_feature_shape(2)

    print *, voxel_feature_shape(1)
    print *, mri_dest%voxel_feature%array

    call mr_io_write_parallel_spacetime(MPI_COMM_WORLD, MPI_INFO_NULL, out_path, mri_dest)

    call MPI_Finalize(err)
        
end program mr_io_test_parallel_reader_writer_space_time
