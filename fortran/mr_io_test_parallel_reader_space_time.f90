program mr_io_test_parallel_reader_space_time
  
    use mr_io_parallel_spacetime
    use mr_io_test_arg_parser

    !include 'mpif.h'

    implicit none

    INTEGER err
!    character(len=40) :: path = "mr_io_test_parallel_space_time.h5"
    type(DistSpacetimeMRI) :: mri_dest

    integer, dimension(5) :: voxel_feature_shape

    call MPI_Init(err)
    
    call mr_io_test_parse_args_parallel_reader()

    call mr_io_read_parallel_spacetime(MPI_COMM_WORLD, MPI_INFO_NULL, mr_io_test_mpi_cart_dims, path, mri_dest)

    print *, SpaceTimeMRI_group_name

    print *, shape(mri_dest%t_coordinates)
    print *, mri_dest%t_coordinates

    print *, shape(mri_dest%x_coordinates)
    print *, mri_dest%x_coordinates

    print *, shape(mri_dest%y_coordinates)
    print *, mri_dest%y_coordinates

    print *, shape(mri_dest%z_coordinates)
    print *, mri_dest%z_coordinates

    voxel_feature_shape = shape(mri_dest%voxel_feature%array)

    print *, mri_dest%voxel_feature%dims
    print *, mri_dest%voxel_feature%offset
    print *, voxel_feature_shape(3:5)

    print *, mri_dest%voxel_feature%time_dim
    print *, mri_dest%voxel_feature%time_offset
    print *, voxel_feature_shape(2)

    print *, voxel_feature_shape(1)
    print *, mri_dest%voxel_feature%array

    call MPI_Finalize(err)
        
end program mr_io_test_parallel_reader_space_time
