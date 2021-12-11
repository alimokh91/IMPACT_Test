program mr_io_test_parallel_reader_space_time_container
  
    use mr_io_parallel_spacetime
    use mr_io_locking_utils
    use mr_io_test_arg_parser
    use mpi

    implicit none

    INTEGER err
!    character(len=40) :: path = "mr_io_test_parallel_space_time.h5"
    type(DistSpacetimeMRI) :: mri_dest

    integer, dimension(5) :: vector_feature_shape
    integer :: st

    ! TODO: lock output file:
    st = mr_io_lock_stdout_stderr()

    call MPI_Init(err)
    
    call mr_io_test_parse_args_parallel_reader()

    ! TODO: with locking (C-version)
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

    vector_feature_shape = shape(mri_dest%vector_feature%array)

    print *, mri_dest%vector_feature%dims
    print *, mri_dest%vector_feature%offset
    print *, vector_feature_shape(3:5)

    print *, mri_dest%vector_feature%time_dim
    print *, mri_dest%vector_feature%time_offset
    print *, vector_feature_shape(2)

    print *, vector_feature_shape(1)
    print *, mri_dest%vector_feature%array

    ! TODO: flush & unlock output file
    st = mr_io_unlock_stdout_stderr()

    call MPI_Finalize(err)
        
end program mr_io_test_parallel_reader_space_time_container
