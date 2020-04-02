program mr_io_test_parallel_reader
  
    use mr_io_parallel
    use mr_io_locking_utils
    use mr_io_test_arg_parser
    use mpi

    implicit none

    integer err
!    character(len=100) :: path != "mr_io_test_parallel.h5"
    type(DistSpatialMRI) :: mri_dest    
    integer :: st

    ! TODO: lock output file:
    st = mr_io_lock_stdout_stderr()

    call MPI_Init(err)
    
    call mr_io_test_parse_args_parallel_reader()

    ! TODO: with locking (C-version)
    call mr_io_read_parallel_spatial(MPI_COMM_WORLD, MPI_INFO_NULL, mr_io_test_mpi_cart_dims, &
                                     path, mri_dest)

    print *, SpatialMRI_group_name
    print *, mri_dest%scalar_feature%dims
    print *, mri_dest%scalar_feature%offset
    print *, shape(mri_dest%scalar_feature%array)
    print *, mri_dest%scalar_feature%array

    ! TODO: flush & unlock output file
    st = mr_io_unlock_stdout_stderr()

    call MPI_Finalize(err)
        
end program mr_io_test_parallel_reader

