program mr_io_test_container
  
    use mr_io
    use mr_io_locking_utils

    use mr_io_test_arg_parser

    implicit none

!    character(len=13) :: path = "mr_io_test.h5"
    type(SpatialMRI) :: mri_dest
    integer :: st

    ! TODO: lock output file:
    st = mr_io_lock_stdout_stderr()

    call mr_io_test_parse_args_reader()

    ! TODO: with locking (C-version)
    call mr_io_read_spatial(path, mri_dest)

    print *, SpatialMRI_group_name
    print *, mri_dest%scalar_feature_dims
    print *, mri_dest%scalar_feature

    ! TODO: flush & unlock output file
    st = mr_io_unlock_stdout_stderr()

end program mr_io_test_container
