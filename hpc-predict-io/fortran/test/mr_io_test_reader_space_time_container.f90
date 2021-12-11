program mr_io_test_reader_space_time
  
    use mr_io
    use mr_io_locking_utils

    use mr_io_test_arg_parser

    implicit none

!    character(len=24) :: path = "mr_io_test_space_time.h5"
    type(SpaceTimeMRI) :: mri_dest    
    integer :: st

    ! TODO: lock output file:
    st = mr_io_lock_stdout_stderr()

    call mr_io_test_parse_args_reader()

    ! TODO: with locking (C-version)
    call mr_io_read_spacetime(path, mri_dest)

    print *, SpaceTimeMRI_group_name

    print *, mri_dest%x_dim
    print *, mri_dest%x_coordinates

    print *, mri_dest%y_dim
    print *, mri_dest%y_coordinates

    print *, mri_dest%z_dim
    print *, mri_dest%z_coordinates

    print *, mri_dest%t_dim
    print *, mri_dest%t_coordinates

    print *, mri_dest%vector_feature_dims    
    print *, mri_dest%vector_feature

    ! TODO: flush & unlock output file
    st = mr_io_unlock_stdout_stderr()

end program mr_io_test_reader_space_time
