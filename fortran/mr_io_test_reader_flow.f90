program mr_io_test_reader_flow
  
    use mr_io

    use mr_io_test_arg_parser

    implicit none

!    character(len=100) :: path = "mr_io_test_flow.h5"
    type(FlowMRI) :: mri_dest

    call mr_io_test_parse_args_reader()

    call mr_io_read_flow(path, mri_dest)

    print *, FlowMRI_group_name

    print *, mri_dest%x_dim
    print *, mri_dest%x_coordinates

    print *, mri_dest%y_dim
    print *, mri_dest%y_coordinates

    print *, mri_dest%z_dim
    print *, mri_dest%z_coordinates

    print *, mri_dest%t_dim
    print *, mri_dest%t_coordinates

    print *, mri_dest%intensity_dims
    print *, mri_dest%intensity

    print *, mri_dest%velocity_mean_dims
    print *, mri_dest%velocity_mean

    print *, mri_dest%velocity_cov_dims
    print *, mri_dest%velocity_cov

end program mr_io_test_reader_flow
