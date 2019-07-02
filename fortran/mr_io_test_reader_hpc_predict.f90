program mr_io_test_reader_hpc_predict
  
    use mr_io
    use mr_protocol

    implicit none

    character(len=100) :: path = "mr_io_test_hpc_predict.h5"
    type(HPCPredictMRI) :: mri_dest

    call mr_io_read_hpcpredict(path, mri_dest)

    print *, HPCPredictMRI_group_name

    print *, mri_dest%x_dim
    print *, mri_dest%x_coordinates

    print *, mri_dest%y_dim
    print *, mri_dest%y_coordinates

    print *, mri_dest%z_dim
    print *, mri_dest%z_coordinates

    print *, mri_dest%t_dim
    print *, mri_dest%t_coordinates

    print *, mri_dest%velocity_mean_dims
    print *, mri_dest%velocity_mean

    print *, mri_dest%velocity_cov_dims
    print *, mri_dest%velocity_cov

end program mr_io_test_reader_hpc_predict
