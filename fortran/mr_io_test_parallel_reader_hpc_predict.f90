program mr_io_test_parallel_reader_space_time
  
    use mr_io_parallel_spacetime
    use mr_protocol
    use mr_io_test_arg_parser

    implicit none

    INTEGER err
!    character(len=100) :: path = "mr_io_test_parallel_hpc_predict.h5"
    type(DistHPCPredictMRI) :: mri_dest

    integer, dimension(4) :: intensity_shape
    integer, dimension(5) :: velocity_mean_shape
    integer, dimension(6) :: velocity_cov_shape

    call MPI_Init(err)
    
    call mr_io_test_parse_args_parallel_reader()

    call mr_io_read_parallel_hpcpredict(MPI_COMM_WORLD, MPI_INFO_NULL, mr_io_test_mpi_cart_dims, path, mri_dest)

    print *, HPCPredictMRI_group_name

    print *, shape(mri_dest%t_coordinates)
    print *, mri_dest%t_coordinates

    print *, shape(mri_dest%x_coordinates)
    print *, mri_dest%x_coordinates

    print *, shape(mri_dest%y_coordinates)
    print *, mri_dest%y_coordinates

    print *, shape(mri_dest%z_coordinates)
    print *, mri_dest%z_coordinates

    intensity_shape = shape(mri_dest%intensity%array)

    print *, mri_dest%intensity%dims
    print *, mri_dest%intensity%offset
    print *, intensity_shape(2:4)

    print *, mri_dest%intensity%time_dim
    print *, mri_dest%intensity%time_offset
    print *, intensity_shape(1)

    print *, mri_dest%intensity%array


    velocity_mean_shape = shape(mri_dest%velocity_mean%array)

    print *, mri_dest%velocity_mean%dims
    print *, mri_dest%velocity_mean%offset
    print *, velocity_mean_shape(3:5)

    print *, mri_dest%velocity_mean%time_dim
    print *, mri_dest%velocity_mean%time_offset
    print *, velocity_mean_shape(2)

    print *, velocity_mean_shape(1)
    print *, mri_dest%velocity_mean%array


    velocity_cov_shape = shape(mri_dest%velocity_cov%array)

    print *, mri_dest%velocity_cov%dims
    print *, mri_dest%velocity_cov%offset
    print *, velocity_cov_shape(4:6)

    print *, mri_dest%velocity_cov%time_dim
    print *, mri_dest%velocity_cov%time_offset
    print *, velocity_cov_shape(3)

    print *, velocity_cov_shape(1:2)
    print *, mri_dest%velocity_cov%array


    call MPI_Finalize(err)
        
end program mr_io_test_parallel_reader_space_time
