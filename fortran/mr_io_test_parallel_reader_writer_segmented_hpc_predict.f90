program mr_io_test_parallel_reader_writer_segmented_hpc_predict
  
    use mr_io_parallel_spacetime
    use mr_io_test_arg_parser

    !include 'mpif.h'

    implicit none

    INTEGER err
!    character(len=100) :: in_path = "mr_io_test_parallel_segmented_hpc_predict_in.h5"
!    character(len=100) :: out_path = "mr_io_test_parallel_segmented_hpc_predict_out.h5"
    type(DistSegmentedHPCPredictMRI) :: mri_dest

    integer, dimension(5) :: velocity_mean_shape
    integer, dimension(6) :: velocity_cov_shape

    call MPI_Init(err)
    
    call mr_io_test_parse_args_parallel_reader_writer()

    call mr_io_read_parallel_segmentedhpcpredict(MPI_COMM_WORLD, MPI_INFO_NULL, mpi_cart_dims, in_path, mri_dest)

    print *, SegmentedHPCPredictMRI_group_name

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

    call mr_io_write_parallel_segmentedhpcpredict(MPI_COMM_WORLD, MPI_INFO_NULL, out_path, mri_dest)

    call MPI_Finalize(err)

end program mr_io_test_parallel_reader_writer_segmented_hpc_predict

