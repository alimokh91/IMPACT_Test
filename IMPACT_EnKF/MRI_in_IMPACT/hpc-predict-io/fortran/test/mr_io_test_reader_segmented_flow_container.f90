program mr_io_test_reader_segmented_flow
  
    use mr_io
    use mr_io_locking_utils

    use mr_io_test_arg_parser

    implicit none

!    character(len=100) :: path = "mr_io_test_segmented_flow.h5"
    type(SegmentedFlowMRI) :: mri_dest
    integer :: st

    ! TODO: lock output file:
    st = mr_io_lock_stdout_stderr()

    call mr_io_test_parse_args_reader()

    ! TODO: with locking (C-version)
    call mr_io_read_segmentedflow(path, mri_dest)

    print *, SegmentedFlowMRI_group_name

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

    print *, mri_dest%segmentation_prob_dims
    print *, mri_dest%segmentation_prob

    ! TODO: flush & unlock output file
    st = mr_io_unlock_stdout_stderr()

end program mr_io_test_reader_segmented_flow
