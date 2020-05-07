program mr_io_test_reader_writer_flow
  
    use mr_io

    use mr_io_test_arg_parser

    implicit none

!    character(len=100) :: in_path = "mr_io_test_flow_in.h5"
!    character(len=100) :: out_path = "mr_io_test_flow_out.h5"
    type(FlowMRI) :: mri_dest

!    INTEGER :: gdb = 0
!    do while (gdb == 0)
!        call sleep(2)
!    end do
!
!    ! End of debugging

    call mr_io_test_parse_args_reader_writer()

    call mr_io_read_flow(in_path, mri_dest)
    call mr_io_write_flow(out_path, mri_dest)

end program mr_io_test_reader_writer_flow
