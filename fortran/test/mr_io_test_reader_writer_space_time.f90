program mr_io_test_reader_writer_space_time
  
    use mr_io

    use mr_io_test_arg_parser

    implicit none

!    character(len=100) :: in_path = "mr_io_test_space_time_in.h5"
!    character(len=100) :: out_path = "mr_io_test_space_time_out.h5"
    type(SpaceTimeMRI) :: mri_dest    

   call mr_io_test_parse_args_reader_writer()

    call mr_io_read_spacetime(in_path, mri_dest)
    call mr_io_write_spacetime(out_path, mri_dest)

end program mr_io_test_reader_writer_space_time
