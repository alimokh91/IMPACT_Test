program mr_io_test_reader_writer_hpc_predict
  
    use mr_io

    use mr_io_test_arg_parser

    implicit none

!    character(len=100) :: in_path = "mr_io_test_hpc_predict_in.h5"
!    character(len=100) :: out_path = "mr_io_test_hpc_predict_out.h5"
    type(HPCPredictMRI) :: mri_dest

    call mr_io_test_parse_args_reader_writer()

    call mr_io_read_hpcpredict(in_path, mri_dest)
    call mr_io_write_hpcpredict(out_path, mri_dest)

end program mr_io_test_reader_writer_hpc_predict
