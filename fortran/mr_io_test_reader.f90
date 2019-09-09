program mr_io_test_reader
  
    use mr_io
    use mr_protocol

    use mr_io_test_arg_parser

    implicit none

!    character(len=13) :: path = "mr_io_test.h5"
    type(SpatialMRI) :: mri_dest    

    call mr_io_test_parse_args_reader()

    call mr_io_read_spatial(path, mri_dest)

    print *, SpatialMRI_group_name
    print *, mri_dest%voxel_feature_dims    
    print *, mri_dest%voxel_feature

end program mr_io_test_reader
