program mr_io_test_reader
  
    use mr_io
    use mr_protocol

    implicit none

    character(len=24) :: path = "mr_io_test_space_time.h5"
    type(SpaceTimeMRI) :: mri_dest    

    call mr_io_read_space_time(path, mri_dest)

    print *, SpaceTimeMRI_group_name
    print *, mri_dest%voxel_feature_dims    
    print *, mri_dest%voxel_feature

end program mr_io_test_reader
