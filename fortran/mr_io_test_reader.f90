program mr_io_test_reader
  
  use mr_io

  implicit none

    character(len=13) :: path = "mr_io_test.h5 "
    type(mri) :: mri_dest    

    call mr_io_read_hdf5(path, mri_dest)

    print *, mri_dest%voxel_feature_dims    
    print *, mri_dest%voxel_feature

end program mr_io_test_reader
