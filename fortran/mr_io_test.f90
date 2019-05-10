program mr_io_test
  
  use mr_io

  implicit none

    character(len=13) :: path = "test_file.h5 "
    type(mri) :: mri_inst 
    real, dimension(2, 2) :: a

    type(mri) :: mri_dest
    real, dimension(2, 2) :: b
    
    a = reshape( (/ &
      1, 2,         &
      3, 4          &
    /), (/ 2, 2 /) )

    b = reshape( (/ &
      7, 8,         &
      9, 10         &
    /), (/ 2, 2 /) )


    mri_inst = mri(a, (/2,2/))
    

    call mr_io_write_hdf5(path, mri_inst)
    call mr_io_read_hdf5(path, mri_dest)

    print *, mri_inst%voxel_feature_dims    
    print *, mri_inst%voxel_feature

end program mr_io_test
