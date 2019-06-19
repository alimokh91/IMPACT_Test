program mr_io_test
  
  use mr_io

  implicit none

    character(len=13) :: path = "test_file.h5 "
    type(SpatialMRI) :: mri_inst 
    real, dimension(2, 2, 2) :: a

    type(SpatialMRI) :: mri_dest
    real, dimension(2, 2, 2) :: b
    
    a = reshape( (/ &
      1, 2,         &
      3, 4,         &
      5, 6,         &
      7, 8          &
    /), (/ 2, 2, 2 /) )

    b = reshape( (/ &
      9, 10,        &
      11,12,        &
      13,14,        &
      15,16         &
    /), (/ 2, 2, 2 /) )


    mri_inst = SpatialMRI(a, (/2,2,2/))
    mri_dest = SpatialMRI(b, (/2,2,2/))

    call mr_io_write_spatial(path, mri_inst)
    call mr_io_read_spatial(path, mri_dest)

    print *, mri_dest%voxel_feature_dims    
    print *, mri_dest%voxel_feature

end program mr_io_test
