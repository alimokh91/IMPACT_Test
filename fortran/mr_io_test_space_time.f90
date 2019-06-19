program mr_io_test_space_time
  
  use mr_io

  implicit none

    character(len=30) :: path = "test_space_time.h5 "
    type(SpaceTimeMRI) :: mri_inst
    real, dimension(4) :: xa
    real, dimension(3) :: ya
    real, dimension(2) :: za 
    real, dimension(3, 2, 4, 3, 2) :: a

    type(SpaceTimeMRI) :: mri_dest
    real, dimension(4) :: xb
    real, dimension(3) :: yb
    real, dimension(2) :: zb 
    real, dimension(3, 2, 4, 3, 2) :: b
    
    real, dimension(144) :: arr

    xa = (/ 1, 1, 1, 1/)
    ya = (/ 2, 2, 2/)
    za = (/ 3, 3 /)
    arr =  (/ 2, 4, 1, 8, 5, 0, 2, 7, 6, 5, 1, 8, 4, 5, 1, 7, 1, 0, 3, 8, &
                    5, 8, 3, 4, 4, 3, 0, 3, 9, 4, 5, 5, 4, 5, 9, 8, 3, &
                    9, 5, 0, 7, 7, 5, 5, 7, 6, 8, 3, 1, 6, 7, 0, 9, 6, &
                    2, 0, 7, 6, 1, 4, 1, 1, 6, 5, 6, 2, 1, 2, 3, 5, 3, 0, 3, 3, &
                    9, 9, 8, 6, 4, 2, 6, 7, 6, 6, 1, 5, 0, 6, 5, 8, 5, &
                    9, 0, 9, 1, 7, 4, 7, 6, 1, 9, 4, 4, 6, 3, 3, 9, 2, 4, 3, 0, & 
                    0, 9, 0, 1, 0, 9, 1, 3, 7, 9, 5, 9, 7, 3, 0, 2, 4, 2, 4, 3, &
                    8, 7, 7, 8, 0, 0, 0, 2, 4, 8, 5, 6, 1 &
    /)

    a = reshape( arr, (/ 3, 2, 4, 3, 2 /) )


    xb = (/ 4, 4, 4, 4 /)
    yb = (/ 5, 5, 5 /)
    zb = (/ 6, 6 /)

    arr = (/ 4, 0, 8, 7, 1, 6, 3, 1, 2, 7, 6, 9, 4, 7, 1, 9, 2, 4, 5, 4, 4, &
                    1, 0, 4, 6, 9, 5, 9, 4, 0, 5, 6, 5, 1, 8, 5, 7, 5, &
                    3, 1, 2, 5, 1, 0, 6, 1, 9, 8, 1, 9, 8, 3, 4, 5, 7, 3, &
                    3, 3, 2, 1, 7, 4, 5, 8, 3, 2, 4, 9, 5, 8, 4, 0, 5, 9, 5, 8, &
                    1, 7, 8, 9, 1, 9, 9, 8, 5, 8, 5, 1, 6, 7, 6, 0, 8, 0, 7, 7, &
                    4, 0, 8, 4, 4, 9, 2, 2, 0, 4, 6, 0, 2, 9, 4, 2, 1, 1, &
                    3, 3, 4, 7, 1, 6, 8, 0, 8, 6, 0, 5, 7, 5, 5, 4, 0, 6, &
                    0, 7, 7, 0, 7, 5, 6, 7, 3, 3, 0, 2 &
    /)

    b = reshape( arr, (/ 3, 2, 4, 3, 2 /) )


    mri_inst = SpaceTimeMRI( xa, 4, ya, 3, za, 2, a,  (/ 3, 2, 4, 3, 2 /))
    mri_dest = SpaceTimeMRI( xb, 4, yb, 3, zb, 2, b,  (/ 3, 2, 4, 3, 2 /))

    call mr_io_write_spacetime(path, mri_inst)
    call mr_io_read_spacetime(path, mri_dest)

    print *, mri_inst%x_coordinates
    print *, mri_dest%x_coordinates
    print *, mri_inst%y_coordinates
    print *, mri_dest%y_coordinates
    print *, mri_inst%z_coordinates
    print *, mri_dest%z_coordinates

    print *, mri_dest%voxel_feature_dims    
    print *, mri_dest%voxel_feature

end program mr_io_test_space_time
