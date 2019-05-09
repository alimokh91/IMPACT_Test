module mr_io

USE HDF5 ! This module contains all necessary modules

implicit none

type MRI
     integer, dimension(2,2) :: voxel_feature
end type

integer(HSIZE_T), dimension(2) :: voxel_feature_dims = (/2,2/)

contains


subroutine mr_io_read_hdf5(path, mri_inst)     
  character(len=*), intent(in) :: path
  type(MRI), intent(out) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HSIZE_T), DIMENSION(2) :: dims = (/2, 2/) ! Dataset dimensions
  INTEGER     ::   rank = 2                         ! Dataset rank

  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)

  ! Open an existing file.
  CALL h5fopen_f (trim(path), H5F_ACC_RDWR_F, file_id, error)

  ! Open an existing dataset.
  CALL h5dopen_f(file_id, "voxel_feature", dset_id, error)

  ! Read the dataset.
  CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, mri_inst%voxel_feature, voxel_feature_dims, error)

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)

end subroutine mr_io_read_hdf5


subroutine mr_io_write_hdf5(path, mri_inst)
  character(len=*), intent(in) :: path
  type(MRI), intent(in) :: mri_inst
  
  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

!  INTEGER, DIMENSION(2,2) :: dset_data              ! Data buffers  
  INTEGER(HSIZE_T), DIMENSION(2) :: dims = (/2, 2/) ! Dataset dimensions
  INTEGER     ::   rank = 2                         ! Dataset rank

  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)

  ! Create a new file using default properties.
  CALL h5fcreate_f(trim(path), H5F_ACC_TRUNC_F, file_id, error)

  ! Create the dataspace.
  CALL h5screate_simple_f(rank, dims, dspace_id, error)

  ! Create the dataset with default properties.
  CALL h5dcreate_f(file_id, "voxel_feature", H5T_NATIVE_INTEGER, dspace_id, &
       dset_id, error)

  ! Write the dataset.
  CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, mri_inst%voxel_feature, voxel_feature_dims, error)
  
  ! End access to the dataset and release resources used by it.
  CALL h5dclose_f(dset_id, error)

  ! Terminate access to the data space.
  CALL h5sclose_f(dspace_id, error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)

end subroutine

end module mr_io
