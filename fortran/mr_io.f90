module mr_io

use hdf5
use mr_protocol

implicit none

contains

subroutine mr_io_read_hdf5(path, mri_inst)     
  character(len=*), intent(in) :: path
  type(SpatialMRI), intent(out) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HSIZE_T), DIMENSION(3) :: dims = (/-1, -1, -1/) ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(3) :: max_dims = (/-1, -1, -1/) ! Dataset dimensions

  INTEGER     ::   rank = 3                         ! Dataset rank
  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)

  ! Open an existing file.
  CALL h5fopen_f (trim(path), H5F_ACC_RDWR_F, file_id, error)

  ! Open an existing group
  CALL h5gopen_f(file_id, SpatialMRI_group_name, grp_id, error)

  ! Open an existing dataset.
  CALL h5dopen_f(grp_id, "voxel_feature", dset_id, error)

  ! Read dims by getting dimension from data space
  CALL h5dget_space_f(dset_id, dspace_id, error)
  CALL h5sget_simple_extent_dims_f(dspace_id, dims, max_dims, error)

  ! Allocate MRI array 
  mri_inst%voxel_feature_dims = dims
  allocate(mri_inst%voxel_feature(mri_inst%voxel_feature_dims(1), mri_inst%voxel_feature_dims(2), mri_inst%voxel_feature_dims(3)))

  ! Read the dataset.
  CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, mri_inst%voxel_feature, dims, error)

  ! Close the data space
  CALL h5sclose_f(dspace_id, error)

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)

  ! Close the group
  CALL h5gclose_f(grp_id, error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)

end subroutine mr_io_read_hdf5


subroutine mr_io_write_hdf5(path, mri_inst)
  character(len=*), intent(in) :: path
  type(SpatialMRI), intent(in) :: mri_inst
  
  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HSIZE_T), DIMENSION(3) :: dims            ! Dataset dimensions
  INTEGER     ::   rank = 3                         ! Dataset rank

  INTEGER     ::   error ! Error flag

  dims = mri_inst%voxel_feature_dims

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)

  ! Create a new file using default properties.
  CALL h5fcreate_f(trim(path), H5F_ACC_TRUNC_F, file_id, error)

  ! Create a new group
  CALL h5gcreate_f(file_id, SpatialMRI_group_name, grp_id, error)

  ! Create the dataspace.
  CALL h5screate_simple_f(rank, dims, dspace_id, error)

  ! Create the dataset with default properties.
  CALL h5dcreate_f(grp_id, "voxel_feature", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)

  ! Write the dataset.
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mri_inst%voxel_feature, dims, error)
  
  ! End access to the dataset and release resources used by it.
  CALL h5dclose_f(dset_id, error)

  ! Terminate access to the data space.
  CALL h5sclose_f(dspace_id, error)

  ! Close the group
  CALL h5gclose_f(grp_id, error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)

end subroutine


subroutine mr_io_read_space_time(path, mri_inst)
  character(len=*), intent(in) :: path
  type(SpaceTimeMRI), intent(out) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HSIZE_T), DIMENSION(1) :: coord_dim = (/-1/) ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(1) :: max_coord_dim = (/-1/) ! Dataset dimensions

  INTEGER(HSIZE_T), DIMENSION(5) :: dims = (/-1, -1, -1, -1, -1/) ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(5) :: max_dims = (/-1, -1, -1, -1, -1/) ! Dataset dimensions

  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)

  ! Open an existing file.
  CALL h5fopen_f (trim(path), H5F_ACC_RDWR_F, file_id, error)

  ! Open an existing group
  CALL h5gopen_f(file_id, SpaceTimeMRI_group_name, grp_id, error)


  ! Read geometry data

  ! Open an existing dataset
  CALL h5dopen_f(grp_id, "x_coordinates", dset_id, error)

  ! Read dims by getting dimension from data space
  CALL h5dget_space_f(dset_id, dspace_id, error)
  CALL h5sget_simple_extent_dims_f(dspace_id, coord_dim, max_coord_dim, error)

  ! Allocate MRI array 
  mri_inst%x_dim = coord_dim(1)
  allocate(mri_inst%x_coordinates(mri_inst%x_dim))

  ! Read the dataset.
  CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, mri_inst%x_coordinates, coord_dim, error)

  ! Close the data space
  CALL h5sclose_f(dspace_id, error)

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)

  ! Open an existing dataset
  CALL h5dopen_f(grp_id, "y_coordinates", dset_id, error)

  ! Read dims by getting dimension from data space
  CALL h5dget_space_f(dset_id, dspace_id, error)
  CALL h5sget_simple_extent_dims_f(dspace_id, coord_dim, max_coord_dim, error)

  ! Allocate MRI array 
  mri_inst%y_dim = coord_dim(1)
  allocate(mri_inst%y_coordinates(mri_inst%y_dim))

  ! Read the dataset.
  CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, mri_inst%y_coordinates, coord_dim, error)

  ! Close the data space
  CALL h5sclose_f(dspace_id, error)

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)

  ! Open an existing dataset
  CALL h5dopen_f(grp_id, "z_coordinates", dset_id, error)

  ! Read dims by getting dimension from data space
  CALL h5dget_space_f(dset_id, dspace_id, error)
  CALL h5sget_simple_extent_dims_f(dspace_id, coord_dim, max_coord_dim, error)

  ! Allocate MRI array 
  mri_inst%z_dim = coord_dim(1)
  allocate(mri_inst%z_coordinates(mri_inst%z_dim))

  ! Read the dataset.
  CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, mri_inst%z_coordinates, coord_dim, error)

  ! Close the data space
  CALL h5sclose_f(dspace_id, error)

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)



  ! Read voxel_feature data

  ! Open an existing dataset
  CALL h5dopen_f(grp_id, "voxel_feature", dset_id, error)

  ! Read dims by getting dimension from data space
  CALL h5dget_space_f(dset_id, dspace_id, error)
  CALL h5sget_simple_extent_dims_f(dspace_id, dims, max_dims, error)

  ! Allocate MRI array 
  mri_inst%voxel_feature_dims = dims
  allocate(mri_inst%voxel_feature(mri_inst%voxel_feature_dims(1), &
                                  mri_inst%voxel_feature_dims(2), &
                                  mri_inst%voxel_feature_dims(3), &
                                  mri_inst%voxel_feature_dims(4), & 
                                  mri_inst%voxel_feature_dims(5) ))

  ! Read the dataset.
  CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, mri_inst%voxel_feature, dims, error)

  ! Close the data space
  CALL h5sclose_f(dspace_id, error)

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)



  ! Close the group
  CALL h5gclose_f(grp_id, error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)

end subroutine mr_io_read_space_time


subroutine mr_io_write_space_time(path, mri_inst)
  character(len=*), intent(in) :: path
  type(SpaceTimeMRI), intent(in) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
 
  INTEGER(HSIZE_T), DIMENSION(1) :: coord_dims            ! Coordinate dataset dimensions

  INTEGER(HSIZE_T), DIMENSION(5) :: dims            ! Dataset dimensions
  INTEGER     ::   rank = 5                         ! Dataset rank

  INTEGER     ::   error ! Error flag

  dims = mri_inst%voxel_feature_dims

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)

  ! Create a new file using default properties.
  CALL h5fcreate_f(trim(path), H5F_ACC_TRUNC_F, file_id, error)

  ! Create a new group
  CALL h5gcreate_f(file_id, SpaceTimeMRI_group_name, grp_id, error)


  ! Write geometry data
 
  coord_dims(1) =  mri_inst%x_dim

  ! Create the dataspace.
  CALL h5screate_simple_f(rank, coord_dims, dspace_id, error)

  ! Create the dataset with default properties.
  CALL h5dcreate_f(grp_id, "x_coordinates", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)

  ! Write the dataset.
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mri_inst%x_coordinates, coord_dims, error)

  ! End access to the dataset and release resources used by it.
  CALL h5dclose_f(dset_id, error)

  ! Terminate access to the data space.
  CALL h5sclose_f(dspace_id, error)

  coord_dims(1) =  mri_inst%y_dim
  
  ! Create the dataspace.
  CALL h5screate_simple_f(rank, coord_dims, dspace_id, error)

  ! Create the dataset with default properties.
  CALL h5dcreate_f(grp_id, "y_coordinates", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)

  ! Write the dataset.
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mri_inst%y_coordinates, coord_dims, error)

  ! End access to the dataset and release resources used by it.
  CALL h5dclose_f(dset_id, error)

  ! Terminate access to the data space.
  CALL h5sclose_f(dspace_id, error)

  coord_dims(1) =  mri_inst%z_dim

  ! Create the dataspace.
  CALL h5screate_simple_f(rank, coord_dims, dspace_id, error)

  ! Create the dataset with default properties.
  CALL h5dcreate_f(grp_id, "z_coordinates", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)

  ! Write the dataset.
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mri_inst%z_coordinates, coord_dims, error)

  ! End access to the dataset and release resources used by it.
  CALL h5dclose_f(dset_id, error)

  ! Terminate access to the data space.
  CALL h5sclose_f(dspace_id, error)



  ! Write voxel_feature data

  dims = mri_inst%voxel_feature_dims

  ! Create the dataspace.
  CALL h5screate_simple_f(rank, dims, dspace_id, error)

  ! Create the dataset with default properties.
  CALL h5dcreate_f(grp_id, "voxel_feature", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)

  ! Write the dataset.
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mri_inst%voxel_feature, dims, error)

  ! End access to the dataset and release resources used by it.
  CALL h5dclose_f(dset_id, error)

  ! Terminate access to the data space.
  CALL h5sclose_f(dspace_id, error)



  ! Close the group
  CALL h5gclose_f(grp_id, error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)

end subroutine


end module mr_io
