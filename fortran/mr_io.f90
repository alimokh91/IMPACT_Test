module mr_io

use hdf5
use mr_protocol

!include 'mpif.h'

! TODO: private/public API

contains


! ************************ error handling ************************

#ifdef __GFORTRAN__

#define mr_io_handle_error(err) if(err /= 0) CALL mr_io_handle_hdf5_error(err);
#define mr_io_handle_arg_error(err) if(err /= 0) CALL mr_io_handle_argument_error(err);
#define mr_io_sget_simple_extent_dims_handle_error(err) if(err == -1) CALL mr_io_handle_hdf5_error(err);

subroutine mr_io_handle_hdf5_error(error)

  implicit none

  INTEGER     ::   error          ! Error flag

  if(error /= 0) then
    write(*,*) "HDF5 error ",error,"- printing backtrace and aborting..."
    call h5eprint_f(error)
    call backtrace
    !call abort
  end if

end subroutine mr_io_handle_hdf5_error


subroutine mr_io_handle_argument_error(error)

  implicit none

  INTEGER     ::   error          ! Error flag

  if(error /= 0) then
    write(*,*) "Argument error ",error,"- printing backtrace and aborting..."
    call backtrace
    !call abort
  end if

end subroutine mr_io_handle_argument_error

#endif


! ************************ SpatialMRI ************************

subroutine mr_io_read_spatial_feature(grp_id, feature_name, feature_array)

  implicit none

  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name
  real*8, dimension(:,:,:), allocatable, intent(out) :: feature_array

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HSIZE_T), DIMENSION(3) :: dims = (/-1, -1, -1/) ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(3) :: max_dims = (/-1, -1, -1/) ! Dataset dimensions

  INTEGER     ::   rank = 3       ! Dataset rank
  INTEGER     ::   error          ! Error flag
  
  ! Open an existing dataset.
  CALL h5dopen_f(grp_id, feature_name, dset_id, error)
  mr_io_handle_error(error)

  ! Read dims by getting dimension from data space
  CALL h5dget_space_f(dset_id, dspace_id, error)
  mr_io_handle_error(error)
  CALL h5sget_simple_extent_dims_f(dspace_id, dims, max_dims, error)
  mr_io_sget_simple_extent_dims_handle_error(error)

  ! Allocate MRI array 
  allocate(feature_array(dims(1), dims(2), dims(3)))

  ! Read the dataset.
  CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, feature_array, dims, error)
  mr_io_handle_error(error)

  ! Close the data space
  CALL h5sclose_f(dspace_id, error)
  mr_io_handle_error(error)

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_read_spatial_feature


subroutine mr_io_write_spatial_feature(grp_id, feature_name, feature_array)

  implicit none

  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name
  real*8, dimension(:,:,:), intent(in) :: feature_array

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HSIZE_T), DIMENSION(3) :: dims = (/-1, -1, -1/) ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(3) :: max_dims = (/-1, -1, -1/) ! Dataset dimensions

  INTEGER     ::   rank = 3       ! Dataset rank
  INTEGER     ::   error          ! Error flag

  dims = shape(feature_array)

  ! Create the dataspace.
  CALL h5screate_simple_f(rank, dims, dspace_id, error)
  mr_io_handle_error(error)

  ! Create the dataset with default properties.
  CALL h5dcreate_f(grp_id, feature_name, H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  mr_io_handle_error(error)

  ! Write the dataset.
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, feature_array, dims, error)
  mr_io_handle_error(error)

  ! End access to the dataset and release resources used by it.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

  ! Terminate access to the data space.
  CALL h5sclose_f(dspace_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_write_spatial_feature


! Spatial MRI I/O

subroutine mr_io_read_spatial(path, mri_inst)

  implicit none

  character(len=*), intent(in) :: path
  type(SpatialMRI), intent(out) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier

  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  mr_io_handle_error(error)

  ! Open an existing file.
  CALL h5fopen_f (trim(path), H5F_ACC_RDWR_F, file_id, error)
  mr_io_handle_error(error)

  ! Open an existing group
  CALL h5gopen_f(file_id, SpatialMRI_group_name, grp_id, error)
  mr_io_handle_error(error)

  ! Read spatial feature
  CALL mr_io_read_spatial_feature(grp_id, "voxel_feature", mri_inst%voxel_feature)
  mri_inst%voxel_feature_dims = shape(mri_inst%voxel_feature)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_read_spatial


subroutine mr_io_write_spatial(path, mri_inst)

  implicit none

  character(len=*), intent(in) :: path
  type(SpatialMRI), intent(in) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier

  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  mr_io_handle_error(error)

  ! Create a new file using default properties.
  CALL h5fcreate_f(trim(path), H5F_ACC_TRUNC_F, file_id, error)
  mr_io_handle_error(error)

  ! Create a new group
  CALL h5gcreate_f(file_id, SpatialMRI_group_name, grp_id, error)
  mr_io_handle_error(error)

  ! Write spatial feature
  CALL mr_io_write_spatial_feature(grp_id, "voxel_feature", mri_inst%voxel_feature)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_write_spatial



! ************************ SpaceTimeMRI ************************

subroutine mr_io_read_spacetime_feature(grp_id, feature_name, feature_array)

  implicit none

  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name
  real*8, dimension(:,:,:,:,:), allocatable, intent(out) :: feature_array

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HSIZE_T), DIMENSION(5) :: dims = (/-1, -1, -1, -1, -1/) ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(5) :: max_dims = (/-1, -1, -1, -1, -1/) ! Dataset dimensions

  INTEGER     ::   error          ! Error flag

  ! Open an existing dataset.
  CALL h5dopen_f(grp_id, feature_name, dset_id, error)
  mr_io_handle_error(error)

  ! Read dims by getting dimension from data space
  CALL h5dget_space_f(dset_id, dspace_id, error)
  mr_io_handle_error(error)
  CALL h5sget_simple_extent_dims_f(dspace_id, dims, max_dims, error)
  mr_io_sget_simple_extent_dims_handle_error(error)

  ! Allocate MRI array 
  allocate(feature_array(dims(1), dims(2), dims(3), dims(4), dims(5)))

  ! Read the dataset.
  CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, feature_array, dims, error)
  mr_io_handle_error(error)

  ! Close the data space
  CALL h5sclose_f(dspace_id, error)
  mr_io_handle_error(error)

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)
  
end subroutine


subroutine mr_io_write_spacetime_feature(grp_id, feature_name, feature_array)

  implicit none

  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name
  real*8, dimension(:,:,:,:,:), intent(in) :: feature_array

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HSIZE_T), DIMENSION(5) :: dims = (/-1, -1, -1, -1, -1/) ! Dataset dimensions

  INTEGER     ::   rank = 5       ! Dataset rank
  INTEGER     ::   error          ! Error flag

  dims = shape(feature_array)

  ! Create the dataspace.
  CALL h5screate_simple_f(rank, dims, dspace_id, error)
  mr_io_handle_error(error)

  ! Create the dataset with default properties.
  CALL h5dcreate_f(grp_id, feature_name, H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  mr_io_handle_error(error)

  ! Write the dataset.
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, feature_array, dims, error)
  mr_io_handle_error(error)

  ! End access to the dataset and release resources used by it.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

  ! Terminate access to the data space.
  CALL h5sclose_f(dspace_id, error)
  mr_io_handle_error(error)

end subroutine


! Coordinate I/O

subroutine mr_io_read_coordinates(grp_id, coordinate, coordinate_array)

  implicit none

  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: coordinate
  real*8, dimension(:), allocatable, intent(out) :: coordinate_array

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HSIZE_T), DIMENSION(1) :: coord_dim = (/-1/) ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(1) :: max_coord_dim = (/-1/) ! Dataset dimensions

  INTEGER     ::   error ! Error flag

  ! Open an existing dataset
  CALL h5dopen_f(grp_id, coordinate//"_coordinates", dset_id, error)
  mr_io_handle_error(error)

  ! Read dims by getting dimension from data space
  CALL h5dget_space_f(dset_id, dspace_id, error)
  mr_io_handle_error(error)
  CALL h5sget_simple_extent_dims_f(dspace_id, coord_dim, max_coord_dim, error)
  mr_io_sget_simple_extent_dims_handle_error(error)

  ! Allocate coordinate array
  allocate(coordinate_array(coord_dim(1)))

  ! Read the dataset.
  CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, coordinate_array, coord_dim, error)
  mr_io_handle_error(error)

  ! Close the data space
  CALL h5sclose_f(dspace_id, error)
  mr_io_handle_error(error)

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_read_coordinates



subroutine mr_io_write_coordinates(grp_id, coordinate, coordinate_array)

  implicit none

  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: coordinate
  real*8, dimension(:), allocatable, intent(in) :: coordinate_array

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HSIZE_T), DIMENSION(1) :: coord_dim = (/-1/) ! Dataset dimensions
  INTEGER     ::   rank = 1

  INTEGER     ::   error ! Error flag

  coord_dim =  shape(coordinate_array)

  ! Create the dataspace.
  CALL h5screate_simple_f(rank, coord_dim, dspace_id, error)
  mr_io_handle_error(error)

  ! Create the dataset with default properties.
  CALL h5dcreate_f(grp_id, coordinate//"_coordinates", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  mr_io_handle_error(error)

  ! Write the dataset.
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, coordinate_array, coord_dim, error)
  mr_io_handle_error(error)

  ! End access to the dataset and release resources used by it.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

  ! Terminate access to the data space.
  CALL h5sclose_f(dspace_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_write_coordinates


! Space-time MRI I/O

subroutine mr_io_read_spacetime(path, mri_inst)

  implicit none

  character(len=*), intent(in) :: path
  type(SpaceTimeMRI), intent(out) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier

  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  mr_io_handle_error(error)

  ! Open an existing file.
  CALL h5fopen_f (trim(path), H5F_ACC_RDWR_F, file_id, error)
  mr_io_handle_error(error)

  ! Open an existing group
  CALL h5gopen_f(file_id, SpaceTimeMRI_group_name, grp_id, error)
  mr_io_handle_error(error)

  ! Read time data
  CALL mr_io_read_coordinates(grp_id, "t", mri_inst%t_coordinates)
  mri_inst%t_dim = size(mri_inst%t_coordinates)

  ! Read geometry data
  CALL mr_io_read_coordinates(grp_id, "x", mri_inst%x_coordinates)
  mri_inst%x_dim = size(mri_inst%x_coordinates)
  CALL mr_io_read_coordinates(grp_id, "y", mri_inst%y_coordinates)
  mri_inst%y_dim = size(mri_inst%y_coordinates)
  CALL mr_io_read_coordinates(grp_id, "z", mri_inst%z_coordinates)
  mri_inst%z_dim = size(mri_inst%z_coordinates)

  ! Read voxel_feature data
  CALL mr_io_read_spacetime_feature(grp_id, "voxel_feature", mri_inst%voxel_feature)
  mri_inst%voxel_feature_dims =  shape(mri_inst%voxel_feature)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_read_spacetime



subroutine mr_io_write_spacetime(path, mri_inst)

  implicit none

  character(len=*), intent(in) :: path
  type(SpaceTimeMRI), intent(in) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier
 
  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  mr_io_handle_error(error)

  ! Create a new file using default properties.
  CALL h5fcreate_f(trim(path), H5F_ACC_TRUNC_F, file_id, error)
  mr_io_handle_error(error)

  ! Create a new group
  CALL h5gcreate_f(file_id, SpaceTimeMRI_group_name, grp_id, error)
  mr_io_handle_error(error)

  ! Write time data
  CALL mr_io_write_coordinates(grp_id, "t", mri_inst%t_coordinates)

  ! Write geometry data
  CALL mr_io_write_coordinates(grp_id, "x", mri_inst%x_coordinates)
  CALL mr_io_write_coordinates(grp_id, "y", mri_inst%y_coordinates)
  CALL mr_io_write_coordinates(grp_id, "z", mri_inst%z_coordinates)

  ! Write voxel_feature data
  CALL mr_io_write_spacetime_feature(grp_id, "voxel_feature", mri_inst%voxel_feature)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_write_spacetime



! ************************ HPCPredictMRI ************************

subroutine mr_io_read_spacetime_scalar_feature(grp_id, feature_name, feature_array)

  implicit none

  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name
  real*8, dimension(:,:,:,:), allocatable, intent(out) :: feature_array

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HSIZE_T), DIMENSION(4) :: dims = (/-1, -1, -1, -1/) ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(4) :: max_dims = (/-1, -1, -1, -1/) ! Dataset dimensions

  INTEGER     ::   error          ! Error flag

  ! Open an existing dataset.
  CALL h5dopen_f(grp_id, feature_name, dset_id, error)
  mr_io_handle_error(error)

  ! Read dims by getting dimension from data space
  CALL h5dget_space_f(dset_id, dspace_id, error)
  mr_io_handle_error(error)
  CALL h5sget_simple_extent_dims_f(dspace_id, dims, max_dims, error)
  mr_io_sget_simple_extent_dims_handle_error(error)

  ! Allocate MRI array
  allocate(feature_array(dims(1), dims(2), dims(3), dims(4)))

  ! Read the dataset.
  CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, feature_array, dims, error)
  mr_io_handle_error(error)

  ! Close the data space
  CALL h5sclose_f(dspace_id, error)
  mr_io_handle_error(error)

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_read_spacetime_scalar_feature


subroutine mr_io_write_spacetime_scalar_feature(grp_id, feature_name, feature_array)

  implicit none

  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name
  real*8, dimension(:,:,:,:), intent(in) :: feature_array

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HSIZE_T), DIMENSION(4) :: dims = (/-1, -1, -1, -1/) ! Dataset dimensions

  INTEGER     ::   rank = 4       ! Dataset rank
  INTEGER     ::   error          ! Error flag

  dims = shape(feature_array)

  ! Create the dataspace.
  CALL h5screate_simple_f(rank, dims, dspace_id, error)
  mr_io_handle_error(error)

  ! Create the dataset with default properties.
  CALL h5dcreate_f(grp_id, feature_name, H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  mr_io_handle_error(error)

  ! Write the dataset.
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, feature_array, dims, error)
  mr_io_handle_error(error)

  ! End access to the dataset and release resources used by it.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

  ! Terminate access to the data space.
  CALL h5sclose_f(dspace_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_write_spacetime_scalar_feature



subroutine mr_io_read_spacetime_matrix_feature(grp_id, feature_name, feature_array)

  implicit none

  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name
  real*8, dimension(:,:,:,:,:,:), allocatable, intent(out) :: feature_array

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HSIZE_T), DIMENSION(6) :: dims = (/-1, -1, -1, -1, -1, -1/) ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(6) :: max_dims = (/-1, -1, -1, -1, -1, -1/) ! Dataset dimensions

  INTEGER     ::   error          ! Error flag

  ! Open an existing dataset.
  CALL h5dopen_f(grp_id, feature_name, dset_id, error)
  mr_io_handle_error(error)

  ! Read dims by getting dimension from data space
  CALL h5dget_space_f(dset_id, dspace_id, error)
  mr_io_handle_error(error)
  CALL h5sget_simple_extent_dims_f(dspace_id, dims, max_dims, error)
  mr_io_sget_simple_extent_dims_handle_error(error)

  ! Allocate MRI array
  allocate(feature_array(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))

  ! Read the dataset.
  CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, feature_array, dims, error)
  mr_io_handle_error(error)

  ! Close the data space
  CALL h5sclose_f(dspace_id, error)
  mr_io_handle_error(error)

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_read_spacetime_matrix_feature


subroutine mr_io_write_spacetime_matrix_feature(grp_id, feature_name, feature_array)

  implicit none

  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name
  real*8, dimension(:,:,:,:,:,:), intent(in) :: feature_array

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HSIZE_T), DIMENSION(6) :: dims = (/-1, -1, -1, -1, -1, -1/) ! Dataset dimensions

  INTEGER     ::   rank = 6       ! Dataset rank
  INTEGER     ::   error          ! Error flag

  dims = shape(feature_array)

  ! Create the dataspace.
  CALL h5screate_simple_f(rank, dims, dspace_id, error)
  mr_io_handle_error(error)

  ! Create the dataset with default properties.
  CALL h5dcreate_f(grp_id, feature_name, H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  mr_io_handle_error(error)

  ! Write the dataset.
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, feature_array, dims, error)
  mr_io_handle_error(error)

  ! End access to the dataset and release resources used by it.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

  ! Terminate access to the data space.
  CALL h5sclose_f(dspace_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_write_spacetime_matrix_feature


subroutine mr_io_read_hpcpredict(path, mri_inst)

  implicit none

  character(len=*), intent(in) :: path
  type(HPCPredictMRI), intent(out) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier

  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  mr_io_handle_error(error)

  ! Open an existing file.
  CALL h5fopen_f (trim(path), H5F_ACC_RDWR_F, file_id, error)
  mr_io_handle_error(error)

  ! Open an existing group
  CALL h5gopen_f(file_id, HPCPredictMRI_group_name, grp_id, error)
  mr_io_handle_error(error)

  ! Read time data
  CALL mr_io_read_coordinates(grp_id, "t", mri_inst%t_coordinates)
  mri_inst%t_dim = size(mri_inst%t_coordinates)

  ! Read geometry data
  CALL mr_io_read_coordinates(grp_id, "x", mri_inst%x_coordinates)
  mri_inst%x_dim = size(mri_inst%x_coordinates)
  CALL mr_io_read_coordinates(grp_id, "y", mri_inst%y_coordinates)
  mri_inst%y_dim = size(mri_inst%y_coordinates)
  CALL mr_io_read_coordinates(grp_id, "z", mri_inst%z_coordinates)
  mri_inst%z_dim = size(mri_inst%z_coordinates)

  ! Read voxel_feature data
  CALL mr_io_read_spacetime_scalar_feature(grp_id, "intensity", mri_inst%intensity)
  mri_inst%intensity_dims =  shape(mri_inst%intensity)
  CALL mr_io_read_spacetime_feature(grp_id, "velocity_mean", mri_inst%velocity_mean)
  mri_inst%velocity_mean_dims =  shape(mri_inst%velocity_mean)
  CALL mr_io_read_spacetime_matrix_feature(grp_id, "velocity_cov", mri_inst%velocity_cov)
  mri_inst%velocity_cov_dims =  shape(mri_inst%velocity_cov)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_read_hpcpredict



subroutine mr_io_write_hpcpredict(path, mri_inst)

  implicit none

  character(len=*), intent(in) :: path
  type(HPCPredictMRI), intent(in) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier

  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  mr_io_handle_error(error)

  ! Create a new file using default properties.
  CALL h5fcreate_f(trim(path), H5F_ACC_TRUNC_F, file_id, error)
  mr_io_handle_error(error)

  ! Create a new group
  CALL h5gcreate_f(file_id, HPCPredictMRI_group_name, grp_id, error)
  mr_io_handle_error(error)

  ! Write time data
  CALL mr_io_write_coordinates(grp_id, "t", mri_inst%t_coordinates)

  ! Write geometry data
  CALL mr_io_write_coordinates(grp_id, "x", mri_inst%x_coordinates)
  CALL mr_io_write_coordinates(grp_id, "y", mri_inst%y_coordinates)
  CALL mr_io_write_coordinates(grp_id, "z", mri_inst%z_coordinates)

  ! Write voxel_feature data
  CALL mr_io_write_spacetime_scalar_feature(grp_id, "intensity", mri_inst%intensity)
  CALL mr_io_write_spacetime_feature(grp_id, "velocity_mean", mri_inst%velocity_mean)
  CALL mr_io_write_spacetime_matrix_feature(grp_id, "velocity_cov", mri_inst%velocity_cov)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_write_hpcpredict


! ************************ SegmentedHPCPredictMRI ************************

subroutine mr_io_read_segmentedhpcpredict(path, mri_inst)

  implicit none

  character(len=*), intent(in) :: path
  type(SegmentedHPCPredictMRI), intent(out) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier

  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  mr_io_handle_error(error)

  ! Open an existing file.
  CALL h5fopen_f (trim(path), H5F_ACC_RDWR_F, file_id, error)
  mr_io_handle_error(error)

  ! Open an existing group
  CALL h5gopen_f(file_id, SegmentedHPCPredictMRI_group_name, grp_id, error)
  mr_io_handle_error(error)

  ! Read time data
  CALL mr_io_read_coordinates(grp_id, "t", mri_inst%t_coordinates)
  mri_inst%t_dim = size(mri_inst%t_coordinates)

  ! Read geometry data
  CALL mr_io_read_coordinates(grp_id, "x", mri_inst%x_coordinates)
  mri_inst%x_dim = size(mri_inst%x_coordinates)
  CALL mr_io_read_coordinates(grp_id, "y", mri_inst%y_coordinates)
  mri_inst%y_dim = size(mri_inst%y_coordinates)
  CALL mr_io_read_coordinates(grp_id, "z", mri_inst%z_coordinates)
  mri_inst%z_dim = size(mri_inst%z_coordinates)

  ! Read voxel_feature data
  CALL mr_io_read_spacetime_scalar_feature(grp_id, "intensity", mri_inst%intensity)
  mri_inst%intensity_dims =  shape(mri_inst%intensity)
  CALL mr_io_read_spacetime_feature(grp_id, "velocity_mean", mri_inst%velocity_mean)
  mri_inst%velocity_mean_dims =  shape(mri_inst%velocity_mean)
  CALL mr_io_read_spacetime_matrix_feature(grp_id, "velocity_cov", mri_inst%velocity_cov)
  mri_inst%velocity_cov_dims =  shape(mri_inst%velocity_cov)
  CALL mr_io_read_spacetime_scalar_feature(grp_id, "segmentation_prob", mri_inst%segmentation_prob)
  mri_inst%segmentation_prob_dims =  shape(mri_inst%segmentation_prob)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_read_segmentedhpcpredict



subroutine mr_io_write_segmentedhpcpredict(path, mri_inst)

  implicit none

  character(len=*), intent(in) :: path
  type(SegmentedHPCPredictMRI), intent(in) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier

  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  mr_io_handle_error(error)

  ! Create a new file using default properties.
  CALL h5fcreate_f(trim(path), H5F_ACC_TRUNC_F, file_id, error)
  mr_io_handle_error(error)

  ! Create a new group
  CALL h5gcreate_f(file_id, SegmentedHPCPredictMRI_group_name, grp_id, error)
  mr_io_handle_error(error)

  ! Write time data
  CALL mr_io_write_coordinates(grp_id, "t", mri_inst%t_coordinates)

  ! Write geometry data
  CALL mr_io_write_coordinates(grp_id, "x", mri_inst%x_coordinates)
  CALL mr_io_write_coordinates(grp_id, "y", mri_inst%y_coordinates)
  CALL mr_io_write_coordinates(grp_id, "z", mri_inst%z_coordinates)

  ! Write voxel_feature data
  CALL mr_io_write_spacetime_scalar_feature(grp_id, "intensity", mri_inst%intensity)
  CALL mr_io_write_spacetime_feature(grp_id, "velocity_mean", mri_inst%velocity_mean)
  CALL mr_io_write_spacetime_matrix_feature(grp_id, "velocity_cov", mri_inst%velocity_cov)
  CALL mr_io_write_spacetime_scalar_feature(grp_id, "segmentation_prob", mri_inst%segmentation_prob)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_write_segmentedhpcpredict


end module mr_io
