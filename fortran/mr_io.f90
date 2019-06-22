module mr_io

use hdf5
use mr_protocol


include 'mpif.h'

contains


#ifdef __GFORTRAN__

#define mr_io_handle_error(err) if(err /= 0) CALL mr_io_handle_hdf5_error(err);
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

#endif


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

! TODO: Parallel version

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



subroutine mr_io_read_parallel_spatial_feature(mpi_comm, grp_id, feature_name, feature_array)

  implicit none

  INTEGER, intent(in) :: mpi_comm
  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name
  real*8, dimension(:,:,:), allocatable, intent(out) :: feature_array

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier   
  
  INTEGER(HSIZE_T), DIMENSION(3) :: dims_file = (/-1, -1, -1/) ! Shape of full dataset
  INTEGER(HSIZE_T), DIMENSION(3) :: max_dims_file = (/-1, -1, -1/) ! Max shape of full dataset
  
  INTEGER(HSIZE_T), DIMENSION(3) :: offset_file = (/-1, -1, -1/) ! Offset of data subset to read/write
  INTEGER(HSIZE_T), DIMENSION(3) :: dims_mem = (/-1, -1, -1/) ! Shape of data subset to read/write

  INTEGER(HSIZE_T), DIMENSION(3) :: count_blocks = (/-1, -1, -1/) ! number of blocks to read/write
  INTEGER(HSIZE_T), DIMENSION(3) :: stride = (/-1, -1, -1/) ! stride between adjacent elements/blocks?? (unused when single block is read)

  INTEGER     ::   rank = 3       ! Dataset rank
  INTEGER     ::   error          ! Error flag

  INTEGER     ::   mpi_rank, mpi_size, mpi_error
  
  INTEGER     ::   i
  
  ! Open an existing dataset.
  CALL h5dopen_f(grp_id, feature_name, dset_id, error)
  mr_io_handle_error(error)
  
  ! Get file space
  CALL h5dget_space_f(dset_id, filespace, error)
  mr_io_handle_error(error)
  CALL h5sget_simple_extent_dims_f(filespace, dims_file, max_dims_file, error)
  mr_io_sget_simple_extent_dims_handle_error(error)

  ! Compute hyperslab offset and shape
  CALL MPI_Comm_rank(mpi_comm, mpi_rank, mpi_error)
  mr_io_handle_error(mpi_error) ! FIXME: MPI error handling
  CALL MPI_Comm_size(mpi_comm, mpi_size, mpi_error)
  mr_io_handle_error(mpi_error)
  if (mpi_rank + 1 /= mpi_size) then
      dims_mem(3) = (dims_file(3) + mpi_size - 1)/mpi_size
  else
      dims_mem(3) =  modulo(dims_file(3),mpi_size)
      if (dims_mem(3) == 0) then
          dims_mem(3) = dims_file(3)/mpi_size
      endif
  endif
  dims_mem(1:2) = dims_file(1:2)
  
  offset_file(3) = (dims_file(3) + mpi_size - 1)/mpi_size*mpi_rank
  offset_file(1:2) = 0
  
  count_blocks = (/1,1,1/)
  stride = (/1,1,1/)
  
  !write (*,*) "offset_file"
  !do i=1,3
  !  write (*,*) offset_file(i)
  !end do
  !write (*,*) "dims_mem"
  !do i=1,3
  !  write (*,*) dims_mem(i)
  !end do
  !write (*,*) "dims_file"
  !do i=1,3
  !  write (*,*) dims_file(i)
  !end do
  
  ! Create memory space
  CALL h5screate_simple_f(rank, dims_mem, memspace, error)
  mr_io_handle_error(error)

  ! Select hyperslab in the file.
  CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_file, count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)

  ! TODO: Select hyperslab in the memspace if it feature_array doesn't agree with layout to read
  
  ! Allocate MRI array 
  allocate(feature_array(dims_mem(1), dims_mem(2), dims_mem(3)))
  
  ! Create property list for collective dataset read
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  mr_io_handle_error(error)
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  mr_io_handle_error(error)

  ! Read the dataset collectively:
  CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, feature_array, dims_file, &
                 error, mem_space_id=memspace, file_space_id=filespace, &
                 xfer_prp=plist_id)
  mr_io_handle_error(error)

  ! Close the file space, memory space and property list
  CALL h5sclose_f(filespace, error)
  mr_io_handle_error(error)
  CALL h5sclose_f(memspace , error)
  mr_io_handle_error(error)
  CALL h5pclose_f(plist_id , error)
  mr_io_handle_error(error)
  
!   ! Read dims by getting dimension from data space
!   CALL h5dget_space_f(dset_id, dspace_id, error)
!   mr_io_handle_error(error)
!   CALL h5sget_simple_extent_dims_f(dspace_id, dims, max_dims, error)
!   mr_io_sget_simple_extent_dims_handle_error(error)
! 
!   ! Allocate MRI array 
!   allocate(feature_array(dims(1), dims(2), dims(3)))
! 
!   ! Read the dataset.
!   CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, feature_array, dims, error)
!   mr_io_handle_error(error)
! 
!   ! Close the data space
!   CALL h5sclose_f(dspace_id, error)
!   mr_io_handle_error(error)
  
  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_read_parallel_spatial_feature



subroutine mr_io_read_parallel_spatial(mpi_comm, mpi_info, path, mri_inst)     

  implicit none

  character(len=*), intent(in) :: path
  INTEGER, intent (in) :: mpi_comm, mpi_info
  type(SpatialMRI), intent(out) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier

  INTEGER(HID_T) :: plist_id      ! Property list identifier   
  
  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  mr_io_handle_error(error)

  ! Setup file access property list with parallel I/O access.
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  mr_io_handle_error(error)
  
  CALL h5pset_fapl_mpio_f(plist_id, mpi_comm, mpi_info, error)
  mr_io_handle_error(error)

  ! Open existing file collectively 
  CALL h5fopen_f(trim(path), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
  mr_io_handle_error(error)

  CALL h5pclose_f(plist_id, error)
  mr_io_handle_error(error)
  
  ! Open an existing group
  CALL h5gopen_f(file_id, SpatialMRI_group_name, grp_id, error)
  mr_io_handle_error(error)

  ! Read spatial feature
  CALL mr_io_read_parallel_spatial_feature(mpi_comm, grp_id, "voxel_feature", mri_inst%voxel_feature)
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

end subroutine mr_io_read_parallel_spatial






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
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HSIZE_T), DIMENSION(1) :: coord_dim = (/-1/) ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(1) :: max_coord_dim = (/-1/) ! Dataset dimensions

  INTEGER(HSIZE_T), DIMENSION(5) :: dims = (/-1, -1, -1, -1, -1/) ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(5) :: max_dims = (/-1, -1, -1, -1, -1/) ! Dataset dimensions

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
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
 
  INTEGER(HSIZE_T), DIMENSION(1) :: coord_dims            ! Coordinate dataset dimensions

  INTEGER(HSIZE_T), DIMENSION(5) :: dims            ! Dataset dimensions
  INTEGER     ::   rank = 5                         ! Dataset rank

  INTEGER     ::   error ! Error flag

  dims = mri_inst%voxel_feature_dims

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  mr_io_handle_error(error)

  ! Create a new file using default properties.
  CALL h5fcreate_f(trim(path), H5F_ACC_TRUNC_F, file_id, error)
  mr_io_handle_error(error)

  ! Create a new group
  CALL h5gcreate_f(file_id, SpaceTimeMRI_group_name, grp_id, error)
  mr_io_handle_error(error)

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


end module mr_io
