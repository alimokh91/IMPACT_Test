module mr_io_parallel_spacetime

use hdf5
use mr_io_protocol
use mr_io, only : mr_io_handle_hdf5_error, &
                  mr_io_handle_argument_error
use mr_io_parallel, only : mr_io_parallel_spatial_hyperslap_compute, &
                           mr_io_parallel_spatial_hyperslab_get_type, &
                           mr_io_parallel_spatial_hyperslap_compute, &
                           mr_io_parallel_spatial_hyperslap_compute_padded

use mpi !include 'mpif.h'

implicit none

private

public :: mr_io_read_parallel_spacetime
public :: mr_io_write_parallel_spacetime
public :: mr_io_read_parallel_flow
public :: mr_io_write_parallel_flow
public :: mr_io_read_parallel_flow_padded
public :: mr_io_read_parallel_segmentedflow
public :: mr_io_write_parallel_segmentedflow

public :: DistSpaceTimeMRI
public :: DistFlowMRI
public :: DistFlowMRIPadded
public :: DistSegmentedFlowMRI

public :: SpaceTimeMRI_group_name
public :: FlowMRI_group_name
public :: SegmentedFlowMRI_group_name

public :: mr_io_deallocate_dist_spacetime_mri
public :: mr_io_deallocate_dist_flow_mri
public :: mr_io_deallocate_dist_flow_mri_padded
public :: mr_io_deallocate_dist_segmentedflow_mri


#ifdef __GFORTRAN__
! Spatial hyperslab macro enum
! Compute spatial hyperslab
#define mr_io_parallel_spatial_hyperslab_enum_compute   ( 0)
! Read specified hyperslab
#define mr_io_parallel_spatial_hyperslab_enum_specified ( 1)
! Inconsistent arguments (shape/offset/feature_array)
#define mr_io_parallel_spatial_hyperslab_enum_error     (-1)
#endif

contains

#ifdef __GFORTRAN__

#define mr_io_handle_error(err) if(err /= 0) CALL mr_io_handle_hdf5_error(err);
#define mr_io_handle_arg_error(err) if(err /= 0) CALL mr_io_handle_argument_error(err);
#define mr_io_sget_simple_extent_dims_handle_error(err) if(err == -1) CALL mr_io_handle_hdf5_error(err);

#define mr_io_handle_mpi_error(err) if(err /= 0) CALL mr_io_handle_mpi_error_(err);

subroutine mr_io_handle_mpi_error_(error)

  implicit none

  INTEGER     ::   error          ! Error flag

  if(error /= 0) then
    write(*,*) "MPI error ",error,"- printing backtrace and aborting..."
    call backtrace
    !call abort
  end if

end subroutine mr_io_handle_mpi_error_

#endif



! ************************ DistSpaceTimeMRI ************************

subroutine mr_io_read_parallel_spacetime_feature(mr_io_mpi_comm, mr_io_mpi_cart_dims, grp_id, feature_name, &
                                                 feature_array, feature_offset, feature_shape, &
                                                 time_offset, time_dim) !, feature_halo_shape)

  implicit none

  INTEGER, intent(in) :: mr_io_mpi_comm
  integer, dimension(3), intent(in) :: mr_io_mpi_cart_dims
  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name

  real*8, dimension(:,:,:,:,:), allocatable, intent(inout) :: feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(inout) :: feature_shape  ! Shape of data subset to read/write in feature_array
  !INTEGER, DIMENSION(3), OPTIONAL, intent(in)    :: feature_halo_shape ! Shape of halo around data to read/write in feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(inout) :: feature_offset ! Offset of data subset to read/write in feature_array

  INTEGER, intent(inout) :: time_offset
  INTEGER, intent(inout) :: time_dim

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier

  INTEGER(HSIZE_T), DIMENSION(5) :: dims_file = (/-1, -1, -1, -1, -1/) ! Shape of full HDF5 dataset
  INTEGER(HSIZE_T), DIMENSION(5) :: max_dims_file = (/-1, -1, -1, -1, -1/) ! Max shape of full dataset

  INTEGER(HSIZE_T), DIMENSION(5) :: offset_file = (/-1, -1, -1, -1, -1/) ! Offset of data subset to read/write in HDF5 dataset

  INTEGER(HSIZE_T), DIMENSION(5) :: dims_mem = (/-1, -1, -1, -1, -1/) ! Shape of data subset to read/write in feature_array
  INTEGER(HSIZE_T), DIMENSION(5) :: offset_mem = (/-1, -1, -1, -1, -1/) ! Offset of data subset to read/write in feature_array

  INTEGER(HSIZE_T), DIMENSION(5) :: count_blocks = (/-1, -1, -1, -1, -1/) ! number of blocks to read/write
  INTEGER(HSIZE_T), DIMENSION(5) :: stride = (/-1, -1, -1, -1, -1/) ! stride between adjacent elements/blocks?? (unused when single block is read)

  INTEGER     ::   rank = 5       ! Dataset rank
  INTEGER     ::   error          ! Error flag

  INTEGER     ::   i

  !INTEGER     ::   spatial_hyperslab_type = mr_io_parallel_spatial_hyperslab_enum_error
  INTEGER, DIMENSION(5) :: feature_array_shape_tmp = (/-1,-1,-1,-1,-1/)
  INTEGER, DIMENSION(3) :: feature_array_shape = (/-1,-1,-1/)

  ! Open an existing dataset.
  CALL h5dopen_f(grp_id, feature_name, dset_id, error)
  mr_io_handle_error(error)

  ! Get file space
  CALL h5dget_space_f(dset_id, filespace, error)
  mr_io_handle_error(error)
  CALL h5sget_simple_extent_dims_f(filespace, dims_file, max_dims_file, error)
  mr_io_sget_simple_extent_dims_handle_error(error)


  if ( allocated(feature_array) ) then
    feature_array_shape_tmp = shape(feature_array)
    feature_array_shape = feature_array_shape_tmp(3:5)
  else ! feature_array not allocated - compute domain decomposition
    feature_array_shape = (/-1,-1,-1/)
  end if


  ! Includes argument check
  select case (mr_io_parallel_spatial_hyperslab_get_type(feature_shape, feature_offset, feature_array_shape))

  case ( mr_io_parallel_spatial_hyperslab_enum_specified )
    dims_mem = shape(feature_array)
    ! FIXME: could be imporved
    offset_file(1) = 0                   ! vectorial
    offset_file(2) = time_offset         ! temporal
    offset_file(3:5) = feature_offset    ! spatial

    if( any(dims_file(3:5) /= feature_shape) .or. (dims_file(2) /= time_dim) ) then
      write(*,*) "Array dimensions read from file and supplied in feature_shape/time_shape disagree."
      mr_io_handle_arg_error(-1)
    endif

  case ( mr_io_parallel_spatial_hyperslab_enum_compute )
    ! vectorial_dim skipped
    time_dim = dims_file(2)
    feature_shape = dims_file(3:5)

    ! vectorial
    offset_file(1) = 0
    dims_mem(1) = dims_file(1)

    ! temporal
    if (time_offset < 0) then
      time_offset = 0
    endif
    offset_file(2) = time_offset

    dims_mem(2) = dims_file(2)-offset_file(2)

    ! spatial
    call mr_io_parallel_spatial_hyperslap_compute(mr_io_mpi_comm, mr_io_mpi_cart_dims, feature_shape, feature_offset, dims_mem(3:5))
    ! feature_offset set in mr_io_parallel_spatial_hyperslap_compute
    offset_file(3:5) = feature_offset

    ! MRI array - TODO: handle feature_halo_shape
    allocate(feature_array(dims_mem(1), dims_mem(2), dims_mem(3), dims_mem(4), dims_mem(5)))

  case default
    mr_io_handle_arg_error(-1)

  end select

  ! TODO: handle feature_halo_shape
  offset_mem = (/ 0, 0, 0, 0 ,0/)

  count_blocks = (/ 1, 1, 1, 1, 1/)
  stride = (/ 1, 1, 1, 1, 1/)

!  write (*,*) "dims_file"
!  do i=1,5
!    write (*,*) dims_file(i)
!  end do
!  write (*,*) "offset_file"
!  do i=1,5
!    write (*,*) offset_file(i)
!  end do
!  write (*,*) "dims_mem"
!  do i=1,5
!    write (*,*) dims_mem(i)
!  end do

  ! Create memory space
  CALL h5screate_simple_f(rank, dims_mem, memspace, error)
  mr_io_handle_error(error)

  ! Select hyperslab in the file.
  CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_file, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)
  CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset_mem, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)

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

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_read_parallel_spacetime_feature


subroutine mr_io_read_parallel_spacetime_feature_padded(mr_io_mpi_comm, mr_io_mpi_cart_dims, grp_id, &
                                                 domain_padding, &
                                                 feature_name, &
                                                 feature_array, feature_offset, feature_shape, &
                                                 time_offset, time_dim) !, feature_halo_shape)

  implicit none

  INTEGER, intent(in) :: mr_io_mpi_comm
  integer, dimension(3), intent(in) :: mr_io_mpi_cart_dims
  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  type(DomainPadding) :: domain_padding
  character(len=*), intent(in) :: feature_name

  real*8, dimension(:,:,:,:,:), allocatable, intent(inout) :: feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(inout) :: feature_shape  ! Shape of data subset to read/write in feature_array
  !INTEGER, DIMENSION(3), OPTIONAL, intent(in)    :: feature_halo_shape ! Shape of halo around data to read/write in feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(inout) :: feature_offset ! Offset of data subset to read/write in feature_array

  INTEGER, intent(inout) :: time_offset
  INTEGER, intent(inout) :: time_dim

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier

  INTEGER(HSIZE_T), DIMENSION(5) :: dims_file = (/-1, -1, -1, -1, -1/) ! Shape of full HDF5 dataset
  INTEGER(HSIZE_T), DIMENSION(5) :: max_dims_file = (/-1, -1, -1, -1, -1/) ! Max shape of full dataset

  INTEGER(HSIZE_T), DIMENSION(5) :: offset_file = (/-1, -1, -1, -1, -1/) ! Offset of data subset to read/write in HDF5 dataset

  INTEGER(HSIZE_T), DIMENSION(5) :: dims_mem = (/-1, -1, -1, -1, -1/) ! Shape of data subset to read/write in feature_array
  INTEGER(HSIZE_T), DIMENSION(5) :: offset_mem = (/-1, -1, -1, -1, -1/) ! Offset of data subset to read/write in feature_array
  INTEGER(HSIZE_T), DIMENSION(3) :: offset_local_hyperslab = (/0, 0, 0/) ! Offset of local hyperslab

  INTEGER(HSIZE_T), DIMENSION(5) :: count_blocks = (/-1, -1, -1, -1, -1/) ! number of blocks to read/write
  INTEGER(HSIZE_T), DIMENSION(5) :: stride = (/-1, -1, -1, -1, -1/) ! stride between adjacent elements/blocks?? (unused when single block is read)

  INTEGER     ::   rank = 5       ! Dataset rank
  INTEGER     ::   error          ! Error flag

  INTEGER     ::   i

  !INTEGER     ::   spatial_hyperslab_type = mr_io_parallel_spatial_hyperslab_enum_error
  INTEGER, DIMENSION(5) :: feature_array_shape_tmp = (/-1,-1,-1,-1,-1/)
  INTEGER, DIMENSION(3) :: feature_array_shape = (/-1,-1,-1/)

  ! Open an existing dataset.
  CALL h5dopen_f(grp_id, feature_name, dset_id, error)
  mr_io_handle_error(error)

  ! Get file space
  CALL h5dget_space_f(dset_id, filespace, error)
  mr_io_handle_error(error)
  CALL h5sget_simple_extent_dims_f(filespace, dims_file, max_dims_file, error)
  mr_io_sget_simple_extent_dims_handle_error(error)


  if ( allocated(feature_array) ) then
    feature_array_shape_tmp = shape(feature_array)
    feature_array_shape = feature_array_shape_tmp(3:5)
  else ! feature_array not allocated - compute domain decomposition
    feature_array_shape = (/-1,-1,-1/)
  end if


  ! Includes argument check
  select case (mr_io_parallel_spatial_hyperslab_get_type(feature_shape, feature_offset, feature_array_shape))

  case ( mr_io_parallel_spatial_hyperslab_enum_specified )
    dims_mem = shape(feature_array)
    ! FIXME: could be imporved
    offset_file(1) = 0                   ! vectorial
    offset_file(2) = time_offset         ! temporal
    offset_file(3:5) = feature_offset    ! spatial

    if( any(dims_file(3:5) /= feature_shape) .or. (dims_file(2) /= time_dim) ) then
      write(*,*) "Array dimensions read from file and supplied in feature_shape/time_shape disagree."
      mr_io_handle_arg_error(-1)
    endif

  case ( mr_io_parallel_spatial_hyperslab_enum_compute )
    ! vectorial_dim skipped
    time_dim = dims_file(2)
    feature_shape = dims_file(3:5)

    ! vectorial
    offset_file(1) = 0
    dims_mem(1) = dims_file(1)

    ! temporal
    if (time_offset < 0) then
      time_offset = 0
    endif
    offset_file(2) = time_offset

    dims_mem(2) = dims_file(2)-offset_file(2)

    ! spatial
    call mr_io_parallel_spatial_hyperslap_compute_padded(mr_io_mpi_comm, mr_io_mpi_cart_dims, domain_padding, &
                                                         feature_shape, feature_offset, dims_mem(3:5), offset_local_hyperslab)
    ! feature_offset set in mr_io_parallel_spatial_hyperslap_compute
    offset_file(3:5) = feature_offset

    ! MRI array - TODO: handle feature_halo_shape
    allocate(feature_array(dims_mem(1), dims_mem(2), &
                           offset_local_hyperslab(1)+1:offset_local_hyperslab(1)+dims_mem(3), &
                           offset_local_hyperslab(2)+1:offset_local_hyperslab(2)+dims_mem(4), &
                           offset_local_hyperslab(3)+1:offset_local_hyperslab(3)+dims_mem(5)))
!    write(0,*) "shape of local feature array: ",shape(feature_array)

  case default
    mr_io_handle_arg_error(-1)

  end select

  ! TODO: handle feature_halo_shape
  offset_mem = (/ 0, 0, 0, 0 ,0/)

  count_blocks = (/ 1, 1, 1, 1, 1/)
  stride = (/ 1, 1, 1, 1, 1/)

!  write (*,*) "dims_file"
!  do i=1,5
!    write (*,*) dims_file(i)
!  end do
!  write (*,*) "offset_file"
!  do i=1,5
!    write (*,*) offset_file(i)
!  end do
!  write (*,*) "dims_mem"
!  do i=1,5
!    write (*,*) dims_mem(i)
!  end do

  ! Create memory space
  CALL h5screate_simple_f(rank, dims_mem, memspace, error)
  mr_io_handle_error(error)

  ! Select hyperslab in the file.
  CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_file, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)
  CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset_mem, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)

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

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_read_parallel_spacetime_feature_padded


subroutine mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, coordinate, coordinate_array)

  implicit none

  INTEGER, intent(in) :: mr_io_mpi_comm
  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: coordinate
  real*8, dimension(:), allocatable, intent(out) :: coordinate_array

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier in file

  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier

  INTEGER(HSIZE_T), DIMENSION(1) :: coord_dim = (/-1/) ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(1) :: max_coord_dim = (/-1/) ! Dataset dimensions

  INTEGER     ::   error          ! Error flag

  ! Open an existing dataset.
  CALL h5dopen_f(grp_id, coordinate//"_coordinates", dset_id, error)
  mr_io_handle_error(error)

  ! Get file space
  CALL h5dget_space_f(dset_id, dspace_id, error)
  mr_io_handle_error(error)
  CALL h5sget_simple_extent_dims_f(dspace_id, coord_dim, max_coord_dim, error)
  mr_io_sget_simple_extent_dims_handle_error(error)

  ! Allocate coordinate array
  allocate(coordinate_array(coord_dim(1)))

  ! Create property list for collective dataset read
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  mr_io_handle_error(error)
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  mr_io_handle_error(error)

  ! Read the dataset collectively:
  CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, coordinate_array, coord_dim, error, xfer_prp=plist_id)
  mr_io_handle_error(error)

  ! Close the file space, memory space and property list
  CALL h5sclose_f(dspace_id, error)
  mr_io_handle_error(error)
  CALL h5pclose_f(plist_id , error)
  mr_io_handle_error(error)

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_read_parallel_coordinates


subroutine mr_io_read_parallel_spacetime(mr_io_mpi_comm, mr_io_mpi_info, mr_io_mpi_cart_dims, path, mri_inst)

  implicit none

  character(len=*), intent(in) :: path
  INTEGER, intent (in) :: mr_io_mpi_comm, mr_io_mpi_info
  integer, dimension(3), intent(in) :: mr_io_mpi_cart_dims
  type(DistSpacetimeMRI), intent(out) :: mri_inst

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

  CALL h5pset_fapl_mpio_f(plist_id, mr_io_mpi_comm, mr_io_mpi_info, error)
  mr_io_handle_error(error)

  ! Open existing file collectively
  CALL h5fopen_f(trim(path), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
  mr_io_handle_error(error)

  CALL h5pclose_f(plist_id, error)
  mr_io_handle_error(error)

  ! Open an existing group
  CALL h5gopen_f(file_id, SpaceTimeMRI_group_name, grp_id, error)
  mr_io_handle_error(error)

  ! Read time coordinates
  CALL mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, "t", &
                                       mri_inst%t_coordinates)

  ! Read coordinates
  CALL mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, "x", &
                                       mri_inst%x_coordinates)
  CALL mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, "y", &
                                       mri_inst%y_coordinates)
  CALL mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, "z", &
                                       mri_inst%z_coordinates)


  ! Read spatial feature
  CALL mr_io_read_parallel_spacetime_feature(mr_io_mpi_comm, &
                                             mr_io_mpi_cart_dims, &
                                             grp_id, "vector_feature", &
                                             mri_inst%vector_feature%array, &
                                             mri_inst%vector_feature%offset, &
                                             mri_inst%vector_feature%dims, &
                                             mri_inst%vector_feature%time_offset, &
                                             mri_inst%vector_feature%time_dim)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_read_parallel_spacetime


subroutine mr_io_write_parallel_spacetime_feature(mr_io_mpi_comm, &
                                                  grp_id, feature_name, &
                                                  feature_array, feature_offset, feature_shape, &
                                                  time_offset, time_dim) !, feature_halo_shape)

  implicit none

  INTEGER, intent(in) :: mr_io_mpi_comm
  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name

  real*8, dimension(:,:,:,:,:), allocatable, intent(in) :: feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(in) :: feature_shape  ! Shape of data subset to read/write in feature_array
  !INTEGER, DIMENSION(3), OPTIONAL, intent(in)    :: feature_halo_shape ! Shape of halo around data to read/write in feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(in) :: feature_offset ! Offset of data subset to read/write in feature_array

  INTEGER, intent(in) :: time_offset
  INTEGER, intent(in) :: time_dim

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier

  INTEGER(HSIZE_T), DIMENSION(5) :: dims_file = (/-1, -1, -1, -1, -1/) ! Shape of full HDF5 dataset
  INTEGER(HSIZE_T), DIMENSION(5) :: offset_file = (/-1, -1, -1, -1, -1/) ! Offset of data subset to read/write in HDF5 dataset

  INTEGER(HSIZE_T), DIMENSION(5) :: dims_mem = (/-1, -1, -1, -1, -1/) ! Shape of data subset to read/write in feature_array
  INTEGER(HSIZE_T), DIMENSION(5) :: offset_mem = (/-1, -1, -1, -1, -1/) ! Offset of data subset to read/write in feature_array

  INTEGER(HSIZE_T), DIMENSION(5) :: count_blocks = (/-1, -1, -1, -1, -1/) ! number of blocks to read/write
  INTEGER(HSIZE_T), DIMENSION(5) :: stride = (/-1, -1, -1, -1, -1/) ! stride between adjacent elements/blocks?? (unused when single block is read)

  INTEGER     ::   rank = 5       ! Dataset rank
  INTEGER     ::   error          ! Error flag

  INTEGER     ::   i

  INTEGER, DIMENSION(5) :: feature_array_shape

  integer, dimension(5) :: arr

  feature_array_shape = shape(feature_array)

  ! Check arguments
  if ( any(feature_offset < (/0, 0, 0/)) .or. &
       ( any(feature_shape < feature_offset + feature_array_shape(3:5)) .and. any(feature_array_shape(3:5) /= (/0,0,0/)) ) .or. &
       ( time_offset < 0 ) .or. &
       ( time_dim < time_offset + feature_array_shape(2) ) ) then
    write(0,*) "spatial_global_shape: ",feature_shape
    write(0,*) "spatial_offset:       ",feature_offset
    write(0,*) "spatial_local_shape:  ",feature_array_shape(3:5)
    write(0,*) "time_global_shape:    ",time_dim
    write(0,*) "time_offset:          ",time_offset
    write(0,*) "time_local_shape:     ",feature_array_shape(2)
    write(0,*) "Invalid arguments (hyperslab must provide sufficient space for value array)"
    call flush()
    mr_io_handle_arg_error(-1)
  endif

  ! TODO: handle if( present(feature_halo_shape) )
  dims_mem = shape(feature_array)
  offset_mem = (/0,0,0,0,0/)

  ! Compute file dimensions and offset
  dims_file(1) = dims_mem(1)
  dims_file(2) = time_dim
  dims_file(3:5) = feature_shape

  offset_file(1) = 0
  offset_file(2) = time_offset
  offset_file(3:5) = feature_offset


!   write (*,*) "dims_file"
!   do i=1,3
!    write (*,*) dims_file(i)
!   end do
!   write (*,*) "feature_shape"
!   do i=1,3
!    write (*,*) feature_shape(i)
!   end do
!   write (*,*) "feature_offset"
!   do i=1,3
!    write (*,*) feature_offset(i)
!   end do
!   write (*,*) "offset_file"
!   do i=1,3
!    write (*,*) offset_file(i)
!   end do
!   write (*,*) "dims_mem"
!   do i=1,3
!    write (*,*) dims_mem(i)
!   end do
!   write (*,*) "shape(feature_array)"
!   arr = shape(feature_array)
!   do i=1,3
!     write(*,*) arr(i)
!   end do
!   write (*,*) "offset_mem"
!   do i=1,3
!    write (*,*) offset_mem(i)
!   end do
!

  ! Create the data space for the  dataset.   !
  CALL h5screate_simple_f(rank, dims_file, filespace, error)
  mr_io_handle_error(error)
  CALL h5screate_simple_f(rank, dims_mem, memspace, error)
  mr_io_handle_error(error)

  ! Create the dataset with default properties.
  CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  mr_io_handle_error(error)
  CALL h5dcreate_f(grp_id, feature_name, H5T_NATIVE_DOUBLE, filespace, &
       dset_id, error, dcpl_id=plist_id)
  mr_io_handle_error(error)
!   CALL h5pclose_f(plist_id , error)
!   mr_io_handle_error(error)

  count_blocks = (/1,1,1,1,1/)
  stride = (/1,1,1,1,1/)

  ! Select hyperslab in the file.
  CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_file, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)
  CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset_mem, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)

  ! Create property list for collective dataset read
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  mr_io_handle_error(error)
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  mr_io_handle_error(error)

  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, feature_array, dims_file, &
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

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_write_parallel_spacetime_feature


subroutine mr_io_write_parallel_coordinates(mr_io_mpi_comm, grp_id, coordinate, coordinate_array) !, feature_halo_shape)

  implicit none

  INTEGER, intent(in) :: mr_io_mpi_comm
  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: coordinate
  real*8, dimension(:), allocatable, intent(in) :: coordinate_array

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
  INTEGER(HID_T) :: plist_id      ! Property list identifier

  INTEGER(HSIZE_T), DIMENSION(1) :: coord_dim = (/-1/) ! Dataset dimensions

  INTEGER(HSIZE_T), DIMENSION(1) :: dims_file = (/-1/) ! Shape of full HDF5 dataset
  INTEGER(HSIZE_T), DIMENSION(1) :: offset_file = (/-1/) ! Offset of data subset to read/write in HDF5 dataset

  INTEGER(HSIZE_T), DIMENSION(1) :: dims_mem = (/-1/) ! Shape of data subset to read/write in feature_array
  INTEGER(HSIZE_T), DIMENSION(1) :: offset_mem = (/-1/) ! Offset of data subset to read/write in feature_array

  INTEGER(HSIZE_T), DIMENSION(1) :: count_blocks = (/-1/) ! number of blocks to read/write
  INTEGER(HSIZE_T), DIMENSION(1) :: stride = (/-1/) ! stride between adjacent elements/blocks?? (unused when single block is read)

  INTEGER     ::   rank = 1       ! Dataset rank
  INTEGER     ::   error          ! Error flag

  INTEGER     ::   mr_io_mpi_rank, mr_io_mpi_size, mr_io_mpi_error

  INTEGER     ::   i

  ! Compute hyperslab offset and shape
  CALL mpi_comm_rank(mr_io_mpi_comm, mr_io_mpi_rank, mr_io_mpi_error)
  mr_io_handle_error(mr_io_mpi_error) ! FIXME: MPI error handling
  CALL mpi_comm_size(mr_io_mpi_comm, mr_io_mpi_size, mr_io_mpi_error)
  mr_io_handle_error(mr_io_mpi_error)


  dims_file = shape(coordinate_array)
  if (mr_io_mpi_rank == 0) then
    dims_mem = dims_file
  else
    dims_mem = (/0/)
  end if

  offset_file = (/0/)
  offset_mem = (/0/)

!  dims_file = shape(coordinate_array)
!  if(mr_io_mpi_rank + 1 == mr_io_mpi_size) then
!    if( modulo(dims_file(1), mr_io_mpi_size) == 0 ) then
!      dims_mem(1) = dims_file(1) / mr_io_mpi_size
!    else
!      dims_mem(1) = modulo(dims_file(1), mr_io_mpi_size)
!    end if
!  else
!    dims_mem(1) = (dims_file(1) + mr_io_mpi_size - 1) / mr_io_mpi_size
!  end if
!  offset_file(1) = (dims_file(1) + mr_io_mpi_size - 1) / mr_io_mpi_size * mr_io_mpi_rank
!  offset_mem = offset_file

!   write (*,*) "dims_file"
!   do i=1,3
!    write (*,*) dims_file(i)
!   end do
!   write (*,*) "feature_shape"
!   do i=1,3
!    write (*,*) feature_shape(i)
!   end do
!   write (*,*) "feature_offset"
!   do i=1,3
!    write (*,*) feature_offset(i)
!   end do
!   write (*,*) "offset_file"
!   do i=1,3
!    write (*,*) offset_file(i)
!   end do
!   write (*,*) "dims_mem"
!   do i=1,3
!    write (*,*) dims_mem(i)
!   end do
!   write (*,*) "shape(feature_array)"
!   arr = shape(feature_array)
!   do i=1,3
!     write(*,*) arr(i)
!   end do
!   write (*,*) "offset_mem"
!   do i=1,3
!    write (*,*) offset_mem(i)
!   end do
!

  ! Create the data space for the  dataset.   !
  CALL h5screate_simple_f(rank, dims_file, filespace, error)
  mr_io_handle_error(error)

  ! Create the dataset with default properties.
  CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  mr_io_handle_error(error)
  CALL h5dcreate_f(grp_id,  coordinate//"_coordinates", H5T_NATIVE_DOUBLE, filespace, &
       dset_id, error, dcpl_id=plist_id)
  mr_io_handle_error(error)
!   CALL h5pclose_f(plist_id , error)
!   mr_io_handle_error(error)

  if (mr_io_mpi_rank == 0) then
    count_blocks = (/1/)
  else
    count_blocks = (/1/)
  end if
  stride = (/1/)

  ! Select hyperslab in the file.
  CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_file, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)

  ! Create property list for collective dataset read
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  mr_io_handle_error(error)
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  mr_io_handle_error(error)

  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, coordinate_array, dims_file, &
                  error, file_space_id=filespace, &
                  xfer_prp=plist_id)
  mr_io_handle_error(error)

  ! Close the file space, memory space and property list
  CALL h5sclose_f(filespace, error)
  mr_io_handle_error(error)
  CALL h5pclose_f(plist_id , error)
  mr_io_handle_error(error)

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_write_parallel_coordinates




subroutine mr_io_write_parallel_spacetime(mr_io_mpi_comm, mr_io_mpi_info, path, mri_inst)

  implicit none

  character(len=*), intent(in) :: path
  INTEGER, intent (in) :: mr_io_mpi_comm, mr_io_mpi_info
  type(DistSpacetimeMRI), intent(in) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier

  INTEGER(HID_T) :: plist_id      ! Property list identifier

  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  mr_io_handle_error(error)

  ! Create a file collectively
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  mr_io_handle_error(error)

  CALL h5pset_fapl_mpio_f(plist_id, mr_io_mpi_comm, mr_io_mpi_info, error)
  mr_io_handle_error(error)

  CALL h5fcreate_f(trim(path), H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
  mr_io_handle_error(error)

  CALL h5pclose_f(plist_id, error)
  mr_io_handle_error(error)

  ! Create a new group collectively
  CALL h5pcreate_f(H5P_GROUP_CREATE_F, plist_id, error) ! FIXME: Not completely sure about flag
  mr_io_handle_error(error)

  CALL h5gcreate_f(file_id, SpaceTimeMRI_group_name, grp_id, error, gcpl_id=plist_id)
  mr_io_handle_error(error)

  CALL h5pclose_f(plist_id, error)
  mr_io_handle_error(error)

!   print *, mri_inst%vector_feature%dims
!   print *, mri_inst%vector_feature%offset
!   print *, shape(mri_inst%vector_feature%array)
  CALL mr_io_write_parallel_coordinates(mr_io_mpi_comm, grp_id, "t", &
                                        mri_inst%t_coordinates)
  CALL mr_io_write_parallel_coordinates(mr_io_mpi_comm, grp_id, "x", &
                                        mri_inst%x_coordinates)
  CALL mr_io_write_parallel_coordinates(mr_io_mpi_comm, grp_id, "y", &
                                        mri_inst%y_coordinates)
  CALL mr_io_write_parallel_coordinates(mr_io_mpi_comm, grp_id, "z", &
                                        mri_inst%z_coordinates)

  ! Read spatial feature
  CALL mr_io_write_parallel_spacetime_feature(mr_io_mpi_comm, &
                                           grp_id, "vector_feature", &
                                           mri_inst%vector_feature%array, &
                                           mri_inst%vector_feature%offset, &
                                           mri_inst%vector_feature%dims, &
                                           mri_inst%vector_feature%time_offset, &
                                           mri_inst%vector_feature%time_dim)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_write_parallel_spacetime


! ************************ DistFlowMRI ************************

subroutine mr_io_read_parallel_spacetime_scalar_feature(mr_io_mpi_comm, mr_io_mpi_cart_dims, grp_id, feature_name, &
                                                 feature_array, feature_offset, feature_shape, &
                                                 time_offset, time_dim) !, feature_halo_shape)

  implicit none

  INTEGER, intent(in) :: mr_io_mpi_comm
  integer, dimension(3), intent(in) :: mr_io_mpi_cart_dims
  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name

  real*8, dimension(:,:,:,:), allocatable, intent(inout) :: feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(inout) :: feature_shape  ! Shape of data subset to read/write in feature_array
  !INTEGER, DIMENSION(3), OPTIONAL, intent(in)    :: feature_halo_shape ! Shape of halo around data to read/write in feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(inout) :: feature_offset ! Offset of data subset to read/write in feature_array

  INTEGER, intent(inout) :: time_offset
  INTEGER, intent(inout) :: time_dim

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier

  INTEGER(HSIZE_T), DIMENSION(4) :: dims_file = (/-1, -1, -1, -1/) ! Shape of full HDF5 dataset
  INTEGER(HSIZE_T), DIMENSION(4) :: max_dims_file = (/-1, -1, -1, -1/) ! Max shape of full dataset

  INTEGER(HSIZE_T), DIMENSION(4) :: offset_file = (/-1, -1, -1, -1/) ! Offset of data subset to read/write in HDF5 dataset

  INTEGER(HSIZE_T), DIMENSION(4) :: dims_mem = (/-1, -1, -1, -1/) ! Shape of data subset to read/write in feature_array
  INTEGER(HSIZE_T), DIMENSION(4) :: offset_mem = (/-1, -1, -1, -1/) ! Offset of data subset to read/write in feature_array

  INTEGER(HSIZE_T), DIMENSION(4) :: count_blocks = (/-1, -1, -1, -1/) ! number of blocks to read/write
  INTEGER(HSIZE_T), DIMENSION(4) :: stride = (/-1, -1, -1, -1/) ! stride between adjacent elements/blocks?? (unused when single block is read)

  INTEGER     ::   rank = 4       ! Dataset rank
  INTEGER     ::   error          ! Error flag

  INTEGER     ::   i

  !INTEGER     ::   spatial_hyperslab_type = mr_io_parallel_spatial_hyperslab_enum_error
  INTEGER, DIMENSION(4) :: feature_array_shape_tmp = (/-1,-1,-1,-1/)
  INTEGER, DIMENSION(3) :: feature_array_shape = (/-1,-1,-1/)

  ! Open an existing dataset.
  CALL h5dopen_f(grp_id, feature_name, dset_id, error)
  mr_io_handle_error(error)

  ! Get file space
  CALL h5dget_space_f(dset_id, filespace, error)
  mr_io_handle_error(error)
  CALL h5sget_simple_extent_dims_f(filespace, dims_file, max_dims_file, error)
  mr_io_sget_simple_extent_dims_handle_error(error)


  if ( allocated(feature_array) ) then
    feature_array_shape_tmp = shape(feature_array)
    feature_array_shape = feature_array_shape_tmp(2:4)
  else ! feature_array not allocated - compute domain decomposition
    feature_array_shape = (/-1,-1,-1/)
  end if


  ! Includes argument check
  select case (mr_io_parallel_spatial_hyperslab_get_type(feature_shape, feature_offset, feature_array_shape))

  case ( mr_io_parallel_spatial_hyperslab_enum_specified )
    dims_mem = shape(feature_array)
    ! FIXME: could be imporved
    offset_file(1) = time_offset         ! temporal
    offset_file(2:4) = feature_offset    ! spatial

    if( any(dims_file(2:4) /= feature_shape) .or. (dims_file(1) /= time_dim) ) then
      write(*,*) "Array dimensions read from file and supplied in feature_shape/time_shape disagree."
      mr_io_handle_arg_error(-1)
    endif

  case ( mr_io_parallel_spatial_hyperslab_enum_compute )
    ! vectorial_dim skipped
    time_dim = dims_file(1)
    feature_shape = dims_file(2:4)

    ! temporal
    if (time_offset < 0) then
      time_offset = 0
    endif
    offset_file(1) = time_offset

    dims_mem(1) = dims_file(1)-offset_file(1)

    ! spatial
    call mr_io_parallel_spatial_hyperslap_compute(mr_io_mpi_comm, mr_io_mpi_cart_dims, feature_shape, feature_offset, dims_mem(2:4))
    ! feature_offset set in mr_io_parallel_spatial_hyperslap_compute
    offset_file(2:4) = feature_offset

    ! MRI array - TODO: handle feature_halo_shape
    allocate(feature_array(dims_mem(1), dims_mem(2), dims_mem(3), dims_mem(4)))

  case default
    mr_io_handle_arg_error(-1)

  end select

  ! TODO: handle feature_halo_shape
  offset_mem = (/ 0, 0, 0, 0 /)

  count_blocks = (/ 1, 1, 1, 1/)
  stride = (/ 1, 1, 1, 1/)

!  write (*,*) "dims_file"
!  do i=1,6
!    write (*,*) dims_file(i)
!  end do
!  write (*,*) "offset_file"
!  do i=1,6
!    write (*,*) offset_file(i)
!  end do
!  write (*,*) "dims_mem"
!  do i=1,6
!    write (*,*) dims_mem(i)
!  end do

  ! Create memory space
  CALL h5screate_simple_f(rank, dims_mem, memspace, error)
  mr_io_handle_error(error)

  ! Select hyperslab in the file.
  CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_file, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)
  CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset_mem, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)

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

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_read_parallel_spacetime_scalar_feature


subroutine mr_io_write_parallel_spacetime_scalar_feature(mr_io_mpi_comm, grp_id, feature_name, &
                                                  feature_array, feature_offset, feature_shape, &
                                                  time_offset, time_dim) !, feature_halo_shape)

  implicit none

  INTEGER, intent(in) :: mr_io_mpi_comm
  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name

  real*8, dimension(:,:,:,:), allocatable, intent(in) :: feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(in) :: feature_shape  ! Shape of data subset to read/write in feature_array
  !INTEGER, DIMENSION(3), OPTIONAL, intent(in)    :: feature_halo_shape ! Shape of halo around data to read/write in feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(in) :: feature_offset ! Offset of data subset to read/write in feature_array

  INTEGER, intent(in) :: time_offset
  INTEGER, intent(in) :: time_dim

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier

  INTEGER(HSIZE_T), DIMENSION(4) :: dims_file = (/-1, -1, -1, -1/) ! Shape of full HDF5 dataset
  INTEGER(HSIZE_T), DIMENSION(4) :: offset_file = (/-1, -1, -1, -1/) ! Offset of data subset to read/write in HDF5 dataset

  INTEGER(HSIZE_T), DIMENSION(4) :: dims_mem = (/-1, -1, -1, -1/) ! Shape of data subset to read/write in feature_array
  INTEGER(HSIZE_T), DIMENSION(4) :: offset_mem = (/-1, -1, -1, -1/) ! Offset of data subset to read/write in feature_array

  INTEGER(HSIZE_T), DIMENSION(4) :: count_blocks = (/-1, -1, -1, -1/) ! number of blocks to read/write
  INTEGER(HSIZE_T), DIMENSION(4) :: stride = (/-1, -1, -1, -1/) ! stride between adjacent elements/blocks?? (unused when single block is read)

  INTEGER     ::   rank = 4       ! Dataset rank
  INTEGER     ::   error          ! Error flag

  INTEGER     ::   i

  INTEGER, DIMENSION(4) :: feature_array_shape

  integer, dimension(4) :: arr

  feature_array_shape = shape(feature_array)

  ! Check arguments
  if ( any(feature_offset < (/0, 0, 0/)) .or. &
       ( any(feature_shape < feature_offset + feature_array_shape(2:4)) .and. any(feature_array_shape(2:4) /= (/0,0,0/)) ) .or. &
       ( time_offset < 0 ) .or. &
       ( time_dim < time_offset + feature_array_shape(1) ) ) then
    write(0,*) "spatial_global_shape: ",feature_shape
    write(0,*) "spatial_offset:       ",feature_offset
    write(0,*) "spatial_local_shape:  ",feature_array_shape(2:4)
    write(0,*) "time_global_shape:    ",time_dim
    write(0,*) "time_offset:          ",time_offset
    write(0,*) "time_local_shape:     ",feature_array_shape(1)
    write(0,*) "Invalid arguments (hyperslab must provide sufficient space for value array)"
    call flush()
    mr_io_handle_arg_error(-1)
  endif

  ! TODO: handle if( present(feature_halo_shape) )
  dims_mem = shape(feature_array)
  offset_mem = (/0,0,0,0/)

  ! Compute file dimensions and offset
  dims_file(1) = time_dim
  dims_file(2:4) = feature_shape

  offset_file(1) = time_offset
  offset_file(2:4) = feature_offset


!   write (*,*) "dims_file"
!   do i=1,3
!    write (*,*) dims_file(i)
!   end do
!   write (*,*) "feature_shape"
!   do i=1,3
!    write (*,*) feature_shape(i)
!   end do
!   write (*,*) "feature_offset"
!   do i=1,3
!    write (*,*) feature_offset(i)
!   end do
!   write (*,*) "offset_file"
!   do i=1,3
!    write (*,*) offset_file(i)
!   end do
!   write (*,*) "dims_mem"
!   do i=1,3
!    write (*,*) dims_mem(i)
!   end do
!   write (*,*) "shape(feature_array)"
!   arr = shape(feature_array)
!   do i=1,3
!     write(*,*) arr(i)
!   end do
!   write (*,*) "offset_mem"
!   do i=1,3
!    write (*,*) offset_mem(i)
!   end do
!

  ! Create the data space for the  dataset.   !
  CALL h5screate_simple_f(rank, dims_file, filespace, error)
  mr_io_handle_error(error)
  CALL h5screate_simple_f(rank, dims_mem, memspace, error)
  mr_io_handle_error(error)

  ! Create the dataset with default properties.
  CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  mr_io_handle_error(error)
  CALL h5dcreate_f(grp_id, feature_name, H5T_NATIVE_DOUBLE, filespace, &
       dset_id, error, dcpl_id=plist_id)
  mr_io_handle_error(error)
!   CALL h5pclose_f(plist_id , error)
!   mr_io_handle_error(error)

  count_blocks = (/1,1,1,1/)
  stride = (/1,1,1,1/)

  ! Select hyperslab in the file.
  CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_file, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)
  CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset_mem, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)

  ! Create property list for collective dataset read
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  mr_io_handle_error(error)
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  mr_io_handle_error(error)

  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, feature_array, dims_file, &
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

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_write_parallel_spacetime_scalar_feature


subroutine mr_io_read_parallel_spacetime_scalar_feature_padded(mr_io_mpi_comm, mr_io_mpi_cart_dims, grp_id, &
                                                 domain_padding, &
                                                 feature_name, &
                                                 feature_array, feature_offset, feature_shape, &
                                                 time_offset, time_dim) !, feature_halo_shape)

  implicit none

  INTEGER, intent(in) :: mr_io_mpi_comm
  integer, dimension(3), intent(in) :: mr_io_mpi_cart_dims
  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  type(DomainPadding) :: domain_padding
  character(len=*), intent(in) :: feature_name

  real*8, dimension(:,:,:,:), allocatable, intent(inout) :: feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(inout) :: feature_shape  ! Shape of data subset to read/write in feature_array
  !INTEGER, DIMENSION(3), OPTIONAL, intent(in)    :: feature_halo_shape ! Shape of halo around data to read/write in feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(inout) :: feature_offset ! Offset of data subset to read/write in feature_array

  INTEGER, intent(inout) :: time_offset
  INTEGER, intent(inout) :: time_dim

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier

  INTEGER(HSIZE_T), DIMENSION(4) :: dims_file = (/-1, -1, -1, -1/) ! Shape of full HDF5 dataset
  INTEGER(HSIZE_T), DIMENSION(4) :: max_dims_file = (/-1, -1, -1, -1/) ! Max shape of full dataset

  INTEGER(HSIZE_T), DIMENSION(4) :: offset_file = (/-1, -1, -1, -1/) ! Offset of data subset to read/write in HDF5 dataset

  INTEGER(HSIZE_T), DIMENSION(4) :: dims_mem = (/-1, -1, -1, -1/) ! Shape of data subset to read/write in feature_array
  INTEGER(HSIZE_T), DIMENSION(4) :: offset_mem = (/-1, -1, -1, -1/) ! Offset of data subset to read/write in feature_array
  INTEGER(HSIZE_T), DIMENSION(3) :: offset_local_hyperslab = (/0, 0, 0/) ! Offset of local hyperslab

  INTEGER(HSIZE_T), DIMENSION(4) :: count_blocks = (/-1, -1, -1, -1/) ! number of blocks to read/write
  INTEGER(HSIZE_T), DIMENSION(4) :: stride = (/-1, -1, -1, -1/) ! stride between adjacent elements/blocks?? (unused when single block is read)

  INTEGER     ::   rank = 4       ! Dataset rank
  INTEGER     ::   error          ! Error flag

  INTEGER     ::   i

  !INTEGER     ::   spatial_hyperslab_type = mr_io_parallel_spatial_hyperslab_enum_error
  INTEGER, DIMENSION(4) :: feature_array_shape_tmp = (/-1,-1,-1,-1/)
  INTEGER, DIMENSION(3) :: feature_array_shape = (/-1,-1,-1/)

  ! Open an existing dataset.
  CALL h5dopen_f(grp_id, feature_name, dset_id, error)
  mr_io_handle_error(error)

  ! Get file space
  CALL h5dget_space_f(dset_id, filespace, error)
  mr_io_handle_error(error)
  CALL h5sget_simple_extent_dims_f(filespace, dims_file, max_dims_file, error)
  mr_io_sget_simple_extent_dims_handle_error(error)


  if ( allocated(feature_array) ) then
    feature_array_shape_tmp = shape(feature_array)
    feature_array_shape = feature_array_shape_tmp(2:4)
  else ! feature_array not allocated - compute domain decomposition
    feature_array_shape = (/-1,-1,-1/)
  end if


  ! Includes argument check
  select case (mr_io_parallel_spatial_hyperslab_get_type(feature_shape, feature_offset, feature_array_shape))

  case ( mr_io_parallel_spatial_hyperslab_enum_specified )
    dims_mem = shape(feature_array)
    ! FIXME: could be imporved
    offset_file(1) = time_offset         ! temporal
    offset_file(2:4) = feature_offset    ! spatial

    if( any(dims_file(2:4) /= feature_shape) .or. (dims_file(1) /= time_dim) ) then
      write(*,*) "Array dimensions read from file and supplied in feature_shape/time_shape disagree."
      mr_io_handle_arg_error(-1)
    endif

  case ( mr_io_parallel_spatial_hyperslab_enum_compute )
    ! vectorial_dim skipped
    time_dim = dims_file(1)
    feature_shape = dims_file(2:4)

    ! temporal
    if (time_offset < 0) then
      time_offset = 0
    endif
    offset_file(1) = time_offset

    dims_mem(1) = dims_file(1)-offset_file(1)

    ! spatial
    call mr_io_parallel_spatial_hyperslap_compute_padded(mr_io_mpi_comm, mr_io_mpi_cart_dims, domain_padding, &
                                                         feature_shape, feature_offset, dims_mem(2:4), offset_local_hyperslab)
    ! feature_offset set in mr_io_parallel_spatial_hyperslap_compute
    offset_file(2:4) = feature_offset

    ! MRI array - TODO: handle feature_halo_shape
    allocate(feature_array(dims_mem(1), &
                           offset_local_hyperslab(1)+1:offset_local_hyperslab(1)+dims_mem(2), &
                           offset_local_hyperslab(2)+1:offset_local_hyperslab(2)+dims_mem(3), &
                           offset_local_hyperslab(3)+1:offset_local_hyperslab(3)+dims_mem(4)))
!    write(0,*) "shape of local feature array: ",shape(feature_array)

  case default
    mr_io_handle_arg_error(-1)

  end select

  ! TODO: handle feature_halo_shape
  offset_mem = (/ 0, 0, 0, 0 /)

  count_blocks = (/ 1, 1, 1, 1/)
  stride = (/ 1, 1, 1, 1/)

!  write (*,*) "dims_file"
!  do i=1,6
!    write (*,*) dims_file(i)
!  end do
!  write (*,*) "offset_file"
!  do i=1,6
!    write (*,*) offset_file(i)
!  end do
!  write (*,*) "dims_mem"
!  do i=1,6
!    write (*,*) dims_mem(i)
!  end do

  ! Create memory space
  CALL h5screate_simple_f(rank, dims_mem, memspace, error)
  mr_io_handle_error(error)

  ! Select hyperslab in the file.
  CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_file, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)
  CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset_mem, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)

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

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_read_parallel_spacetime_scalar_feature_padded


subroutine mr_io_read_parallel_spacetime_matrix_feature(mr_io_mpi_comm, mr_io_mpi_cart_dims, grp_id, feature_name, &
                                                 feature_array, feature_offset, feature_shape, &
                                                 time_offset, time_dim) !, feature_halo_shape)

  implicit none

  INTEGER, intent(in) :: mr_io_mpi_comm
  integer, dimension(3), intent(in) :: mr_io_mpi_cart_dims
  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name

  real*8, dimension(:,:,:,:,:,:), allocatable, intent(inout) :: feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(inout) :: feature_shape  ! Shape of data subset to read/write in feature_array
  !INTEGER, DIMENSION(3), OPTIONAL, intent(in)    :: feature_halo_shape ! Shape of halo around data to read/write in feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(inout) :: feature_offset ! Offset of data subset to read/write in feature_array

  INTEGER, intent(inout) :: time_offset
  INTEGER, intent(inout) :: time_dim

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier

  INTEGER(HSIZE_T), DIMENSION(6) :: dims_file = (/-1, -1, -1, -1, -1, -1/) ! Shape of full HDF5 dataset
  INTEGER(HSIZE_T), DIMENSION(6) :: max_dims_file = (/-1, -1, -1, -1, -1, -1/) ! Max shape of full dataset

  INTEGER(HSIZE_T), DIMENSION(6) :: offset_file = (/-1, -1, -1, -1, -1, -1/) ! Offset of data subset to read/write in HDF5 dataset

  INTEGER(HSIZE_T), DIMENSION(6) :: dims_mem = (/-1, -1, -1, -1, -1, -1/) ! Shape of data subset to read/write in feature_array
  INTEGER(HSIZE_T), DIMENSION(6) :: offset_mem = (/-1, -1, -1, -1, -1, -1/) ! Offset of data subset to read/write in feature_array

  INTEGER(HSIZE_T), DIMENSION(6) :: count_blocks = (/-1, -1, -1, -1, -1, -1/) ! number of blocks to read/write
  INTEGER(HSIZE_T), DIMENSION(6) :: stride = (/-1, -1, -1, -1, -1, -1/) ! stride between adjacent elements/blocks?? (unused when single block is read)

  INTEGER     ::   rank = 6       ! Dataset rank
  INTEGER     ::   error          ! Error flag

  INTEGER     ::   i

  !INTEGER     ::   spatial_hyperslab_type = mr_io_parallel_spatial_hyperslab_enum_error
  INTEGER, DIMENSION(6) :: feature_array_shape_tmp = (/-1,-1,-1,-1,-1,-1/)
  INTEGER, DIMENSION(3) :: feature_array_shape = (/-1,-1,-1/)

  ! Open an existing dataset.
  CALL h5dopen_f(grp_id, feature_name, dset_id, error)
  mr_io_handle_error(error)

  ! Get file space
  CALL h5dget_space_f(dset_id, filespace, error)
  mr_io_handle_error(error)
  CALL h5sget_simple_extent_dims_f(filespace, dims_file, max_dims_file, error)
  mr_io_sget_simple_extent_dims_handle_error(error)


  if ( allocated(feature_array) ) then
    feature_array_shape_tmp = shape(feature_array)
    feature_array_shape = feature_array_shape_tmp(4:6)
  else ! feature_array not allocated - compute domain decomposition
    feature_array_shape = (/-1,-1,-1/)
  end if


  ! Includes argument check
  select case (mr_io_parallel_spatial_hyperslab_get_type(feature_shape, feature_offset, feature_array_shape))

  case ( mr_io_parallel_spatial_hyperslab_enum_specified )
    dims_mem = shape(feature_array)
    ! FIXME: could be imporved
    offset_file(1:2) = 0                 ! vectorial
    offset_file(3) = time_offset         ! temporal
    offset_file(4:6) = feature_offset    ! spatial

    if( any(dims_file(4:6) /= feature_shape) .or. (dims_file(3) /= time_dim) ) then
      write(*,*) "Array dimensions read from file and supplied in feature_shape/time_shape disagree."
      mr_io_handle_arg_error(-1)
    endif

  case ( mr_io_parallel_spatial_hyperslab_enum_compute )
    ! vectorial_dim skipped
    time_dim = dims_file(3)
    feature_shape = dims_file(4:6)

    ! vectorial
    offset_file(1:2) = 0
    dims_mem(1:2) = dims_file(1:2)

    ! temporal
    if (time_offset < 0) then
      time_offset = 0
    endif
    offset_file(3) = time_offset

    dims_mem(3) = dims_file(3)-offset_file(3)

    ! spatial
    call mr_io_parallel_spatial_hyperslap_compute(mr_io_mpi_comm, mr_io_mpi_cart_dims, feature_shape, feature_offset, dims_mem(4:6))
    ! feature_offset set in mr_io_parallel_spatial_hyperslap_compute
    offset_file(4:6) = feature_offset

    ! MRI array - TODO: handle feature_halo_shape
    allocate(feature_array(dims_mem(1), dims_mem(2), dims_mem(3), dims_mem(4), dims_mem(5), dims_mem(6)))

  case default
    mr_io_handle_arg_error(-1)

  end select

  ! TODO: handle feature_halo_shape
  offset_mem = (/ 0, 0, 0, 0 ,0 ,0/)

  count_blocks = (/ 1, 1, 1, 1, 1, 1/)
  stride = (/ 1, 1, 1, 1, 1, 1/)

!  write (*,*) "dims_file"
!  do i=1,6
!    write (*,*) dims_file(i)
!  end do
!  write (*,*) "offset_file"
!  do i=1,6
!    write (*,*) offset_file(i)
!  end do
!  write (*,*) "dims_mem"
!  do i=1,6
!    write (*,*) dims_mem(i)
!  end do

  ! Create memory space
  CALL h5screate_simple_f(rank, dims_mem, memspace, error)
  mr_io_handle_error(error)

  ! Select hyperslab in the file.
  CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_file, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)
  CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset_mem, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)

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

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_read_parallel_spacetime_matrix_feature


subroutine mr_io_write_parallel_spacetime_matrix_feature(mr_io_mpi_comm, grp_id, feature_name, &
                                                  feature_array, feature_offset, feature_shape, &
                                                  time_offset, time_dim) !, feature_halo_shape)

  implicit none

  INTEGER, intent(in) :: mr_io_mpi_comm
  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name

  real*8, dimension(:,:,:,:,:,:), allocatable, intent(in) :: feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(in) :: feature_shape  ! Shape of data subset to read/write in feature_array
  !INTEGER, DIMENSION(3), OPTIONAL, intent(in)    :: feature_halo_shape ! Shape of halo around data to read/write in feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(in) :: feature_offset ! Offset of data subset to read/write in feature_array

  INTEGER, intent(in) :: time_offset
  INTEGER, intent(in) :: time_dim

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier

  INTEGER(HSIZE_T), DIMENSION(6) :: dims_file = (/-1, -1, -1, -1, -1, -1/) ! Shape of full HDF5 dataset
  INTEGER(HSIZE_T), DIMENSION(6) :: offset_file = (/-1, -1, -1, -1, -1, -1/) ! Offset of data subset to read/write in HDF5 dataset

  INTEGER(HSIZE_T), DIMENSION(6) :: dims_mem = (/-1, -1, -1, -1, -1, -1/) ! Shape of data subset to read/write in feature_array
  INTEGER(HSIZE_T), DIMENSION(6) :: offset_mem = (/-1, -1, -1, -1, -1, -1/) ! Offset of data subset to read/write in feature_array

  INTEGER(HSIZE_T), DIMENSION(6) :: count_blocks = (/-1, -1, -1, -1, -1, -1/) ! number of blocks to read/write
  INTEGER(HSIZE_T), DIMENSION(6) :: stride = (/-1, -1, -1, -1, -1, -1/) ! stride between adjacent elements/blocks?? (unused when single block is read)

  INTEGER     ::   rank = 6       ! Dataset rank
  INTEGER     ::   error          ! Error flag

  INTEGER     ::   i

  INTEGER, DIMENSION(6) :: feature_array_shape

  integer, dimension(6) :: arr

  feature_array_shape = shape(feature_array)

  ! Check arguments
  if ( any(feature_offset < (/0, 0, 0/)) .or. &
       ( any(feature_shape < feature_offset + feature_array_shape(4:6)) .and. any(feature_array_shape(4:6) /= (/0,0,0/)) ) .or. &
       ( time_offset < 0 ) .or. &
       ( time_dim < time_offset + feature_array_shape(3) ) ) then
    write(0,*) "spatial_global_shape: ",feature_shape
    write(0,*) "spatial_offset:       ",feature_offset
    write(0,*) "spatial_local_shape:  ",feature_array_shape(4:6)
    write(0,*) "time_global_shape:    ",time_dim
    write(0,*) "time_offset:          ",time_offset
    write(0,*) "time_local_shape:     ",feature_array_shape(3)
    write(0,*) "Invalid arguments (hyperslab must provide sufficient space for value array)"
    call flush()
    mr_io_handle_arg_error(-1)
  endif

  ! TODO: handle if( present(feature_halo_shape) )
  dims_mem = shape(feature_array)
  offset_mem = (/0,0,0,0,0,0/)

  ! Compute file dimensions and offset
  dims_file(1:2) = dims_mem(1:2)
  dims_file(3) = time_dim
  dims_file(4:6) = feature_shape

  offset_file(1:2) = 0
  offset_file(3) = time_offset
  offset_file(4:6) = feature_offset


!   write (*,*) "dims_file"
!   do i=1,3
!    write (*,*) dims_file(i)
!   end do
!   write (*,*) "feature_shape"
!   do i=1,3
!    write (*,*) feature_shape(i)
!   end do
!   write (*,*) "feature_offset"
!   do i=1,3
!    write (*,*) feature_offset(i)
!   end do
!   write (*,*) "offset_file"
!   do i=1,3
!    write (*,*) offset_file(i)
!   end do
!   write (*,*) "dims_mem"
!   do i=1,3
!    write (*,*) dims_mem(i)
!   end do
!   write (*,*) "shape(feature_array)"
!   arr = shape(feature_array)
!   do i=1,3
!     write(*,*) arr(i)
!   end do
!   write (*,*) "offset_mem"
!   do i=1,3
!    write (*,*) offset_mem(i)
!   end do
!

  ! Create the data space for the  dataset.   !
  CALL h5screate_simple_f(rank, dims_file, filespace, error)
  mr_io_handle_error(error)
  CALL h5screate_simple_f(rank, dims_mem, memspace, error)
  mr_io_handle_error(error)

  ! Create the dataset with default properties.
  CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  mr_io_handle_error(error)
  CALL h5dcreate_f(grp_id, feature_name, H5T_NATIVE_DOUBLE, filespace, &
       dset_id, error, dcpl_id=plist_id)
  mr_io_handle_error(error)
!   CALL h5pclose_f(plist_id , error)
!   mr_io_handle_error(error)

  count_blocks = (/1,1,1,1,1,1/)
  stride = (/1,1,1,1,1,1/)

  ! Select hyperslab in the file.
  CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_file, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)
  CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset_mem, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)

  ! Create property list for collective dataset read
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  mr_io_handle_error(error)
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  mr_io_handle_error(error)

  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, feature_array, dims_file, &
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

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_write_parallel_spacetime_matrix_feature


subroutine mr_io_read_parallel_spacetime_matrix_feature_padded(mr_io_mpi_comm, mr_io_mpi_cart_dims, grp_id, &
                                                 domain_padding, &
                                                 feature_name, &
                                                 feature_array, feature_offset, feature_shape, &
                                                 time_offset, time_dim) !, feature_halo_shape)

  implicit none

  INTEGER, intent(in) :: mr_io_mpi_comm
  integer, dimension(3), intent(in) :: mr_io_mpi_cart_dims
  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  type(DomainPadding) :: domain_padding
  character(len=*), intent(in) :: feature_name

  real*8, dimension(:,:,:,:,:,:), allocatable, intent(inout) :: feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(inout) :: feature_shape  ! Shape of data subset to read/write in feature_array
  !INTEGER, DIMENSION(3), OPTIONAL, intent(in)    :: feature_halo_shape ! Shape of halo around data to read/write in feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(inout) :: feature_offset ! Offset of data subset to read/write in feature_array

  INTEGER, intent(inout) :: time_offset
  INTEGER, intent(inout) :: time_dim

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier

  INTEGER(HSIZE_T), DIMENSION(6) :: dims_file = (/-1, -1, -1, -1, -1, -1/) ! Shape of full HDF5 dataset
  INTEGER(HSIZE_T), DIMENSION(6) :: max_dims_file = (/-1, -1, -1, -1, -1, -1/) ! Max shape of full dataset

  INTEGER(HSIZE_T), DIMENSION(6) :: offset_file = (/-1, -1, -1, -1, -1, -1/) ! Offset of data subset to read/write in HDF5 dataset

  INTEGER(HSIZE_T), DIMENSION(6) :: dims_mem = (/-1, -1, -1, -1, -1, -1/) ! Shape of data subset to read/write in feature_array
  INTEGER(HSIZE_T), DIMENSION(6) :: offset_mem = (/-1, -1, -1, -1, -1, -1/) ! Offset of data subset to read/write in feature_array
  INTEGER(HSIZE_T), DIMENSION(3) :: offset_local_hyperslab = (/0, 0, 0/) ! Offset of local hyperslab

  INTEGER(HSIZE_T), DIMENSION(6) :: count_blocks = (/-1, -1, -1, -1, -1, -1/) ! number of blocks to read/write
  INTEGER(HSIZE_T), DIMENSION(6) :: stride = (/-1, -1, -1, -1, -1, -1/) ! stride between adjacent elements/blocks?? (unused when single block is read)

  INTEGER     ::   rank = 6       ! Dataset rank
  INTEGER     ::   error          ! Error flag

  INTEGER     ::   i

  !INTEGER     ::   spatial_hyperslab_type = mr_io_parallel_spatial_hyperslab_enum_error
  INTEGER, DIMENSION(6) :: feature_array_shape_tmp = (/-1,-1,-1,-1,-1,-1/)
  INTEGER, DIMENSION(3) :: feature_array_shape = (/-1,-1,-1/)

  ! Open an existing dataset.
  CALL h5dopen_f(grp_id, feature_name, dset_id, error)
  mr_io_handle_error(error)

  ! Get file space
  CALL h5dget_space_f(dset_id, filespace, error)
  mr_io_handle_error(error)
  CALL h5sget_simple_extent_dims_f(filespace, dims_file, max_dims_file, error)
  mr_io_sget_simple_extent_dims_handle_error(error)


  if ( allocated(feature_array) ) then
    feature_array_shape_tmp = shape(feature_array)
    feature_array_shape = feature_array_shape_tmp(4:6)
  else ! feature_array not allocated - compute domain decomposition
    feature_array_shape = (/-1,-1,-1/)
  end if


  ! Includes argument check
  select case (mr_io_parallel_spatial_hyperslab_get_type(feature_shape, feature_offset, feature_array_shape))

  case ( mr_io_parallel_spatial_hyperslab_enum_specified )
    dims_mem = shape(feature_array)
    ! FIXME: could be imporved
    offset_file(1:2) = 0                 ! vectorial
    offset_file(3) = time_offset         ! temporal
    offset_file(4:6) = feature_offset    ! spatial

    if( any(dims_file(4:6) /= feature_shape) .or. (dims_file(3) /= time_dim) ) then
      write(*,*) "Array dimensions read from file and supplied in feature_shape/time_shape disagree."
      mr_io_handle_arg_error(-1)
    endif

  case ( mr_io_parallel_spatial_hyperslab_enum_compute )
    ! vectorial_dim skipped
    time_dim = dims_file(3)
    feature_shape = dims_file(4:6)

    ! vectorial
    offset_file(1:2) = 0
    dims_mem(1:2) = dims_file(1:2)

    ! temporal
    if (time_offset < 0) then
      time_offset = 0
    endif
    offset_file(3) = time_offset

    dims_mem(3) = dims_file(3)-offset_file(3)

    ! spatial
    call mr_io_parallel_spatial_hyperslap_compute_padded(mr_io_mpi_comm, mr_io_mpi_cart_dims, domain_padding, &
                                                         feature_shape, feature_offset, dims_mem(4:6), offset_local_hyperslab)
    ! feature_offset set in mr_io_parallel_spatial_hyperslap_compute
    offset_file(4:6) = feature_offset

    ! MRI array - TODO: handle feature_halo_shape
    allocate(feature_array(dims_mem(1), dims_mem(2), dims_mem(3), &
                           offset_local_hyperslab(1)+1:offset_local_hyperslab(1)+dims_mem(4), &
                           offset_local_hyperslab(2)+1:offset_local_hyperslab(2)+dims_mem(5), &
                           offset_local_hyperslab(3)+1:offset_local_hyperslab(3)+dims_mem(6)))
!    write(0,*) "shape of local feature array: ",shape(feature_array)

  case default
    mr_io_handle_arg_error(-1)

  end select

  ! TODO: handle feature_halo_shape
  offset_mem = (/ 0, 0, 0, 0 ,0 ,0/)

  count_blocks = (/ 1, 1, 1, 1, 1, 1/)
  stride = (/ 1, 1, 1, 1, 1, 1/)

!  write (*,*) "dims_file"
!  do i=1,6
!    write (*,*) dims_file(i)
!  end do
!  write (*,*) "offset_file"
!  do i=1,6
!    write (*,*) offset_file(i)
!  end do
!  write (*,*) "dims_mem"
!  do i=1,6
!    write (*,*) dims_mem(i)
!  end do

  ! Create memory space
  CALL h5screate_simple_f(rank, dims_mem, memspace, error)
  mr_io_handle_error(error)

  ! Select hyperslab in the file.
  CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_file, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)
  CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset_mem, &
                             count_blocks, error, stride, dims_mem) ! default values for stride, block
  mr_io_handle_error(error)

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

  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  mr_io_handle_error(error)

end subroutine mr_io_read_parallel_spacetime_matrix_feature_padded


subroutine mr_io_read_parallel_flow(mr_io_mpi_comm, mr_io_mpi_info, mr_io_mpi_cart_dims, path, mri_inst)

  implicit none

  character(len=*), intent(in) :: path
  INTEGER, intent (in) :: mr_io_mpi_comm, mr_io_mpi_info
  integer, dimension(3), intent(in) :: mr_io_mpi_cart_dims
  type(DistFlowMRI), intent(out) :: mri_inst

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

  CALL h5pset_fapl_mpio_f(plist_id, mr_io_mpi_comm, mr_io_mpi_info, error)
  mr_io_handle_error(error)

  ! Open existing file collectively
  CALL h5fopen_f(trim(path), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
  mr_io_handle_error(error)

  CALL h5pclose_f(plist_id, error)
  mr_io_handle_error(error)

  ! Open an existing group
  CALL h5gopen_f(file_id, FlowMRI_group_name, grp_id, error)
  mr_io_handle_error(error)

  ! Read time coordinates
  CALL mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, "t", &
                                       mri_inst%t_coordinates)

  ! Read coordinates
  CALL mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, "x", &
                                       mri_inst%x_coordinates)
  CALL mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, "y", &
                                       mri_inst%y_coordinates)
  CALL mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, "z", &
                                       mri_inst%z_coordinates)


  ! Read spatial feature
  CALL mr_io_read_parallel_spacetime_scalar_feature(mr_io_mpi_comm, &
                                             mr_io_mpi_cart_dims, &
                                             grp_id, "intensity", &
                                             mri_inst%intensity%array, &
                                             mri_inst%intensity%offset, &
                                             mri_inst%intensity%dims, &
                                             mri_inst%intensity%time_offset, &
                                             mri_inst%intensity%time_dim)
  CALL mr_io_read_parallel_spacetime_feature(mr_io_mpi_comm, &
                                             mr_io_mpi_cart_dims, &
                                             grp_id, "velocity_mean", &
                                             mri_inst%velocity_mean%array, &
                                             mri_inst%velocity_mean%offset, &
                                             mri_inst%velocity_mean%dims, &
                                             mri_inst%velocity_mean%time_offset, &
                                             mri_inst%velocity_mean%time_dim)
  CALL mr_io_read_parallel_spacetime_matrix_feature(mr_io_mpi_comm, &
                                             mr_io_mpi_cart_dims, &
                                             grp_id, "velocity_cov", &
                                             mri_inst%velocity_cov%array, &
                                             mri_inst%velocity_cov%offset, &
                                             mri_inst%velocity_cov%dims, &
                                             mri_inst%velocity_cov%time_offset, &
                                             mri_inst%velocity_cov%time_dim)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_read_parallel_flow



subroutine mr_io_write_parallel_flow(mr_io_mpi_comm, mr_io_mpi_info, path, mri_inst)

  implicit none

  character(len=*), intent(in) :: path
  INTEGER, intent (in) :: mr_io_mpi_comm, mr_io_mpi_info
  type(DistFlowMRI), intent(in) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier

  INTEGER(HID_T) :: plist_id      ! Property list identifier

  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  mr_io_handle_error(error)

  ! Create a file collectively
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  mr_io_handle_error(error)

  CALL h5pset_fapl_mpio_f(plist_id, mr_io_mpi_comm, mr_io_mpi_info, error)
  mr_io_handle_error(error)

  CALL h5fcreate_f(trim(path), H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
  mr_io_handle_error(error)

  CALL h5pclose_f(plist_id, error)
  mr_io_handle_error(error)

  ! Create a new group collectively
  CALL h5pcreate_f(H5P_GROUP_CREATE_F, plist_id, error) ! FIXME: Not completely sure about flag
  mr_io_handle_error(error)

  CALL h5gcreate_f(file_id, FlowMRI_group_name, grp_id, error, gcpl_id=plist_id)
  mr_io_handle_error(error)

  CALL h5pclose_f(plist_id, error)
  mr_io_handle_error(error)

!   print *, mri_inst%vector_feature%dims
!   print *, mri_inst%vector_feature%offset
!   print *, shape(mri_inst%vector_feature%array)
  CALL mr_io_write_parallel_coordinates(mr_io_mpi_comm, grp_id, "t", &
                                        mri_inst%t_coordinates)
  CALL mr_io_write_parallel_coordinates(mr_io_mpi_comm, grp_id, "x", &
                                        mri_inst%x_coordinates)
  CALL mr_io_write_parallel_coordinates(mr_io_mpi_comm, grp_id, "y", &
                                        mri_inst%y_coordinates)
  CALL mr_io_write_parallel_coordinates(mr_io_mpi_comm, grp_id, "z", &
                                        mri_inst%z_coordinates)

  ! Read spatial feature
  CALL mr_io_write_parallel_spacetime_scalar_feature(mr_io_mpi_comm, grp_id, "intensity", &
                                           mri_inst%intensity%array, &
                                           mri_inst%intensity%offset, &
                                           mri_inst%intensity%dims, &
                                           mri_inst%intensity%time_offset, &
                                           mri_inst%intensity%time_dim)
  CALL mr_io_write_parallel_spacetime_feature(mr_io_mpi_comm, grp_id, "velocity_mean", &
                                           mri_inst%velocity_mean%array, &
                                           mri_inst%velocity_mean%offset, &
                                           mri_inst%velocity_mean%dims, &
                                           mri_inst%velocity_mean%time_offset, &
                                           mri_inst%velocity_mean%time_dim)
  CALL mr_io_write_parallel_spacetime_matrix_feature(mr_io_mpi_comm, grp_id, "velocity_cov", &
                                           mri_inst%velocity_cov%array, &
                                           mri_inst%velocity_cov%offset, &
                                           mri_inst%velocity_cov%dims, &
                                           mri_inst%velocity_cov%time_offset, &
                                           mri_inst%velocity_cov%time_dim)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_write_parallel_flow


subroutine mr_io_read_parallel_flow_padded(mr_io_mpi_comm, mr_io_mpi_info, mr_io_mpi_cart_dims, path, mri_inst_padded)

  implicit none

  character(len=*), intent(in) :: path
  INTEGER, intent (in) :: mr_io_mpi_comm, mr_io_mpi_info
  integer, dimension(3), intent(in) :: mr_io_mpi_cart_dims
  type(DistFlowMRIPadded), intent(inout) :: mri_inst_padded
!  type(DistFlowMRI) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier

  INTEGER(HID_T) :: plist_id      ! Property list identifier

  INTEGER     ::   error ! Error flag

!  mri_inst = mri_inst_padded%mri

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  mr_io_handle_error(error)

  ! Setup file access property list with parallel I/O access.
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  mr_io_handle_error(error)

  CALL h5pset_fapl_mpio_f(plist_id, mr_io_mpi_comm, mr_io_mpi_info, error)
  mr_io_handle_error(error)

  ! Open existing file collectively
  CALL h5fopen_f(trim(path), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
  mr_io_handle_error(error)

  CALL h5pclose_f(plist_id, error)
  mr_io_handle_error(error)

  ! Open an existing group
  CALL h5gopen_f(file_id, FlowMRI_group_name, grp_id, error)
  mr_io_handle_error(error)

  ! Read time coordinates
  CALL mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, "t", &
                                       mri_inst_padded%mri%t_coordinates)

  ! Read coordinates
  CALL mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, "x", &
                                       mri_inst_padded%mri%x_coordinates)
  CALL mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, "y", &
                                       mri_inst_padded%mri%y_coordinates)
  CALL mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, "z", &
                                       mri_inst_padded%mri%z_coordinates)


  ! Read spatial feature
  CALL mr_io_read_parallel_spacetime_scalar_feature_padded(mr_io_mpi_comm, &
                                             mr_io_mpi_cart_dims, &
                                             grp_id, &
                                             mri_inst_padded%domain_padding, &
                                             "intensity", &
                                             mri_inst_padded%mri%intensity%array, &
                                             mri_inst_padded%mri%intensity%offset, &
                                             mri_inst_padded%mri%intensity%dims, &
                                             mri_inst_padded%mri%intensity%time_offset, &
                                             mri_inst_padded%mri%intensity%time_dim)
  CALL mr_io_read_parallel_spacetime_feature_padded(mr_io_mpi_comm, &
                                             mr_io_mpi_cart_dims, &
                                             grp_id, &
                                             mri_inst_padded%domain_padding, &
                                             "velocity_mean", &
                                             mri_inst_padded%mri%velocity_mean%array, &
                                             mri_inst_padded%mri%velocity_mean%offset, &
                                             mri_inst_padded%mri%velocity_mean%dims, &
                                             mri_inst_padded%mri%velocity_mean%time_offset, &
                                             mri_inst_padded%mri%velocity_mean%time_dim)
  CALL mr_io_read_parallel_spacetime_matrix_feature_padded(mr_io_mpi_comm, &
                                             mr_io_mpi_cart_dims, &
                                             grp_id, &
                                             mri_inst_padded%domain_padding, &
                                             "velocity_cov", &
                                             mri_inst_padded%mri%velocity_cov%array, &
                                             mri_inst_padded%mri%velocity_cov%offset, &
                                             mri_inst_padded%mri%velocity_cov%dims, &
                                             mri_inst_padded%mri%velocity_cov%time_offset, &
                                             mri_inst_padded%mri%velocity_cov%time_dim)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_read_parallel_flow_padded


! ************************ DistSegmentedFlowMRI ************************

subroutine mr_io_read_parallel_segmentedflow(mr_io_mpi_comm, mr_io_mpi_info, mr_io_mpi_cart_dims, path, mri_inst)

  implicit none

  character(len=*), intent(in) :: path
  INTEGER, intent (in) :: mr_io_mpi_comm, mr_io_mpi_info
  integer, dimension(3), intent(in) :: mr_io_mpi_cart_dims
  type(DistSegmentedFlowMRI), intent(out) :: mri_inst

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

  CALL h5pset_fapl_mpio_f(plist_id, mr_io_mpi_comm, mr_io_mpi_info, error)
  mr_io_handle_error(error)

  ! Open existing file collectively
  CALL h5fopen_f(trim(path), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
  mr_io_handle_error(error)

  CALL h5pclose_f(plist_id, error)
  mr_io_handle_error(error)

  ! Open an existing group
  CALL h5gopen_f(file_id, SegmentedFlowMRI_group_name, grp_id, error)
  mr_io_handle_error(error)

  ! Read time coordinates
  CALL mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, "t", &
                                       mri_inst%t_coordinates)

  ! Read coordinates
  CALL mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, "x", &
                                       mri_inst%x_coordinates)
  CALL mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, "y", &
                                       mri_inst%y_coordinates)
  CALL mr_io_read_parallel_coordinates(mr_io_mpi_comm, grp_id, "z", &
                                       mri_inst%z_coordinates)


  ! Read spatial feature
  CALL mr_io_read_parallel_spacetime_scalar_feature(mr_io_mpi_comm, &
                                             mr_io_mpi_cart_dims, &
                                             grp_id, "intensity", &
                                             mri_inst%intensity%array, &
                                             mri_inst%intensity%offset, &
                                             mri_inst%intensity%dims, &
                                             mri_inst%intensity%time_offset, &
                                             mri_inst%intensity%time_dim)
  CALL mr_io_read_parallel_spacetime_feature(mr_io_mpi_comm, &
                                             mr_io_mpi_cart_dims, &
                                             grp_id, "velocity_mean", &
                                             mri_inst%velocity_mean%array, &
                                             mri_inst%velocity_mean%offset, &
                                             mri_inst%velocity_mean%dims, &
                                             mri_inst%velocity_mean%time_offset, &
                                             mri_inst%velocity_mean%time_dim)
  CALL mr_io_read_parallel_spacetime_matrix_feature(mr_io_mpi_comm, &
                                             mr_io_mpi_cart_dims, &
                                             grp_id, "velocity_cov", &
                                             mri_inst%velocity_cov%array, &
                                             mri_inst%velocity_cov%offset, &
                                             mri_inst%velocity_cov%dims, &
                                             mri_inst%velocity_cov%time_offset, &
                                             mri_inst%velocity_cov%time_dim)
  CALL mr_io_read_parallel_spacetime_scalar_feature(mr_io_mpi_comm, &
                                             mr_io_mpi_cart_dims, &
                                             grp_id, "segmentation_prob", &
                                             mri_inst%segmentation_prob%array, &
                                             mri_inst%segmentation_prob%offset, &
                                             mri_inst%segmentation_prob%dims, &
                                             mri_inst%segmentation_prob%time_offset, &
                                             mri_inst%segmentation_prob%time_dim)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_read_parallel_segmentedflow



subroutine mr_io_write_parallel_segmentedflow(mr_io_mpi_comm, mr_io_mpi_info, path, mri_inst)

  implicit none

  character(len=*), intent(in) :: path
  INTEGER, intent (in) :: mr_io_mpi_comm, mr_io_mpi_info
  type(DistSegmentedFlowMRI), intent(in) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier

  INTEGER(HID_T) :: plist_id      ! Property list identifier

  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  mr_io_handle_error(error)

  ! Create a file collectively
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  mr_io_handle_error(error)

  CALL h5pset_fapl_mpio_f(plist_id, mr_io_mpi_comm, mr_io_mpi_info, error)
  mr_io_handle_error(error)

  CALL h5fcreate_f(trim(path), H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
  mr_io_handle_error(error)

  CALL h5pclose_f(plist_id, error)
  mr_io_handle_error(error)

  ! Create a new group collectively
  CALL h5pcreate_f(H5P_GROUP_CREATE_F, plist_id, error) ! FIXME: Not completely sure about flag
  mr_io_handle_error(error)

  CALL h5gcreate_f(file_id, SegmentedFlowMRI_group_name, grp_id, error, gcpl_id=plist_id)
  mr_io_handle_error(error)

  CALL h5pclose_f(plist_id, error)
  mr_io_handle_error(error)

!   print *, mri_inst%vector_feature%dims
!   print *, mri_inst%vector_feature%offset
!   print *, shape(mri_inst%vector_feature%array)
  CALL mr_io_write_parallel_coordinates(mr_io_mpi_comm, grp_id, "t", &
                                        mri_inst%t_coordinates)
  CALL mr_io_write_parallel_coordinates(mr_io_mpi_comm, grp_id, "x", &
                                        mri_inst%x_coordinates)
  CALL mr_io_write_parallel_coordinates(mr_io_mpi_comm, grp_id, "y", &
                                        mri_inst%y_coordinates)
  CALL mr_io_write_parallel_coordinates(mr_io_mpi_comm, grp_id, "z", &
                                        mri_inst%z_coordinates)

  ! Read spatial feature
  CALL mr_io_write_parallel_spacetime_scalar_feature(mr_io_mpi_comm, grp_id, "intensity", &
                                           mri_inst%intensity%array, &
                                           mri_inst%intensity%offset, &
                                           mri_inst%intensity%dims, &
                                           mri_inst%intensity%time_offset, &
                                           mri_inst%intensity%time_dim)
  CALL mr_io_write_parallel_spacetime_feature(mr_io_mpi_comm, grp_id, "velocity_mean", &
                                           mri_inst%velocity_mean%array, &
                                           mri_inst%velocity_mean%offset, &
                                           mri_inst%velocity_mean%dims, &
                                           mri_inst%velocity_mean%time_offset, &
                                           mri_inst%velocity_mean%time_dim)
  CALL mr_io_write_parallel_spacetime_matrix_feature(mr_io_mpi_comm, grp_id, "velocity_cov", &
                                           mri_inst%velocity_cov%array, &
                                           mri_inst%velocity_cov%offset, &
                                           mri_inst%velocity_cov%dims, &
                                           mri_inst%velocity_cov%time_offset, &
                                           mri_inst%velocity_cov%time_dim)
  CALL mr_io_write_parallel_spacetime_scalar_feature(mr_io_mpi_comm, grp_id, "segmentation_prob", &
                                           mri_inst%segmentation_prob%array, &
                                           mri_inst%segmentation_prob%offset, &
                                           mri_inst%segmentation_prob%dims, &
                                           mri_inst%segmentation_prob%time_offset, &
                                           mri_inst%segmentation_prob%time_dim)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
  CALL h5fclose_f(file_id, error)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_write_parallel_segmentedflow


end module mr_io_parallel_spacetime
