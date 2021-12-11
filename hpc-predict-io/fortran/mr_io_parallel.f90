module mr_io_parallel

use hdf5
use mr_io_protocol
use mr_io, only : mr_io_handle_hdf5_error, mr_io_handle_argument_error
use mr_io_locking_utils, only: mr_io_h5_parallel_reader_open_f, &
                               mr_io_h5_parallel_reader_close_f, &
                               mr_io_h5_parallel_writer_open_f, &
                               mr_io_h5_parallel_writer_close_f, &
                               stdin, stdout, stderr

use mpi !include 'mpif.h'

implicit none
private

public :: mr_io_read_parallel_spatial
public :: mr_io_write_parallel_spatial

public :: mr_io_parallel_spatial_hyperslap_compute
public :: mr_io_parallel_spatial_hyperslap_compute_padded
public :: mr_io_parallel_spatial_hyperslab_get_type

public :: DistSpatialMRI
public :: SpatialMRI_group_name
public :: mr_io_deallocate_dist_spatial_mri

!#ifdef __GFORTRAN__
!! Spatial hyperslab macro enum
!! Compute spatial hyperslab
!#define mr_io_parallel_spatial_hyperslab_enum_compute   ( 0)
!! Read specified hyperslab
!#define mr_io_parallel_spatial_hyperslab_enum_specified ( 1)
!! Inconsistent arguments (shape/offset/feature_array)
!#define mr_io_parallel_spatial_hyperslab_enum_error     (-1)
!#endif

#include "mr_io_hyperslab_enum.f90"

contains

#include "mr_io_mpi_error_handler.f90"


! ************************ Domain/Hyperslab decomposition ************************

! This does unsymmetric positioning with padding at the RHS of the domain. For symmetric alignment use the padded version
! with 0 domain padding.
subroutine mr_io_parallel_spatial_hyperslap_compute(mr_io_mpi_comm, mr_io_mpi_cart_dims, dims_file, offset_file, dims_mem)

  implicit none

  INTEGER, intent(in) :: mr_io_mpi_comm
  integer, dimension(3), intent(in) :: mr_io_mpi_cart_dims

  INTEGER, DIMENSION(3), intent(in)  :: dims_file ! Shape of full HDF5 dataset
  INTEGER, DIMENSION(3), intent(out) :: offset_file  ! Offset of data subset to read/write in HDF5 dataset
  INTEGER(HSIZE_T), DIMENSION(3), intent(out) :: dims_mem  ! Shape of data subset to read/write in feature_array

  INTEGER, DIMENSION(3) :: block_id = (/-1, -1, -1/) ! block index of this MPI process

  INTEGER     ::   mr_io_mpi_rank, mr_io_mpi_size, mr_io_mpi_error

  INTEGER     ::   i

  offset_file = (/-1, -1, -1/)
  dims_mem = (/-1, -1, -1/)

  ! Compute hyperslab offset and shape
  CALL mpi_comm_rank(mr_io_mpi_comm, mr_io_mpi_rank, mr_io_mpi_error)
  mr_io_handle_error(mr_io_mpi_error) ! FIXME: MPI error handling

  ! The IMPACT block layout is row-major
  block_id(1) = mr_io_mpi_rank/(mr_io_mpi_cart_dims(2)*mr_io_mpi_cart_dims(3))
  block_id(2) = modulo( mr_io_mpi_rank / mr_io_mpi_cart_dims(3) , mr_io_mpi_cart_dims(2) )
  block_id(3) = modulo( mr_io_mpi_rank , mr_io_mpi_cart_dims(3) )

  do i=1,3
    dims_mem(i) = (dims_file(i) + mr_io_mpi_cart_dims(i) - 1)/mr_io_mpi_cart_dims(i)
    offset_file(i) = dims_mem(i)*block_id(i)
    if ((block_id(i) + 1 == mr_io_mpi_cart_dims(i)) .and. (modulo(dims_file(i), mr_io_mpi_cart_dims(i)) /= 0)) then ! boundary correction
        dims_mem(i) =  modulo(dims_file(i), dims_mem(i))
    endif
  end do

  ! Debug-logging
  !write(stderr,*) "mr_io_parallel_spatial_hyperslap_compute (MPI rank/size, dims_file, offset_file, dims_mem)"
  !write(stderr,*) "mr_io_mpi_rank =      ",mr_io_mpi_rank
  !write(stderr,*) "mr_io_mpi_cart_dims = ",mr_io_mpi_cart_dims
  !write(stderr,*) "dims_file =           ",dims_file
  !write(stderr,*) "offset_file =         ",offset_file
  !write(stderr,*) "dim_mem =             ",dims_mem 

end subroutine mr_io_parallel_spatial_hyperslap_compute


subroutine mr_io_parallel_spatial_hyperslap_compute_padded(mr_io_mpi_comm, mr_io_mpi_cart_dims, domain_padding, &
                                                           dims_file, offset_file, dims_mem, offset_hyperslab)

  implicit none

  INTEGER, intent(in) :: mr_io_mpi_comm
  integer, dimension(3), intent(in) :: mr_io_mpi_cart_dims
  type(DomainPadding) :: domain_padding

  INTEGER, DIMENSION(3), intent(in)  :: dims_file ! Shape of full HDF5 dataset
  INTEGER, DIMENSION(3), intent(out) :: offset_file  ! Offset of data subset to read/write in HDF5 dataset
  INTEGER(HSIZE_T), DIMENSION(3), intent(out) :: dims_mem  ! Shape of data subset to read/write in feature_array
  INTEGER(HSIZE_T), DIMENSION(3), intent(out) :: offset_hyperslab  ! Offset of data within local hyperslab

  INTEGER, DIMENSION(3) :: dims_file_padded ! Shape of full HDF5 dataset

  INTEGER, DIMENSION(3) :: block_id = (/-1, -1, -1/) ! block index of this MPI process

  INTEGER     ::   mr_io_mpi_rank, mr_io_mpi_size, mr_io_mpi_error

  INTEGER     ::   i

  offset_file = (/-1, -1, -1/)
  dims_mem = (/-1, -1, -1/)
  offset_hyperslab = (/0, 0, 0/)

  ! Compute hyperslab offset and shape
  CALL mpi_comm_rank(mr_io_mpi_comm, mr_io_mpi_rank, mr_io_mpi_error)
  mr_io_handle_error(mr_io_mpi_error) ! FIXME: MPI error handling

  ! The IMPACT block layout is row-major
  block_id(1) = mr_io_mpi_rank/(mr_io_mpi_cart_dims(2)*mr_io_mpi_cart_dims(3))
  block_id(2) = modulo( mr_io_mpi_rank / mr_io_mpi_cart_dims(3) , mr_io_mpi_cart_dims(2) )
  block_id(3) = modulo( mr_io_mpi_rank , mr_io_mpi_cart_dims(3) )

  dims_file_padded = dims_file + domain_padding%lhs + domain_padding%rhs

  do i=1,3
    dims_mem(i) = (dims_file_padded(i) + mr_io_mpi_cart_dims(i) - 1)/mr_io_mpi_cart_dims(i)
    offset_file(i) = dims_mem(i)*block_id(i)
    if ((block_id(i) + 1 == mr_io_mpi_cart_dims(i)) .and. (modulo(dims_file_padded(i), mr_io_mpi_cart_dims(i)) /= 0)) then ! boundary correction
        dims_mem(i) =  modulo(dims_file_padded(i), dims_mem(i))
    endif
  end do

!  write(0,*) "mr_io_mpi_cart_rank      ", block_id
!  write(0,*) "mr_io_mpi_cart_dims      ", mr_io_mpi_cart_dims
!
!  write(0,*) "dims_file          ", dims_file
!  write(0,*) "dims_file_padded   ", dims_file_padded
!  write(0,*) "domain_padding_lhs ", domain_padding%lhs
!  write(0,*) "domain_padding_rhs ", domain_padding%rhs
!
!  write(0,*) "Extended grid hyperslab:"
!  write(0,*) "offset:            ",offset_file
!  write(0,*) "shape:             ",dims_mem

  do i=1,3
    if( (offset_file(i) < domain_padding%lhs(i)) .and. (offset_file(i) + dims_mem(i) >= domain_padding%lhs(i)) ) then
      offset_hyperslab(i) = modulo(domain_padding%lhs(i), &
                                   (dims_file_padded(i) + mr_io_mpi_cart_dims(i) - 1)/mr_io_mpi_cart_dims(i))
    end if
  end do

  dims_mem = dims_mem - max(offset_file + dims_mem - (dims_file_padded - domain_padding%rhs), (/0,0,0/))
  dims_mem = dims_mem - max(domain_padding%lhs - offset_file, (/0,0,0/))
  if(any(dims_mem <= (/0,0,0/))) then
    dims_mem = (/0,0,0/)
  end if

  offset_file = max(offset_file - domain_padding%lhs ,(/0,0,0/))

!  write(0,*) "Subdomain hyperslab:"
!  write(0,*) "offset:             ",offset_file
!  write(0,*) "local_offset:       ",offset_hyperslab
!  write(0,*) "shape:              ",dims_mem

end subroutine mr_io_parallel_spatial_hyperslap_compute_padded


INTEGER function mr_io_parallel_spatial_hyperslab_get_type(feature_global_shape, feature_offset, feature_array_shape)

  implicit none

  !real*8, dimension(:,:,:), allocatable, intent(in) :: feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(in) :: feature_array_shape  ! Shape of data subset to read/write in feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(in) :: feature_global_shape  ! Shape of data subset to read/write in feature_array
  !INTEGER, DIMENSION(3), OPTIONAL, intent(in)    :: feature_halo_shape ! Shape of halo around data to read/write in feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(in) :: feature_offset ! Offset of data subset to read/write in feature_array

  if ( all(feature_global_shape <= (/0, 0, 0/)) .and. all(feature_offset < (/0, 0, 0/)) &
       .and. all(feature_array_shape <= (/0,0,0/)) ) then
    !write(*,*) "Computing hyperslab and allocating it"
    mr_io_parallel_spatial_hyperslab_get_type = mr_io_parallel_spatial_hyperslab_enum_compute

  else if ( all(feature_global_shape > (/0, 0, 0/)) .and. &
            all(feature_offset >= (/0, 0, 0/)) .and. &
            all(feature_array_shape > (/0,0,0/)) ) then
    !write(*,*) "Reading specified hyperslab"
    if (all(feature_offset + feature_array_shape <= feature_global_shape)) then
      mr_io_parallel_spatial_hyperslab_get_type = mr_io_parallel_spatial_hyperslab_enum_specified
    else
      write(*,*) "Inconsistent arguments (feature_offset + feature_array_shape > feature_global_shape)."
      mr_io_handle_arg_error(-1)
      mr_io_parallel_spatial_hyperslab_get_type = mr_io_parallel_spatial_hyperslab_enum_error
    end if

  else
    write(*,*) "Inconsistent arguments feature array, (global) shape and (global) file offset."
    mr_io_handle_arg_error(-1)
    mr_io_parallel_spatial_hyperslab_get_type = mr_io_parallel_spatial_hyperslab_enum_error
  endif

end function mr_io_parallel_spatial_hyperslab_get_type


! ************************ DistSpatialMRI ************************

subroutine mr_io_read_parallel_spatial_feature(mr_io_mpi_comm, mr_io_mpi_cart_dims, grp_id, feature_name, &
                                               feature_array, feature_offset, feature_shape) !, feature_halo_shape)

  implicit none

  INTEGER, intent(in) :: mr_io_mpi_comm
  integer, dimension(3), intent(in) :: mr_io_mpi_cart_dims
  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name
  
  real*8, dimension(:,:,:), allocatable, intent(inout) :: feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(inout) :: feature_shape  ! Shape of data subset to read/write in feature_array
  !INTEGER, DIMENSION(3), OPTIONAL, intent(in)    :: feature_halo_shape ! Shape of halo around data to read/write in feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(inout) :: feature_offset ! Offset of data subset to read/write in feature_array

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier   
  
  INTEGER(HSIZE_T), DIMENSION(3) :: dims_file = (/-1, -1, -1/) ! Shape of full HDF5 dataset
  INTEGER(HSIZE_T), DIMENSION(3) :: max_dims_file = (/-1, -1, -1/) ! Max shape of full dataset
  
  INTEGER(HSIZE_T), DIMENSION(3) :: offset_file = (/-1, -1, -1/) ! Offset of data subset to read/write in HDF5 dataset

  INTEGER(HSIZE_T), DIMENSION(3) :: dims_mem = (/-1, -1, -1/) ! Shape of data subset to read/write in feature_array
  INTEGER(HSIZE_T), DIMENSION(3) :: offset_mem = (/-1, -1, -1/) ! Offset of data subset to read/write in feature_array

  INTEGER(HSIZE_T), DIMENSION(3) :: count_blocks = (/-1, -1, -1/) ! number of blocks to read/write
  INTEGER(HSIZE_T), DIMENSION(3) :: stride = (/-1, -1, -1/) ! stride between adjacent elements/blocks?? (unused when single block is read)

  INTEGER     ::   rank = 3       ! Dataset rank
  INTEGER     ::   error          ! Error flag
  
  INTEGER     ::   i

  !INTEGER     ::   spatial_hyperslab_type = mr_io_parallel_spatial_hyperslab_enum_error
  INTEGER, DIMENSION(3) :: feature_array_shape
  
  ! Open an existing dataset.
  CALL h5dopen_f(grp_id, feature_name, dset_id, error)
  mr_io_handle_error(error)
  
  ! Get file space
  CALL h5dget_space_f(dset_id, filespace, error)
  mr_io_handle_error(error)
  CALL h5sget_simple_extent_dims_f(filespace, dims_file, max_dims_file, error)
  mr_io_sget_simple_extent_dims_handle_error(error)


  if ( allocated(feature_array) ) then
    feature_array_shape = shape(feature_array)
  else ! feature_array not allocated - compute domain decomposition
    feature_array_shape = (/-1,-1,-1/)
  end if

  ! Includes argument check
  select case ( mr_io_parallel_spatial_hyperslab_get_type(feature_shape, feature_offset, feature_array_shape) ) !spatial_hyperslab_type )

  case ( mr_io_parallel_spatial_hyperslab_enum_specified )
    dims_mem = shape(feature_array)
    offset_file = feature_offset

    if( any(dims_file /= feature_shape) ) then
      write(*,*) "Array dimensions read from file and supplied in feature_shape disagree."
      mr_io_handle_arg_error(-1)
    endif

  case ( mr_io_parallel_spatial_hyperslab_enum_compute )
    feature_shape = dims_file

    call mr_io_parallel_spatial_hyperslap_compute(mr_io_mpi_comm, mr_io_mpi_cart_dims, feature_shape, feature_offset, dims_mem)
    offset_file = feature_offset
    ! MRI array - TODO: handle feature_halo_shape
    allocate(feature_array(dims_mem(1), dims_mem(2), dims_mem(3)))
    
  case default
    mr_io_handle_arg_error(-1)

  end select

  ! TODO: handle feature_halo_shape 
  offset_mem = (/ 0, 0, 0/)

  count_blocks = (/1,1,1/)
  stride = (/1,1,1/)

  !write (*,*) "dims_file"
  !do i=1,3
  !  write (*,*) dims_file(i)
  !end do
  !write(*,*) "dims_blocks"
  !do i=1,3
  !  write (*,*) dims_blocks(i)
  !end do
  !write (*,*) "block_id"
  !do i=1,3
  !  write (*,*) block_id(i)
  !end do
  !write (*,*) "dims_mem"
  !do i=1,3
  !  write (*,*) dims_mem(i)
  !end do
  !write (*,*) "offset_file"
  !do i=1,3
  !  write (*,*) offset_file(i)
  !end do
  
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

end subroutine mr_io_read_parallel_spatial_feature



subroutine mr_io_read_parallel_spatial(mr_io_mpi_comm, mr_io_mpi_info, mr_io_mpi_cart_dims, path, mri_inst)

  implicit none

  integer, intent (in) :: mr_io_mpi_comm, mr_io_mpi_info
  integer, dimension(3), intent(in) :: mr_io_mpi_cart_dims
  character(len=*), intent(in) :: path
  type(DistSpatialMRI), intent(out) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier

  INTEGER(HID_T) :: plist_id      ! Property list identifier   
  
  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  mr_io_handle_error(error)

!  ! Setup file access property list with parallel I/O access.
!  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
!  mr_io_handle_error(error)
!
!  CALL h5pset_fapl_mpio_f(plist_id, mr_io_mpi_comm, mr_io_mpi_info, error)
!  mr_io_handle_error(error)
!
!  ! Open existing file collectively
!  CALL h5fopen_f(trim(path), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
!  mr_io_handle_error(error)
!
!  CALL h5pclose_f(plist_id, error)
!  mr_io_handle_error(error)
  file_id = mr_io_h5_parallel_reader_open_f(mr_io_mpi_comm, mr_io_mpi_info, path)
  
  ! Open an existing group
  CALL h5gopen_f(file_id, SpatialMRI_group_name, grp_id, error)
  mr_io_handle_error(error)

  ! Read spatial feature
  CALL mr_io_read_parallel_spatial_feature(mr_io_mpi_comm, &
                                           mr_io_mpi_cart_dims, &
                                           grp_id, "scalar_feature", &
                                           mri_inst%scalar_feature%array, &
                                           mri_inst%scalar_feature%offset, &
                                           mri_inst%scalar_feature%dims)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
!  CALL h5fclose_f(file_id, error)
  error = mr_io_h5_parallel_reader_close_f(file_id)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_read_parallel_spatial


subroutine mr_io_write_parallel_spatial_feature(mr_io_mpi_comm, grp_id, feature_name, &
                                               feature_array, feature_offset, feature_shape) !, feature_halo_shape)

  implicit none

  INTEGER, intent(in) :: mr_io_mpi_comm
  INTEGER(HID_T), intent(in) :: grp_id                 ! Group identifier
  character(len=*), intent(in) :: feature_name
  
  real*8, dimension(:,:,:), allocatable, intent(in) :: feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(in) :: feature_shape  ! Shape of data subset to read/write in feature_array
  !INTEGER, DIMENSION(3), OPTIONAL, intent(in)    :: feature_halo_shape ! Shape of halo around data to read/write in feature_array
  INTEGER, DIMENSION(3), OPTIONAL, intent(in) :: feature_offset ! Offset of data subset to read/write in feature_array

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier   
  
  INTEGER(HSIZE_T), DIMENSION(3) :: dims_file = (/-1, -1, -1/) ! Shape of full HDF5 dataset  
  INTEGER(HSIZE_T), DIMENSION(3) :: offset_file = (/-1, -1, -1/) ! Offset of data subset to read/write in HDF5 dataset

  INTEGER(HSIZE_T), DIMENSION(3) :: dims_mem = (/-1, -1, -1/) ! Shape of data subset to read/write in feature_array
  INTEGER(HSIZE_T), DIMENSION(3) :: offset_mem = (/-1, -1, -1/) ! Offset of data subset to read/write in feature_array

  INTEGER(HSIZE_T), DIMENSION(3) :: count_blocks = (/-1, -1, -1/) ! number of blocks to read/write
  INTEGER(HSIZE_T), DIMENSION(3) :: stride = (/-1, -1, -1/) ! stride between adjacent elements/blocks?? (unused when single block is read)

  INTEGER     ::   rank = 3       ! Dataset rank
  INTEGER     ::   error          ! Error flag
  
  INTEGER     ::   i
  
  integer, dimension(3) :: arr
  
  ! Check arguments
  if ( any(feature_shape < shape(feature_array)) .or. any(feature_offset < (/0, 0, 0/)) ) then
    write(*,*) "Invalid arguments (hyperslab must provide sufficient space for value array)"
    mr_io_handle_arg_error(-1)
  endif

  ! TODO: handle if( present(feature_halo_shape) )
  dims_mem = shape(feature_array)
  offset_mem = (/0,0,0/)


  ! Compute dimensions  
  dims_file = feature_shape
  offset_file = feature_offset

  
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

  count_blocks = (/1,1,1/)
  stride = (/1,1,1/)
  
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

end subroutine mr_io_write_parallel_spatial_feature


subroutine mr_io_write_parallel_spatial(mr_io_mpi_comm, mr_io_mpi_info, path, mri_inst)

  implicit none

  character(len=*), intent(in) :: path
  INTEGER, intent (in) :: mr_io_mpi_comm, mr_io_mpi_info
  type(DistSpatialMRI), intent(in) :: mri_inst

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: grp_id        ! Group identifier

  INTEGER(HID_T) :: plist_id      ! Property list identifier   
  
  INTEGER     ::   error ! Error flag

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  mr_io_handle_error(error)

!  ! Create a file collectively
!  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
!  mr_io_handle_error(error)
!
!  CALL h5pset_fapl_mpio_f(plist_id, mr_io_mpi_comm, mr_io_mpi_info, error)
!  mr_io_handle_error(error)
!
!  CALL h5fcreate_f(trim(path), H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
!  mr_io_handle_error(error)
!
!  CALL h5pclose_f(plist_id, error)
!  mr_io_handle_error(error)
  file_id = mr_io_h5_parallel_writer_open_f(mr_io_mpi_comm, mr_io_mpi_info, path)

  ! Create a new group collectively  
  CALL h5pcreate_f(H5P_GROUP_CREATE_F, plist_id, error) ! FIXME: Not completely sure about flag
  mr_io_handle_error(error)

  CALL h5gcreate_f(file_id, SpatialMRI_group_name, grp_id, error, gcpl_id=plist_id)
  mr_io_handle_error(error)
  
  CALL h5pclose_f(plist_id, error)
  mr_io_handle_error(error)
  
!   print *, mri_inst%vector_feature%dims
!   print *, mri_inst%vector_feature%offset   
!   print *, shape(mri_inst%vector_feature%array)   

  ! Read spatial feature
  CALL mr_io_write_parallel_spatial_feature(mr_io_mpi_comm, &
                                           grp_id, "scalar_feature", &
                                           mri_inst%scalar_feature%array, &
                                           mri_inst%scalar_feature%offset, &
                                           mri_inst%scalar_feature%dims)

  ! Close the group
  CALL h5gclose_f(grp_id, error)
  mr_io_handle_error(error)

  ! Close the file.
!  CALL h5fclose_f(file_id, error)
  error = mr_io_h5_parallel_writer_close_f(file_id)
  mr_io_handle_error(error)

  ! Close FORTRAN interface.
  CALL h5close_f(error)
  mr_io_handle_error(error)

end subroutine mr_io_write_parallel_spatial


end module mr_io_parallel
