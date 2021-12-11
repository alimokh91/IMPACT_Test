program mr_io_test_impact_input

  USE mod_impact_core
  USE mod_vars
  USE usr_vars
  USE mr_io_parallel_spacetime
  USE mr_io_protocol

  implicit none

  integer, dimension(3) :: impact_offset = (/-1,-1,-1/)
  integer, dimension(3) :: impact_shape = (/-1,-1,-1/)

  ! For interactive testing
!  INTEGER :: mr_io_test_mpi_rank = -1
!  type(DistFlowMRI) :: mri_dest
!  integer, dimension(5) :: velocity_mean_shape = (/-1,-1,-1,-1,-1/)

  ! For interactive testing
!  INTEGER :: gdb = 0
!  do while (gdb == 0)
!      call sleep(2)
!  end do

  !--- Set configuration / topology --------------------------------------------------------------------------
  CALL configuration

  CALL impact_core_init

  ! Test domain decomposition (via coordinates, includes NB{1..3}, L{1..3} from config.txt), time not yet considered

  call find_impact_hyperslab(x1p, y1p, impact_offset(1), impact_shape(1))
  call find_impact_hyperslab(x2p, y2p, impact_offset(2), impact_shape(2))
  call find_impact_hyperslab(x3p, y3p, impact_offset(3), impact_shape(3))

  write(*,*) "***** mr_io_test_impact_input output *****"
!  write(*,*) "Offset/shape/coordinates owned by this rank (incl. halo):"
!  write(*,*) "*********************************"
!  write(*,*) impact_offset(1)
  write(*,*) impact_shape(1)
  write(*,*) x1p(1:impact_shape(1)) ! corresponds to M1
!  write(*,*) "*********************************"
!  write(*,*) impact_offset(2)
  write(*,*) impact_shape(2)
  write(*,*) x2p(1:impact_shape(2)) ! corresponds to M2
!  write(*,*) "*********************************"
!  write(*,*) impact_offset(3)
  write(*,*) impact_shape(3)
  write(*,*) x3p(1:impact_shape(3)) ! corresponds to M3
!  write(*,*) "*********************************"

! Test MRI padding parameters for Kalman filter
  write(*,*) kalman_num_data_voxels_per_process
  write(*,*) kalman_domain_padding%lhs
  write(*,*) kalman_domain_padding%rhs

! Test MRI meta information (including attributes)
  write(*,*) kalman_mri_input_file_path
  write(*,*) kalman_mri_output_file_path
  write(*,*) kalman_num_time_refinements
  write(*,*) kalman_num_spatial_refinements
  write(*,*) kalman_mri_input_attr_t_heart_cycle_period

   ! This is for interactive testing...
!  CALL MPI_COMM_RANK(MPI_COMM_WORLD,mr_io_test_mpi_rank,merror)
!  write(*,*) "MPI rank is "
!  write(*,*) mr_io_test_mpi_rank

!  call mr_io_read_parallel_flow(MPI_COMM_WORLD, MPI_INFO_NULL, kalman_mri_input_file_path, mri_dest)
!  write(*,*) "MRI: Spatial hyperslab offset"
!  write(*,*) mri_dest%velocity_mean%offset
!  write(*,*) "MRI: Spatial hyperslab shape"
!  velocity_mean_shape = shape(mri_dest%velocity_mean%array)
!  write(*,*) velocity_mean_shape(3:5)



  ! Pseudocode for implementing the Kalman filter with power of two number of MRI voxels and refinements and # MPI processes
!  num_voxel_refinements = (impact_shape - (/1,1,1/)) /shape(mri_dest%velocity_mean%array)
!  for (i,j,k)=(1,1,1):shape(mri_dest%velocity_mean%array):
!     for (l,m,n)=(1,1,1):num_voxel_refinements:
!         kalman_filter(mri_dest%velocity_mean%array(i,j,k), mri_dest%velocity_cov%array(i,j,k), impact_velocity_datastructure( num_voxel_refinements*(i,j,k)+(l,m,n) ) )
!     end for
!  end for

  !...

  ! TODO: n_timesteps

  !#--- Extents --------------------------------------------------------------
  !time_start      {{ time_start }}
  !time_end        {{ time_end }}

  CALL impact_core_finalize


  contains
    subroutine find_impact_hyperslab(xp, yp, offset, ishape)

      implicit none

      real*8, dimension(:), intent(in) :: xp
      real*8, dimension(:), intent(in) :: yp
      INTEGER, intent(out) :: offset
      INTEGER, intent(out) :: ishape

      INTEGER, dimension(1) :: tmp
      INTEGER :: lb_x
      INTEGER :: i_y
      INTEGER :: ub_y

      tmp = lbound(yp)-1; offset = tmp(1)
      tmp = lbound(xp); lb_x = tmp(1) - b1L + 1
      tmp = lbound(yp); i_y = tmp(1)
      tmp = ubound(yp); ub_y = tmp(1)
      tmp = shape(xp) - (/b1U - b1L +1/); ishape = tmp(1)

      do while (i_y <= ub_y)
        if (yp(i_y) < xp(lb_x)) then
          i_y = i_y+1
        else
          offset = i_y
          return
        endif
      end do

      call backtrace
      call abort

    end subroutine find_impact_hyperslab

end program mr_io_test_impact_input
