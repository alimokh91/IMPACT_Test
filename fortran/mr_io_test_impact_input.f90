program mr_io_test_impact_input

  USE mod_vars
  USE usr_vars
  USE mr_io_parallel_spacetime
  USE mr_protocol

  implicit none

  integer, dimension(3) :: impact_offset = (/-1,-1,-1/)
  integer, dimension(3) :: impact_shape = (/-1,-1,-1/)

  ! For interactive testing
!  INTEGER :: mpi_rank = -1
!  character(len=100) :: path = "mr_io_test_impact_input.h5"
!  type(DistHPCPredictMRI) :: mri_dest
!  integer, dimension(5) :: velocity_mean_shape = (/-1,-1,-1,-1,-1/)

  ! For interactive testing
!  INTEGER :: gdb = 0
!  do while (gdb == 0)
!      call sleep(2)
!  end do

  call impact_initialization_code

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


   ! This is for interactive testing...
!  CALL MPI_COMM_RANK(MPI_COMM_WORLD,mpi_rank,merror)
!  write(*,*) "MPI rank is "
!  write(*,*) mpi_rank

!  call mr_io_read_parallel_hpcpredict(MPI_COMM_WORLD, MPI_INFO_NULL, path, mri_dest)
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

  call impact_finalization_code


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




! Directly copied from impact.f90
subroutine impact_initialization_code

  USE mod_dims
  USE mod_vars
  USE mod_setup
  USE mod_coeffs
  USE mod_lib
  USE mod_inout
  USE mod_test
  USE mod_geometry
  USE mod_timeint
  USE HDF5
  USE usr_func

  IMPLICIT NONE

  INCLUDE 'mpif.h'


  !----------------------------------------------------------------------------------------------------------!
  ! Annotations / TO-DO list:                                                                                !
  !  - henniger 2011: unite the four init_* (general, parallel, boundaries, limits), already intermingled    !
  !  - mjohn 060513: moved task 3 (eigensolver) to a usr file because mod functions depended themselves on   !
  !                  usr functions, i.e. variables / fields of task 3 were not independent! May be recovered !
  !                  as part of a mod file if new fields are provided and dependence on usr functions lifted !
  !----------------------------------------------------------------------------------------------------------!


  !===========================================================================================================
  !=== Initialization ========================================================================================
  !===========================================================================================================
  !--- Initialize MPI ----------------------------------------------------------------------------------------
  CALL MPI_INIT(merror)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,merror)

  !--- Initialize HDF5 ---------------------------------------------------------------------------------------
  CALL h5open_f(herror)

  !--- Set alarm if queue is being used ----------------------------------------------------------------------
  CALL init_alarm

  !--- Set configuration / topology --------------------------------------------------------------------------
  CALL configuration

  !--- Test of input parameters ------------------------------------------------------------------------------
  CALL test_parameter

   !--- General -------------------------------------------------------------------------------------------
   CALL init_general
   !--- MPI -----------------------------------------------------------------------------------------------
   CALL init_parallel
   !--- Type of boundary conditions -----------------------------------------------------------------------
   CALL init_boundaries
   !--- Limits of indices ---------------------------------------------------------------------------------
   CALL init_limits

  !--- Physical coordinates ----------------------------------------------------------------------------------
  CALL coordinates

  IF (task .EQ. 5) THEN
  write(*,*) 'x1p(s1p) =', x1p(s1p), ', x1p(n1p) =', x1p(n1p), ', s1p = ', s1p,', n1p = ', n1p
  write(*,*) 'y1p(1) =', y1p(1), ', y1p(m1) =', y1p(m1), ', m1 = ', m1
  write(*,*) 's3p = ', s3p,', n3p = ', n3p
  write(*,*) 's33 = ', s33,', n33 = ', n33
  write(*,*) 's33b = ', s33b,', n33b = ', n33b
  write(*,*) 'n3 = ', n3
  write(*,*) 'x3w(s33) = ', x3w(s33),', x3w(n33) = ', x3w(n33)
  write(*,*) 'x3w(s33b) = ', x3w(s33b+1),', x3w(n33b) = ', x3w(n33b)
  write(*,*) 'm1 = ', m1,'m2 = ', m2 ,'m3 = ', m3
  write(*,*) 'ib(1,1) = ', iB(1,1), ' n1 = ', n1 ,';'
  write(*,*) 'ib(2,1) = ', iB(2,1), ' n2 = ', n2 ,';'
  write(*,*) 'ib(3,1) = ', iB(3,1), ' n3 = ', n3 ,';'


  END IF
  !--- reference velocity ------------------------------------------------------------------------------------
  U_ref = (mu_fluid*Re)/(rho_fluid*L_ref)

  !--- Determine differential coefficients -------------------------------------------------------------------
  CALL FD_coeffs

  !--- Get stencil of operator div(grad( )), order of convergence is 2 ---------------------------------------
  CALL get_stencil
  CALL get_stencil_Helm

  !--- Get interpolation coefficients (multigrid, order of convergence is 2) ---------------------------------
  CALL interp_coeffs
  CALL interp_coeffs_Helm
  CALL restr_coeffs
  CALL restr_coeffs_Helm

  !--- Differenzenkoeffizienten testen -----------------------------------------------------------------------
  ! bbecsek 031214: TODO: according to usr_config these routines are deprecated.
  IF (write_test_yes) CALL test_coeffs

  !--- Weights for stop criterion ----------------------------------------------------------------------------
  CALL get_weights

  !--- Beta (for analysis) -----------------------------------------------------------------------------------
  CALL get_beta
  !===========================================================================================================

end subroutine impact_initialization_code

! Directly copied from impact.f90
subroutine impact_finalization_code

  USE mod_dims
  USE mod_vars
  USE mod_setup
  USE mod_coeffs
  USE mod_lib
  USE mod_inout
  USE mod_test
  USE mod_geometry
  USE mod_timeint
  USE HDF5
  USE usr_func

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  !---- close HDF5 -------------------------------------------------------------------------------------------
  CALL h5close_f(herror)

  !---- close MPI --------------------------------------------------------------------------------------------
  CALL MPI_FINALIZE(merror)
  !===========================================================================================================

end subroutine impact_finalization_code

