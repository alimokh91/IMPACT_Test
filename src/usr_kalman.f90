!***********************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Hennniger, Institute of Fluid Dynamics, ETH Zurich (hennniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* created by Dario De Marinis, ARTORG CVE, University of Bern (dario.demarinis@artorg.unibe.ch)             *
!* January 2018                                                                                              *
!*************************************************************************************************************

! file usr_kalman.f90
! file containing subroutines used for data assimilation.

! !pgi$g unroll = n:8
! !pgi$r unroll = n:8
! !pgi$l unroll = n:8

  SUBROUTINE validate_input_kalman()
 
  USE mod_vars, only: refine_2dir_yes
  IMPLICIT NONE

  IF(refine_2dir_yes .NEQV. .FALSE.) THEN
    write(0,*) "refine_2dir_yes = TRUE causes 2D-coordinates that are not aligned with MRI voxels. Exiting..."
    CALL backtrace
    CALL flush()
    CALL abort
  END IF  
  
  END SUBROUTINE validate_input_kalman



  SUBROUTINE open_kalman
  ! (basic subroutine)

  USE mod_dims
  USE mod_vars
  USE mod_lib
  USE mod_exchange
  USE mod_diff
  USE mod_inout
  USE usr_func
  USE usr_vars
  USE MPI
  USE HDF5
  USE mr_io_protocol
  USE mr_io_parallel_spacetime 

  IMPLICIT NONE

  TYPE(kalman_t), pointer :: klmn
  INTEGER :: i,ii,j,jj,k,kk,n,m,pts
  INTEGER, DIMENSION(2,3) :: bounds, bounds2
  CHARACTER(LEN=2) ::  phs

  if (rank.eq.0) write(0,*) "open kalman... "

  call validate_input_kalman() 

  write_kalm_count = 0

  ! vel(:,:,:,i) --> worki(:,:,:)
  IF (task == 1) CALL interpolate_vel(.TRUE.)
 
  !===========================================================================================================
  !=== read the observed data (components, time, x, y, z) ====================================================
  !===========================================================================================================
!  if (associated(mri_inst).eqv..false.) then
     nullify(mri_inst); allocate(mri_inst)
     mri_inst%domain_padding = kalman_domain_padding
     CALL mr_io_read_parallel_segmentedflow_padded(MPI_COMM_WORLD, MPI_INFO_NULL, (/NB1,NB2,NB3/), &
          trim(kalman_mri_input_file_path), mri_inst)
     CALL h5open_f(herror) ! Required as hpc-predict-io closes HDF5 environment with h5close_f
!  end if
  mri_inst%mri%velocity_mean%array = mri_inst%mri%velocity_mean%array/U_ref
  mri_inst%mri%velocity_cov%array  = mri_inst%mri%velocity_cov%array/U_ref/U_ref
  !===========================================================================================================

  ALLOCATE(dtime_kalm_phases(1:intervals)); dtime_kalm_phases = 0.
  dtime_kalm_phases(1:intervals-1) = mri_inst%mri%t_coordinates(2:intervals) - mri_inst%mri%t_coordinates(1:intervals-1)
  dtime_kalm_phases(intervals) = kalman_mri_input_attr_t_heart_cycle_period - &
                                 mri_inst%mri%t_coordinates(intervals) + mri_inst%mri%t_coordinates(1)
  !===========================================================================================================

  !===========================================================================================================
  !=== construct matrices x_f, H, R, P_f =====================================================================
  !===========================================================================================================
  bounds(1,1) = lbound(mri_inst%mri%velocity_mean%array,3)
  bounds(2,1) = ubound(mri_inst%mri%velocity_mean%array,3)
  bounds(1,2) = lbound(mri_inst%mri%velocity_mean%array,4)
  bounds(2,2) = ubound(mri_inst%mri%velocity_mean%array,4)
  bounds(1,3) = lbound(mri_inst%mri%velocity_mean%array,5)
  bounds(2,3) = ubound(mri_inst%mri%velocity_mean%array,5)

  write(*,*) 'sonoqui'

  if (size(mri_inst%mri%velocity_mean%array,3)*size(mri_inst%mri%velocity_mean%array,4)*size(mri_inst%mri%velocity_mean%array,5).ne.0) then
     allocate(kalman_first)
     klmn => kalman_first

     klmn%n =         size(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)+1
     klmn%n = klmn%n*(size(mri_inst%mri%velocity_mean%array,4)*kalman_num_spatial_refinements(2)+1)
     klmn%n = klmn%n*(size(mri_inst%mri%velocity_mean%array,5)*kalman_num_spatial_refinements(3)+1)

     klmn%m =        size(mri_inst%mri%velocity_mean%array,3)
     klmn%m = klmn%m*size(mri_inst%mri%velocity_mean%array,4)
     klmn%m = klmn%m*size(mri_inst%mri%velocity_mean%array,5)

     NULLIFY(klmn%K);  ALLOCATE(klmn%K(1:3*klmn%n,1)); klmn%K = 0.0

     pts = (kalman_num_spatial_refinements(1)+1)*(kalman_num_spatial_refinements(2)+1)*(kalman_num_spatial_refinements(3)+1)

     klmn%m = 1
     klmn%n = pts

     NULLIFY(klmn%obs_data);  ALLOCATE(klmn%obs_data(1:3*klmn%m,1)); klmn%obs_data = 0.0
     NULLIFY(klmn%obs_covar); ALLOCATE(klmn%obs_covar(1:3*klmn%m,1:3*klmn%m)); klmn%obs_covar = 0.0
     NULLIFY(klmn%K_vx);      ALLOCATE(klmn%K_vx(1:3*klmn%n,1)); klmn%K_vx = 0.0
     NULLIFY(klmn%obs_oper);  ALLOCATE(klmn%obs_oper(1:3*klmn%m,1:3*klmn%n)); klmn%obs_oper = 0.0
     NULLIFY(klmn%PfHt);      ALLOCATE(klmn%PfHt(1:3*klmn%n,1:3*klmn%m)); klmn%PfHt = 0.0
 
     m = 0
     n = 0
     do n = 0,klmn%n-1
       klmn%obs_oper(3*m+1,3*n+1) = 1.0/pts
       klmn%obs_oper(3*m+2,3*n+2) = 1.0/pts
       klmn%obs_oper(3*m+3,3*n+3) = 1.0/pts
     end do

     bounds2(1,1) = (bounds(1,1)-1)*kalman_num_spatial_refinements(1)+1
     bounds2(2,1) =  bounds(2,1)   *kalman_num_spatial_refinements(1)+1
     bounds2(1,2) = (bounds(1,2)-1)*kalman_num_spatial_refinements(2)+1
     bounds2(2,2) =  bounds(2,2)   *kalman_num_spatial_refinements(2)+1
     bounds2(1,3) = (bounds(1,3)-1)*kalman_num_spatial_refinements(3)+1
     bounds2(2,3) =  bounds(2,3)   *kalman_num_spatial_refinements(3)+1

     allocate(wgt_interp(bounds2(1,1):bounds2(2,1)+1,bounds2(1,2):bounds2(2,2)+1,bounds2(1,3):bounds2(2,3)+1))
     wgt_interp = 0
     DO k = bounds(1,3), bounds(2,3)
       DO j = bounds(1,2), bounds(2,2)
         DO i = bounds(1,1), bounds(2,1)
           DO kk = (k-1)*kalman_num_spatial_refinements(3)+1,k*kalman_num_spatial_refinements(3)+1
             DO jj = (j-1)*kalman_num_spatial_refinements(2)+1,j*kalman_num_spatial_refinements(2)+1
               DO ii = (i-1)*kalman_num_spatial_refinements(1)+1,i*kalman_num_spatial_refinements(1)+1
                 wgt_interp(ii,jj,kk) = wgt_interp(ii,jj,kk) + 1
               END DO
             END DO
           END DO
         END DO
       END DO
     END DO

     ALLOCATE(mean_f (1:intervals,bounds2(1,1):bounds2(2,1),bounds2(1,2):bounds2(2,2),bounds2(1,3):bounds2(2,3),1:3))
     ALLOCATE(covar_f(1:intervals,bounds2(1,1):bounds2(2,1),bounds2(1,2):bounds2(2,2),bounds2(1,3):bounds2(2,3),1:6))
     mean_f (1:intervals,bounds2(1,1):bounds2(2,1),bounds2(1,2):bounds2(2,2),bounds2(1,3):bounds2(2,3),1:3) = 0.0
     covar_f(1:intervals,bounds2(1,1):bounds2(2,1),bounds2(1,2):bounds2(2,2),bounds2(1,3):bounds2(2,3),1:6) = 0.0

     DO k = bounds(1,3), bounds(2,3)
       DO j = bounds(1,2), bounds(2,2)
         DO i = bounds(1,1), bounds(2,1)
           DO kk = (k-1)*kalman_num_spatial_refinements(3)+1,k*kalman_num_spatial_refinements(3)+1
             DO jj = (j-1)*kalman_num_spatial_refinements(2)+1,j*kalman_num_spatial_refinements(2)+1
               DO ii = (i-1)*kalman_num_spatial_refinements(1)+1,i*kalman_num_spatial_refinements(1)+1
                 work1(ii,jj,kk) = mri_inst%mri%velocity_mean%array(1,1,i,j,k)
                 work2(ii,jj,kk) = mri_inst%mri%velocity_mean%array(2,1,i,j,k)
                 work3(ii,jj,kk) = mri_inst%mri%velocity_mean%array(3,1,i,j,k)
               END DO
             END DO
           END DO
         END DO
       END DO
     END DO

     NULLIFY(klmn%next)
  end if
  !===========================================================================================================

  !-----------------------------------------------------------------------------------------------------------
  !--- worki(:,:,:) --> vel(:,:,:,i) -------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  CALL interpolate2_pre_vel(.TRUE.,1,work1,vel(:,:,:,1))
  CALL interpolate2_pre_vel(.TRUE.,2,work2,vel(:,:,:,2))
  CALL interpolate2_pre_vel(.TRUE.,3,work3,vel(:,:,:,3))
  ! exchange across boundaries on vel grid
  CALL exchange_all_all(.TRUE.,vel)

  RETURN

  END SUBROUTINE open_kalman


  SUBROUTINE close_kalman
  ! (basic subroutine)

  USE mod_dims
  USE mod_vars
  USE usr_func
  USE usr_vars
  USE MPI
  USE mr_io_protocol
  USE mr_io_parallel_spacetime

  IMPLICIT NONE

  TYPE(kalman_t), pointer :: klmn

  DEALLOCATE(dtime_kalm_phases)

  if (associated(kalman_first)) then 
     DEALLOCATE(wgt_interp)
     DEALLOCATE(mean_f,covar_f)
     klmn => kalman_first
     CALL destroy_kalman(klmn)
     DEALLOCATE(kalman_first)
  end if

  call mr_io_deallocate_dist_segmentedflow_mri_padded(mri_inst)

  contains 
  RECURSIVE SUBROUTINE destroy_kalman(klmn)

     implicit none

     TYPE(kalman_t), pointer, intent(inout) :: klmn

     ! ---------------------------------

     if (associated(klmn%next)) then
        call destroy_kalman(klmn%next)
        DEALLOCATE(klmn%next)
     end if

     DEALLOCATE(klmn%obs_data,klmn%obs_covar,klmn%obs_oper)
     DEALLOCATE(klmn%K,klmn%K_vx)
     DEALLOCATE(klmn%PfHt)

  END SUBROUTINE destroy_kalman

  END SUBROUTINE close_kalman


  SUBROUTINE compute_kalman
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE mod_lib
  USE mod_exchange
  USE mod_diff
  USE mod_inout
  USE usr_func
  USE usr_vars
  USE MPI
  USE HDF5

  IMPLICIT NONE
  
  TYPE(kalman_t), pointer :: klmn
  INTEGER :: i, j, k, ii, jj, kk, iii,jjj,kkk, n, nnn, m, mm
  INTEGER, DIMENSION(2,3) :: bounds
  CHARACTER(LEN=8) :: count_char
  CHARACTER(LEN=2) :: count_char2

  phase = mod(write_kalm_count,intervals) + 1

  IF (rank == 0) WRITE(*,'(a,i8,a,i8,a)') 'kalman repetition', write_kalm_count/intervals+1, '   for phase', phase,' ...'

  bounds(1,1) = lbound(mri_inst%mri%velocity_mean%array,3)
  bounds(2,1) = ubound(mri_inst%mri%velocity_mean%array,3)
  bounds(1,2) = lbound(mri_inst%mri%velocity_mean%array,4)
  bounds(2,2) = ubound(mri_inst%mri%velocity_mean%array,4)
  bounds(1,3) = lbound(mri_inst%mri%velocity_mean%array,5)
  bounds(2,3) = ubound(mri_inst%mri%velocity_mean%array,5)
 
  !===========================================================================================================
  !=== interpolation of velocity to pressure grid points =====================================================
  !===========================================================================================================
  ! vel(:,:,:,i) --> worki(:,:,:)
  IF (task == 1) CALL interpolate_vel(.TRUE.)
  !===========================================================================================================

!  if (phase.eq.1) then
!     CALL num_to_string(8,write_kalm_count/intervals+1,count_char)
!   
!     if (associated(mri_inst).eqv..true.) call mr_io_deallocate_dist_segmentedflow_mri_padded(mri_inst)
!     nullify(mri_inst); allocate(mri_inst)
!     mri_inst%domain_padding = kalman_domain_padding
!
!     CALL mr_io_read_parallel_segmentedflow_padded(MPI_COMM_WORLD, MPI_INFO_NULL, (/NB1,NB2,NB3/), &
!          trim(kalman_mri_input_file_path)//'_cycle.'//count_char//'.h5', mri_inst)
!     CALL h5open_f(herror) ! Required as hpc-predict-io closes HDF5 environment with h5close_f
!
!     mri_inst%mri%velocity_mean%array = mri_inst%mri%velocity_mean%array/U_ref
!     mri_inst%mri%velocity_cov%array  = mri_inst%mri%velocity_cov%array/U_ref/U_ref
!  end if

  klmn => kalman_first
  do while(associated(klmn))

    !========================================================================================================
    !=== forecast statistics ================================================================================
    !========================================================================================================
    !--------------------------------------------------------------------------------------------------------
    !--- <u_i u_j>(t_k) = <u'_i u'_j>(t_k) + <u_i>(t_k) <u_j>(t_k) ------------------------------------------
    !--------------------------------------------------------------------------------------------------------
    if (write_kalm_count/intervals.ge.1) then
      covar_f(phase,:,:,:,1) = covar_f(phase,:,:,:,1) + mean_f(phase,:,:,:,1)*mean_f(phase,:,:,:,1) &
                                                      - 1.0e1/(write_kalm_count/intervals)
      covar_f(phase,:,:,:,2) = covar_f(phase,:,:,:,2) + mean_f(phase,:,:,:,2)*mean_f(phase,:,:,:,2) &
                                                      - 1.0e1/(write_kalm_count/intervals)
      covar_f(phase,:,:,:,3) = covar_f(phase,:,:,:,3) + mean_f(phase,:,:,:,3)*mean_f(phase,:,:,:,3) &
                                                      - 1.0e1/(write_kalm_count/intervals)
      covar_f(phase,:,:,:,4) = covar_f(phase,:,:,:,4) + mean_f(phase,:,:,:,1)*mean_f(phase,:,:,:,2)
      covar_f(phase,:,:,:,5) = covar_f(phase,:,:,:,5) + mean_f(phase,:,:,:,1)*mean_f(phase,:,:,:,3)
      covar_f(phase,:,:,:,6) = covar_f(phase,:,:,:,6) + mean_f(phase,:,:,:,2)*mean_f(phase,:,:,:,3)
    end if
    !--------------------------------------------------------------------------------------------------------
    !--- <u_i>(t_k) -----------------------------------------------------------------------------------------
    !--- <u_i u_j>(t_k) -------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------
     mean_f(phase,:,:,:,:) =  mean_f(phase,:,:,:,:)*( write_kalm_count/intervals )
    covar_f(phase,:,:,:,:) = covar_f(phase,:,:,:,:)*( write_kalm_count/intervals )
    DO k = (bounds(1,3)-1)*kalman_num_spatial_refinements(3)+1,bounds(2,3)*kalman_num_spatial_refinements(3)+1
       DO j = (bounds(1,2)-1)*kalman_num_spatial_refinements(2)+1,bounds(2,2)*kalman_num_spatial_refinements(2)+1
          DO i = (bounds(1,1)-1)*kalman_num_spatial_refinements(1)+1,bounds(2,1)*kalman_num_spatial_refinements(1)+1
              mean_f(phase,i,j,k,1) =  mean_f(phase,i,j,k,1) + work1(i,j,k)
              mean_f(phase,i,j,k,2) =  mean_f(phase,i,j,k,2) + work2(i,j,k)
              mean_f(phase,i,j,k,3) =  mean_f(phase,i,j,k,3) + work3(i,j,k)
             covar_f(phase,i,j,k,1) = covar_f(phase,i,j,k,1) + work1(i,j,k)*work1(i,j,k)
             covar_f(phase,i,j,k,2) = covar_f(phase,i,j,k,2) + work2(i,j,k)*work2(i,j,k)
             covar_f(phase,i,j,k,3) = covar_f(phase,i,j,k,3) + work3(i,j,k)*work3(i,j,k)
             covar_f(phase,i,j,k,4) = covar_f(phase,i,j,k,4) + work1(i,j,k)*work2(i,j,k)
             covar_f(phase,i,j,k,5) = covar_f(phase,i,j,k,5) + work1(i,j,k)*work3(i,j,k)
             covar_f(phase,i,j,k,6) = covar_f(phase,i,j,k,6) + work2(i,j,k)*work3(i,j,k)
          END DO
       END DO
    END DO
     mean_f(phase,:,:,:,:) =  mean_f(phase,:,:,:,:)/( write_kalm_count/intervals + 1 )
    covar_f(phase,:,:,:,:) = covar_f(phase,:,:,:,:)/( write_kalm_count/intervals + 1 )
    !--------------------------------------------------------------------------------------------------------
    !--- <u'_i u'_j>(t_k) -----------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------
    covar_f(phase,:,:,:,1) = covar_f(phase,:,:,:,1) - mean_f(phase,:,:,:,1)*mean_f(phase,:,:,:,1) &
                                                    + 1.0e1/(write_kalm_count/intervals+1)
    covar_f(phase,:,:,:,2) = covar_f(phase,:,:,:,2) - mean_f(phase,:,:,:,2)*mean_f(phase,:,:,:,2) &
                                                    + 1.0e1/(write_kalm_count/intervals+1)
    covar_f(phase,:,:,:,3) = covar_f(phase,:,:,:,3) - mean_f(phase,:,:,:,3)*mean_f(phase,:,:,:,3) &
                                                    + 1.0e1/(write_kalm_count/intervals+1)
    covar_f(phase,:,:,:,4) = covar_f(phase,:,:,:,4) - mean_f(phase,:,:,:,1)*mean_f(phase,:,:,:,2)
    covar_f(phase,:,:,:,5) = covar_f(phase,:,:,:,5) - mean_f(phase,:,:,:,1)*mean_f(phase,:,:,:,3)
    covar_f(phase,:,:,:,6) = covar_f(phase,:,:,:,6) - mean_f(phase,:,:,:,2)*mean_f(phase,:,:,:,3)
    !========================================================================================================

    klmn => klmn%next
  end do

  !==========================================================================================================
  !=== apply kalman filtering ===============================================================================
  !==========================================================================================================
  klmn => kalman_first
  do while(associated(klmn))
    !----------------------------------------------------------------------------------------------------------
    !--- Fill the matrices and the vectors involved in the kalman filter algorithm ----------------------------
    !----------------------------------------------------------------------------------------------------------
    klmn%K = 0.0

    DO k = bounds(1,3), bounds(2,3)
      DO j = bounds(1,2), bounds(2,2)
        DO i = bounds(1,1), bounds(2,1)
 
          m = 0
          klmn%obs_covar(3*m+1:3*m+3,3*m+1:3*m+3) = mri_inst%mri%velocity_cov%array(1:3,1:3,phase,i,j,k)
          if (mri_inst%mri%segmentation_prob%array(phase,i,j,k).le.0.5) then
             klmn%obs_covar(3*m+1:3*m+3,3*m+1:3*m+3) = 1.0e6+klmn%obs_covar(3*m+1:3*m+3,3*m+1:3*m+3)
          end if
          klmn%obs_data (3*m+1:3*m+3,1) = mri_inst%mri%velocity_mean%array (1:3,phase,i,j,k)

          n = 0
          klmn%PfHt = 0.0
          DO kk = (k-1)*kalman_num_spatial_refinements(3)+1,k*kalman_num_spatial_refinements(3)+1
            DO jj = (j-1)*kalman_num_spatial_refinements(2)+1,j*kalman_num_spatial_refinements(2)+1
              DO ii = (i-1)*kalman_num_spatial_refinements(1)+1,i*kalman_num_spatial_refinements(1)+1

          !----------------------------------------------------------------------------------------------
          !--- matrix R ---------------------------------------------------------------------------------
          !--- data = data - H*mu_f ---------------------------------------------------------------------
          !----------------------------------------------------------------------------------------------
                klmn%obs_data(3*m+1,1) = klmn%obs_data(3*m+1,1) - klmn%obs_oper(3*m+1,3*n+1)*work1(ii,jj,kk)
                klmn%obs_data(3*m+2,1) = klmn%obs_data(3*m+2,1) - klmn%obs_oper(3*m+2,3*n+2)*work2(ii,jj,kk)
                klmn%obs_data(3*m+3,1) = klmn%obs_data(3*m+3,1) - klmn%obs_oper(3*m+3,3*n+3)*work3(ii,jj,kk)
          !----------------------------------------------------------------------------------------------

          !----------------------------------------------------------------------------------------------
          !--- klmn%PfHt = Pf*transpose(H) valid for three-diagonal Pf ----------------------------------
          !----------------------------------------------------------------------------------------------
                DO mm = 1,3*klmn%m
                  klmn%PfHt(3*n+1,mm) = covar_f(phase,ii,jj,kk,1)*klmn%obs_oper(mm,3*n+1) + &
                                        covar_f(phase,ii,jj,kk,4)*klmn%obs_oper(mm,3*n+2) + &
                                        covar_f(phase,ii,jj,kk,5)*klmn%obs_oper(mm,3*n+3)

                  klmn%PfHt(3*n+2,mm) = covar_f(phase,ii,jj,kk,4)*klmn%obs_oper(mm,3*n+1) + &
                                        covar_f(phase,ii,jj,kk,2)*klmn%obs_oper(mm,3*n+2) + &
                                        covar_f(phase,ii,jj,kk,6)*klmn%obs_oper(mm,3*n+3)

                  klmn%PfHt(3*n+3,mm) = covar_f(phase,ii,jj,kk,5)*klmn%obs_oper(mm,3*n+1) + &
                                        covar_f(phase,ii,jj,kk,6)*klmn%obs_oper(mm,3*n+2) + &
                                        covar_f(phase,ii,jj,kk,3)*klmn%obs_oper(mm,3*n+3)
                END DO

                n = n + 1
              END DO
            END DO
          END DO
          !----------------------------------------------------------------------------------------------

          !----------------------------------------------------------------------------------------------
          !--- R = H*Pf*transpose(H) + R (resulting matrix is symmetric) --------------------------------
          !----------------------------------------------------------------------------------------------
          n = 0
          DO kk = (k-1)*kalman_num_spatial_refinements(3)+1,k*kalman_num_spatial_refinements(3)+1
            DO jj = (j-1)*kalman_num_spatial_refinements(2)+1,j*kalman_num_spatial_refinements(2)+1
              DO ii = (i-1)*kalman_num_spatial_refinements(1)+1,i*kalman_num_spatial_refinements(1)+1

                DO mm = 3*m+1,3*klmn%m
                  klmn%obs_covar(3*m+1,mm) = klmn%obs_covar(3*m+1,mm) + klmn%obs_oper(3*m+1,3*n+1)*klmn%PfHt(3*n+1,mm)
                  klmn%obs_covar(3*m+2,mm) = klmn%obs_covar(3*m+2,mm) + klmn%obs_oper(3*m+2,3*n+2)*klmn%PfHt(3*n+2,mm)
                  klmn%obs_covar(3*m+3,mm) = klmn%obs_covar(3*m+3,mm) + klmn%obs_oper(3*m+3,3*n+3)*klmn%PfHt(3*n+3,mm)
                END DO

                n = n + 1
              END DO
            END DO
          END DO

          DO mm = 1, 3*klmn%m-1
            klmn%obs_covar(mm+1:3*klmn%m,mm) = klmn%obs_covar(mm,mm+1:3*klmn%m)
          END DO
          !----------------------------------------------------------------------------------------------

          !----------------------------------------------------------------------------------------------
          !--- compute Kalman gain K = Pf*transpose(H)*(H*Pf*transpose(H)+R)^-1*(data - H*mu_f) ---------
          !----------------------------------------------------------------------------------------------
          klmn%K_vx = 0.0
          CALL kalman_gain(klmn%PfHt,klmn%obs_covar,klmn%obs_data,klmn%K_vx,size(klmn%obs_oper,2),size(klmn%obs_oper,1))
          !----------------------------------------------------------------------------------------------

          !----------------------------------------------------------------------------------------------
          !--- then add K_vx to K -----------------------------------------------------------------------
          !----------------------------------------------------------------------------------------------
          n = 0
          DO kk = (k-1)*kalman_num_spatial_refinements(3)+1,k*kalman_num_spatial_refinements(3)+1
            DO jj = (j-1)*kalman_num_spatial_refinements(2)+1,j*kalman_num_spatial_refinements(2)+1
              DO ii = (i-1)*kalman_num_spatial_refinements(1)+1,i*kalman_num_spatial_refinements(1)+1
                nnn =      ii-((bounds(1,1)-1)*kalman_num_spatial_refinements(1)+1)
                nnn = nnn + (jj-((bounds(1,2)-1)*kalman_num_spatial_refinements(2)+1))* &
                        (size(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)+1)
                nnn = nnn + (kk-((bounds(1,3)-1)*kalman_num_spatial_refinements(3)+1))* &
                        (size(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)+1) * &
                        (size(mri_inst%mri%velocity_mean%array,4)*kalman_num_spatial_refinements(2)+1)

                klmn%K(3*nnn+1,1) = klmn%K(3*nnn+1,1) + klmn%K_vx(3*n+1,1)/wgt_interp(ii,jj,kk)
                klmn%K(3*nnn+2,1) = klmn%K(3*nnn+2,1) + klmn%K_vx(3*n+2,1)/wgt_interp(ii,jj,kk)
                klmn%K(3*nnn+3,1) = klmn%K(3*nnn+3,1) + klmn%K_vx(3*n+3,1)/wgt_interp(ii,jj,kk)

                n = n + 1
              END DO
            END DO
          END DO
         !----------------------------------------------------------------------------------------------

        END DO
      END DO
    END DO

    !----------------------------------------------------------------------------------------------
    !--- then update u_a = u_f + K ----------------------------------------------------------------
    !----------------------------------------------------------------------------------------------
    DO kk = (bounds(1,3)-1)*kalman_num_spatial_refinements(3)+1,bounds(2,3)*kalman_num_spatial_refinements(3)+1
      DO jj = (bounds(1,2)-1)*kalman_num_spatial_refinements(2)+1,bounds(2,2)*kalman_num_spatial_refinements(2)+1
        DO ii = (bounds(1,1)-1)*kalman_num_spatial_refinements(1)+1,bounds(2,1)*kalman_num_spatial_refinements(1)+1
          n =      ii-((bounds(1,1)-1)*kalman_num_spatial_refinements(1)+1)
          n = n + (jj-((bounds(1,2)-1)*kalman_num_spatial_refinements(2)+1))* &
                  (size(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)+1)
          n = n + (kk-((bounds(1,3)-1)*kalman_num_spatial_refinements(3)+1))* &
                  (size(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)+1) * &
                  (size(mri_inst%mri%velocity_mean%array,4)*kalman_num_spatial_refinements(2)+1)

         work1(ii,jj,kk) = work1(ii,jj,kk) + klmn%K(3*n+1,1)
         work2(ii,jj,kk) = work2(ii,jj,kk) + klmn%K(3*n+2,1)
         work3(ii,jj,kk) = work3(ii,jj,kk) + klmn%K(3*n+3,1)

        END DO
      END DO
    END DO

    klmn => klmn%next
  end do
  !===========================================================================================================

  !----------------------------------------------------------------------------------------------
  !--- worki(:,:,:) --> vel(:,:,:,i) -------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------
  CALL interpolate2_pre_vel(.TRUE.,1,work1,vel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1))
  CALL interpolate2_pre_vel(.TRUE.,2,work2,vel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),2))
  CALL interpolate2_pre_vel(.TRUE.,3,work3,vel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),3))
   ! exchange across boundaries on vel grid
  CALL exchange_all_all(.TRUE.,vel)
  !===========================================================================================================

  write_out_kalm = .FALSE.
  write_kalm_count = write_kalm_count + 1
  dtime_out_kalm = dtime_kalm_phases(phase)
  time_out_kalm = time_out_kalm + dtime_out_kalm

  RETURN

  END SUBROUTINE compute_kalman
 
 
  SUBROUTINE write_restart_kalman

  USE mod_dims
  USE mod_vars
  USE mod_lib
  USE mod_inout
  USE usr_func
  USE usr_vars
  USE MPI
  USE HDF5
  USE mr_io_protocol, only: DistSegmentedFlowMRIPadded
  USE mr_io_parallel_spacetime 

  ! (basic subroutine)

  IMPLICIT NONE

  !--- write time-dependent data after time integration for next restart ---

  CHARACTER(LEN=3) :: next_restart_char
  CHARACTER(LEN=2) :: phase_char
  REAL :: write_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1)
  INTEGER, DIMENSION(2,3) :: bds
  INTEGER :: i
  CHARACTER(LEN=300) :: write_file

  IF (write_restart_yes.eqv..FALSE.) RETURN

  IF (rank == 0) WRITE(*,'(a,i3,a)') 'writing kalman data for restart',restart,' ...'
  CALL num_to_string(3,restart,next_restart_char)
  !===========================================================================================================
  !=== write means_ and covariances_fields for visualization =================================================
  !===========================================================================================================

  if (allocated(mean_f)) then
    bds(1,1) = lbound(mean_f,2)
    bds(2,1) = ubound(mean_f,2)
    bds(1,2) = lbound(mean_f,3)
    bds(2,2) = ubound(mean_f,3)
    bds(1,3) = lbound(mean_f,4)
    bds(2,3) = ubound(mean_f,4)
  end if
 
  do i = 1,intervals 

    CALL num_to_string(2,i,phase_char)

    write_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1) = 0.0
    if (allocated(mean_f)) then 
      write_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1) = mean_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1)
    end if
    CALL write_hdf('./kf_result/meanf_X_phase'//phase_char//'_restart.'//next_restart_char,'meanfX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_f(b1L,b2L,b3L,1))

    write_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1) = 0.0
    if (allocated(mean_f)) then 
      write_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1) = mean_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),2)
    end if
    CALL write_hdf('./kf_result/meanf_Y_phase'//phase_char//'_restart.'//next_restart_char,'meanfY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_f(b1L,b2L,b3L,1))

    write_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1) = 0.0
    if (allocated(mean_f)) then 
      write_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1) = mean_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),3)
    end if
    CALL write_hdf('./kf_result/meanf_Z_phase'//phase_char//'_restart.'//next_restart_char,'meanfZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_f(b1L,b2L,b3L,1))

    write_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1) = 0.0
    if (allocated(covar_f)) then
      write_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1) = covar_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1)
    end if
    CALL write_hdf('./kf_result/covarf_XX_phase'//phase_char//'_restart.'//next_restart_char,'covarfXX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_f(b1L,b2L,b3L,1))

    write_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1) = 0.0
    if (allocated(covar_f)) then
      write_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1) = covar_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),2)
    end if
    CALL write_hdf('./kf_result/covarf_YY_phase'//phase_char//'_restart.'//next_restart_char,'covarfYY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_f(b1L,b2L,b3L,1))

    write_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1) = 0.0
    if (allocated(covar_f)) then
      write_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1) = covar_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),3)
    end if
    CALL write_hdf('./kf_result/covarf_ZZ_phase'//phase_char//'_restart.'//next_restart_char,'covarfZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_f(b1L,b2L,b3L,1))

    write_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1) = 0.0
    if (allocated(covar_f)) then
      write_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1) = covar_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),4)
    end if
    CALL write_hdf('./kf_result/covarf_XY_phase'//phase_char//'_restart.'//next_restart_char,'covarfXY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_f(b1L,b2L,b3L,1))

    write_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1) = 0.0
    if (allocated(covar_f)) then
      write_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1) = covar_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),5)
    end if
    CALL write_hdf('./kf_result/covarf_XZ_phase'//phase_char//'_restart.'//next_restart_char,'covarfXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_f(b1L,b2L,b3L,1))

    write_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1) = 0.0
    if (allocated(covar_f)) then
      write_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1) = covar_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),6)
    end if
    CALL write_hdf('./kf_result/covarf_YZ_phase'//phase_char//'_restart.'//next_restart_char,'covarfYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_f(b1L,b2L,b3L,1))
   !===========================================================================================================

  end do 

  RETURN
 
  END SUBROUTINE write_restart_kalman
 

  SUBROUTINE read_restart_kalman

  USE mod_dims
  USE mod_vars
  USE mod_lib
  USE mod_inout
  USE usr_func
  USE usr_vars
  USE MPI
  USE HDF5
  USE mr_io_protocol, only: DistSegmentedFlowMRIPadded
  USE mr_io_parallel_spacetime 

  !--- read time-dependent data from previous restart before time integration starts ---

  IMPLICIT NONE


  CHARACTER(LEN=2) :: phase_char
  CHARACTER(LEN=300) :: read_file
  INTEGER :: i
  REAL :: read_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  INTEGER, DIMENSION(2,3) :: bds

  IF (rank == 0) WRITE(*,'(a,i3,a)') 'reading kalman data for restart',restart,' ...'
  write(*,*) 'balala', rank

  CALL num_to_string(3,restart,restart_char)

  if (allocated(mean_f)) then
    bds(1,1) = lbound(mean_f,2)
    bds(2,1) = ubound(mean_f,2)
    bds(1,2) = lbound(mean_f,3)
    bds(2,2) = ubound(mean_f,3)
    bds(1,3) = lbound(mean_f,4)
    bds(2,3) = ubound(mean_f,4)
  end if
 
  do i = 1,intervals

    CALL num_to_string(2,i,phase_char)

    CALL read2_hdf('./kf_result/meanf_X_phase'//phase_char//'_restart.'//restart_char,'meanfX',S1p,S2p,S3p,N1p,N2p,N3p,0,read_f(b1L,b2L,b3L,1))
    if (allocated(mean_f)) then 
      mean_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1) = read_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1)
    end if

    CALL read2_hdf('./kf_result/meanf_Y_phase'//phase_char//'_restart.'//restart_char,'meanfY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_f(b1L,b2L,b3L,1))
    if (allocated(mean_f)) then 
      mean_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),2) = read_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1)
    end if

    CALL read2_hdf('./kf_result/meanf_Z_phase'//phase_char//'_restart.'//restart_char,'meanfZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_f(b1L,b2L,b3L,1))
    if (allocated(mean_f)) then 
      mean_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),3) = read_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1)
    end if
    
    CALL read2_hdf('./kf_result/covarf_XX_phase'//phase_char//'_restart.'//restart_char,'covarfXX',S1p,S2p,S3p,N1p,N2p,N3p,0,read_f(b1L,b2L,b3L,1))
    if (allocated(covar_f)) then 
      covar_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1) = read_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1)
    end if

    CALL read2_hdf('./kf_result/covarf_YY_phase'//phase_char//'_restart.'//restart_char,'covarfYY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_f(b1L,b2L,b3L,1))
    if (allocated(covar_f)) then 
      covar_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),2) = read_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1)
    end if

    CALL read2_hdf('./kf_result/covarf_ZZ_phase'//phase_char//'_restart.'//restart_char,'covarfZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_f(b1L,b2L,b3L,1))
    if (allocated(covar_f)) then 
      covar_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),3) = read_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1)
    end if

    CALL read2_hdf('./kf_result/covarf_XY_phase'//phase_char//'_restart.'//restart_char,'covarfXY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_f(b1L,b2L,b3L,1))
    if (allocated(covar_f)) then 
      covar_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),4) = read_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1)
    end if

    CALL read2_hdf('./kf_result/covarf_XZ_phase'//phase_char//'_restart.'//restart_char,'covarfXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_f(b1L,b2L,b3L,1))
    if (allocated(covar_f)) then 
      covar_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),5) = read_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1)
    end if

    CALL read2_hdf('./kf_result/covarf_YZ_phase'//phase_char//'_restart.'//restart_char,'covarfYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_f(b1L,b2L,b3L,1))
    if (allocated(covar_f)) then 
      covar_f(i,bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),6) = read_f(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3),1)
    end if

  end do

  RETURN
 
  END SUBROUTINE read_restart_kalman

  SUBROUTINE kalman_gain(var, R, d, K, n, m)

  USE MPI

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: m
  real, intent(in)  :: var(n,m)
  real, intent(inout) :: R(m,m),d(m,1)
  real :: d2(m,1),R2(m,m)
  real, intent(out) :: K(n,1)
  integer :: i,j,IPIV(m),INFO,ITER
  !real  :: swork(m*(m+1)),work(m,1) !dsgesv
  real  :: work(m) !dsysv

  d2(1:m,1) = d(1:m,1)
  R2 = R

  ! work_2 = (H*Sigma*transpose(H)+R)^-1*(d-H*mu_f)
  !call dsgesv(m,1,R,m,IPIV,d2,m,d,m,work,swork,ITER,INFO)
  !if (INFO.ne.0) write(*,*) 'error in kalman gain info',INFO
  !if (ITER.lt.0) write(*,*) 'error in kalman gain iter',ITER

  ! d = (H*Sigma*transpose(H)+R)^-1*(d-H*mu_f)
  !call dsysv('U',m,1,R,m,IPIV,d,m,work,m,INFO)
  if (INFO.ne.0) write(*,*) 'error in kalman gain',INFO

  ! d = (H*Sigma*transpose(H)+R)^-1*(d-H*mu_f)
  call dsysv_rook('U',m,1,R,m,IPIV,d,m,work,m,INFO)
  if (INFO.ne.0) write(*,*) 'error in kalman gain',INFO

  !call dposv('U',m,1,R,m,d,m,INFO)
  !if (INFO.ne.0) write(*,*) 'error in kalman gain',INFO

  !write(*,*) maxval( matmul(R2(:,:),d(:,1)) - d2(:,1) )

  ! K = Sigma*transpose(H)*(H*Sigma*transpose(H)+R)^-1*(d-H*mu_f)
  call dgemm('N', 'N', n, 1, m, 1.0, var, n, d, m, 0.0, K, n)

  RETURN

  END SUBROUTINE kalman_gain
