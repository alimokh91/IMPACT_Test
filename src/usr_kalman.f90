!***********************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
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
  USE usr_vars
  USE mod_inout
  USE MPI
  USE HDF5
  USE mr_io_protocol
  USE mr_io_parallel_spacetime 

  IMPLICIT NONE

  TYPE(kalman_t), pointer :: klmn
  INTEGER :: i,ii,n
  !added for writing turb_statxz_ with M1 M2 M3 ================================================================
  CHARACTER(LEN=2) ::  phs
  CHARACTER(LEN=3) ::  M1_char
  CHARACTER(LEN=3) ::  M2_char
  CHARACTER(LEN=3) ::  M3_char

  if (rank.eq.0) write(0,*) "open kalman... "

  call validate_input_kalman() 

  ALLOCATE(write_gain (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)); write_gain = 0.
  ALLOCATE(mean_f(1:intervals,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)); mean_f = 0. 
  ALLOCATE(covar_f(1:intervals,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)); covar_f = 0.
  covar_f(1:intervals,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3) = 1.0e1
  covar_f(1:intervals,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),4:6) = 1.0

  write_kalm_count = 0
  phase = 1

  ALLOCATE(dtime_kalm_phases(1:intervals)); dtime_kalm_phases = 0.
 
  !===========================================================================================================
  !added for writing turb_statxz_ with M1 M2 M3===============================================================
  !===========================================================================================================
  IF (rank == 0 .AND. dtime_out_kalm /= 0.) THEN
     CALL num_to_string(3,restart,restart_char)
     CALL num_to_string(3,M1,M1_char)
     CALL num_to_string(3,M2,M2_char)
     CALL num_to_string(3,M3,M3_char)
     CALL system('mkdir -p kf_result')
     OPEN(63,FILE='./kf_result/tke_window_kalm_'//restart_char//'_'//M1_char//'x'//M2_char//'x'//M3_char//'.txt',STATUS='UNKNOWN')
     do n = 1,intervals
        CALL num_to_string(2,n,phs)
        CALL system('mkdir -p kf_result/phase_'//phs )
     end do
  END IF

  !===========================================================================================================
  !=== read the observed data (components, time, x, y, z) ====================================================
  !===========================================================================================================
  if (associated(mri_inst).eqv..false.) then
     nullify(mri_inst); allocate(mri_inst)
     mri_inst%domain_padding = kalman_domain_padding
     CALL mr_io_read_parallel_flow_padded(MPI_COMM_WORLD, MPI_INFO_NULL, (/NB1,NB2,NB3/), &
          kalman_mri_input_file_path, mri_inst)
     CALL h5open_f(herror) ! Required as hpc-predict-io closes HDF5 environment with h5close_f
  end if

  mri_inst%mri%velocity_mean%array = mri_inst%mri%velocity_mean%array/U_ref
  mri_inst%mri%velocity_cov%array  = mri_inst%mri%velocity_cov%array/U_ref/U_ref

  dtime_kalm_phases(1:intervals-1) = mri_inst%mri%t_coordinates(2:intervals) - mri_inst%mri%t_coordinates(1:intervals-1)
  dtime_kalm_phases(intervals) = kalman_mri_input_attr_t_heart_cycle_period - &
                                 mri_inst%mri%t_coordinates(intervals) + mri_inst%mri%t_coordinates(1)
  !===========================================================================================================

  !===========================================================================================================
  !=== construct matrices x_f, H, R, P_f =====================================================================
  !===========================================================================================================
  allocate(kalman_first)
  klmn => kalman_first
  klmn%m = (kalman_num_spatial_refinements(1)+1)*(kalman_num_spatial_refinements(2)+1)*(kalman_num_spatial_refinements(3)+1)
  NULLIFY(klmn%muf);       ALLOCATE(klmn%muf(1:3*klmn%m)); klmn%muf = 0.
  NULLIFY(klmn%pf);        ALLOCATE(klmn%pf(1:3*klmn%m,1:3*klmn%m)); klmn%pf = 0.
  NULLIFY(klmn%obs_data);  ALLOCATE(klmn%obs_data(1:3))
  NULLIFY(klmn%obs_covar); ALLOCATE(klmn%obs_covar(1:3,1:3))
  NULLIFY(klmn%obs_oper);  ALLOCATE(klmn%obs_oper(1:3,1:3*klmn%m)); klmn%obs_oper = 0.
  NULLIFY(klmn%K);         ALLOCATE(klmn%K(1:3*klmn%m,1:3)); klmn%K = 0.
  do n = 1,klmn%m
     do ii = 1,3
        klmn%obs_oper(ii,3*(n-1)+ii) = 1./klmn%m
     end do
  end do
  NULLIFY(klmn%next)
  !===========================================================================================================

  !===========================================================================================================
  !=== construct pressure-grid mri data-structure ============================================================
  !===========================================================================================================
  nullify(mri_dest); allocate(mri_dest)

  !time informations
  mri_dest%t_dim = mri_inst%mri%t_dim*kalman_num_time_refinements
  allocate(mri_dest%t_coordinates(1:size(mri_inst%mri%t_coordinates)*kalman_num_time_refinements))
  do i = 0, size(mri_inst%mri%t_coordinates)-1
      mri_dest%t_coordinates(i*kalman_num_time_refinements+1:(i+1)*kalman_num_time_refinements) = mri_inst%mri%t_coordinates(i+1)
  end do
  mri_dest%vector_feature%time_offset = mri_inst%mri%velocity_mean%time_offset*kalman_num_time_refinements
  mri_dest%vector_feature%time_dim    = mri_inst%mri%velocity_mean%time_dim*kalman_num_time_refinements

  !flow features
  allocate(mri_dest%vector_feature%array( &
            lbound(mri_inst%mri%velocity_mean%array,1)                                       :ubound(mri_inst%mri%velocity_mean%array,1), &
           (lbound(mri_inst%mri%velocity_mean%array,2)-1)*kalman_num_time_refinements+1      :ubound(mri_inst%mri%velocity_mean%array,2)*kalman_num_time_refinements, &
           (lbound(mri_inst%mri%velocity_mean%array,3)-1)*kalman_num_spatial_refinements(1)+1:ubound(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1), &
           (lbound(mri_inst%mri%velocity_mean%array,4)-1)*kalman_num_spatial_refinements(2)+1:ubound(mri_inst%mri%velocity_mean%array,4)*kalman_num_spatial_refinements(2), &
           (lbound(mri_inst%mri%velocity_mean%array,5)-1)*kalman_num_spatial_refinements(3)+1:ubound(mri_inst%mri%velocity_mean%array,5)*kalman_num_spatial_refinements(3) ))
  mri_dest%vector_feature%array = 0.

  mri_dest%vector_feature%offset = (/ mri_inst%mri%velocity_mean%offset(1)*kalman_num_spatial_refinements(1), &
                                     mri_inst%mri%velocity_mean%offset(2)*kalman_num_spatial_refinements(2), &
                                     mri_inst%mri%velocity_mean%offset(3)*kalman_num_spatial_refinements(3) /)
  mri_dest%vector_feature%dims =   (/ mri_inst%mri%velocity_mean%dims(1)*kalman_num_spatial_refinements(1), &
                                     mri_inst%mri%velocity_mean%dims(2)*kalman_num_spatial_refinements(2), &
                                     mri_inst%mri%velocity_mean%dims(3)*kalman_num_spatial_refinements(3) /)

  !spatial informations
  mri_dest%x_dim = mri_inst%mri%x_dim*kalman_num_spatial_refinements(1)
  mri_dest%y_dim = mri_inst%mri%y_dim*kalman_num_spatial_refinements(2)
  mri_dest%z_dim = mri_inst%mri%z_dim*kalman_num_spatial_refinements(3)

  allocate(mri_dest%x_coordinates(1:size(mri_inst%mri%x_coordinates)*kalman_num_spatial_refinements(1)))
  allocate(mri_dest%y_coordinates(1:size(mri_inst%mri%y_coordinates)*kalman_num_spatial_refinements(2)))
  allocate(mri_dest%z_coordinates(1:size(mri_inst%mri%z_coordinates)*kalman_num_spatial_refinements(3)))

  mri_dest%x_coordinates = 0.5*y1p( mri_inst%domain_padding%lhs(1)*kalman_num_spatial_refinements(1)+1:                  &
                                   (mri_inst%domain_padding%lhs(1)+size(mri_inst%mri%x_coordinates))*kalman_num_spatial_refinements(1)+0) + & 
                           0.5*y1p( mri_inst%domain_padding%lhs(1)*kalman_num_spatial_refinements(1)+2:                  &
                                   (mri_inst%domain_padding%lhs(1)+size(mri_inst%mri%x_coordinates))*kalman_num_spatial_refinements(1)+1)

  mri_dest%y_coordinates = 0.5*y2p( mri_inst%domain_padding%lhs(2)*kalman_num_spatial_refinements(2)+1:                  &
                                   (mri_inst%domain_padding%lhs(2)+size(mri_inst%mri%y_coordinates))*kalman_num_spatial_refinements(2)+0) + & 
                           0.5*y2p( mri_inst%domain_padding%lhs(2)*kalman_num_spatial_refinements(2)+2:                  &
                                   (mri_inst%domain_padding%lhs(2)+size(mri_inst%mri%y_coordinates))*kalman_num_spatial_refinements(2)+1)

  mri_dest%z_coordinates = 0.5*y3p( mri_inst%domain_padding%lhs(3)*kalman_num_spatial_refinements(3)+1:                  &
                                   (mri_inst%domain_padding%lhs(3)+size(mri_inst%mri%z_coordinates))*kalman_num_spatial_refinements(3)+0) + & 
                           0.5*y3p( mri_inst%domain_padding%lhs(3)*kalman_num_spatial_refinements(3)+2:                  &
                                   (mri_inst%domain_padding%lhs(3)+size(mri_inst%mri%z_coordinates))*kalman_num_spatial_refinements(3)+1)
  !===========================================================================================================

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

  DEALLOCATE(mean_f,covar_f,write_gain)
  DEALLOCATE(dtime_kalm_phases)

  if (associated(kalman_first)) then 
     klmn => kalman_first
     CALL destroy_kalman(klmn)
     DEALLOCATE(kalman_first)
  end if

  call mr_io_deallocate_dist_flow_mri_padded(mri_inst)
  call mr_io_deallocate_dist_spacetime_mri(mri_dest)


  IF (rank == 0 .AND. dtime_out_kalm /= 0.) THEN
     CLOSE(63)
  END IF

  contains 
  RECURSIVE SUBROUTINE destroy_kalman(klmn)

     implicit none

     TYPE(kalman_t), pointer, intent(inout) :: klmn

     ! ---------------------------------

     if (associated(klmn%next)) then
        call destroy_kalman(klmn%next)
        DEALLOCATE(klmn%next)
     end if

     DEALLOCATE(klmn%muf,klmn%pf)
     DEALLOCATE(klmn%obs_data,klmn%obs_covar,klmn%obs_oper)
     DEALLOCATE(klmn%K)

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
  INTEGER :: i, j, k, ii, jj, kk, n
  INTEGER, DIMENSION(2,3) :: bounds
  INTEGER :: INFO
  REAL    :: work11(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL    :: work22(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL    :: work33(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL    :: TKE(1:3), TKE_global(1:3)

  phase = mod(write_kalm_count,intervals) + 1

  IF (rank == 0) WRITE(*,'(a,i8,a,i8,a)') 'kalman repetition', write_kalm_count/intervals+1, '   for phase', phase,' ...'

  INFO = 0
  !===========================================================================================================
  !=== interpolation of velocity to pressure grid points =====================================================
  !===========================================================================================================
  ! vel(:,:,:,i) --> worki(:,:,:)
  IF (task == 1) CALL interpolate_vel(.TRUE.)
  !===========================================================================================================

  !===========================================================================================================
  !=== forecast statistics ==================================================================================
  !===========================================================================================================
  !--------------------------------------------------------------------------------------------------------
  !--- <u_i u_j>(t_k) = <u'_i u'_j>(t_k) + <u_i>(t_k) <u_j>(t_k) ------------------------------------------
  !--------------------------------------------------------------------------------------------------------
  covar_f(phase,:,:,:,1) = covar_f(phase,:,:,:,1) + mean_f(phase,:,:,:,1)*mean_f(phase,:,:,:,1)
  covar_f(phase,:,:,:,2) = covar_f(phase,:,:,:,2) + mean_f(phase,:,:,:,2)*mean_f(phase,:,:,:,2)
  covar_f(phase,:,:,:,3) = covar_f(phase,:,:,:,3) + mean_f(phase,:,:,:,3)*mean_f(phase,:,:,:,3)
  covar_f(phase,:,:,:,4) = covar_f(phase,:,:,:,4) + mean_f(phase,:,:,:,1)*mean_f(phase,:,:,:,2)
  covar_f(phase,:,:,:,5) = covar_f(phase,:,:,:,5) + mean_f(phase,:,:,:,1)*mean_f(phase,:,:,:,3)
  covar_f(phase,:,:,:,6) = covar_f(phase,:,:,:,6) + mean_f(phase,:,:,:,2)*mean_f(phase,:,:,:,3)
 
  !--------------------------------------------------------------------------------------------------------
  !--- <u_i>(t_k) -----------------------------------------------------------------------------------------
  !--- <u_i u_j>(t_k) -------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
   mean_f(phase,:,:,:,:) =  mean_f(phase,:,:,:,:)*(write_kalm_count/intervals  )
  covar_f(phase,:,:,:,:) = covar_f(phase,:,:,:,:)*(write_kalm_count/intervals+1)
  DO k = b3L,(N3+b3U)
     DO j = b2L,(N2+b2U)
        DO i = b1L,(N1+b1U)
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
   mean_f(phase,:,:,:,:) =  mean_f(phase,:,:,:,:)/(write_kalm_count/intervals+1)
  covar_f(phase,:,:,:,:) = covar_f(phase,:,:,:,:)/(write_kalm_count/intervals+2)
  !-----------------------------------------------------------------------------------------------------------
  !--- <u'_i u'_j>(t_k) --------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  covar_f(phase,:,:,:,1) = covar_f(phase,:,:,:,1) - mean_f(phase,:,:,:,1)*mean_f(phase,:,:,:,1)
  covar_f(phase,:,:,:,2) = covar_f(phase,:,:,:,2) - mean_f(phase,:,:,:,2)*mean_f(phase,:,:,:,2)
  covar_f(phase,:,:,:,3) = covar_f(phase,:,:,:,3) - mean_f(phase,:,:,:,3)*mean_f(phase,:,:,:,3)
  covar_f(phase,:,:,:,4) = covar_f(phase,:,:,:,4) - mean_f(phase,:,:,:,1)*mean_f(phase,:,:,:,2)
  covar_f(phase,:,:,:,5) = covar_f(phase,:,:,:,5) - mean_f(phase,:,:,:,1)*mean_f(phase,:,:,:,3)
  covar_f(phase,:,:,:,6) = covar_f(phase,:,:,:,6) - mean_f(phase,:,:,:,2)*mean_f(phase,:,:,:,3)
  !===========================================================================================================

  bounds(1,1) = lbound(mri_inst%mri%velocity_mean%array,3)
  bounds(2,1) = ubound(mri_inst%mri%velocity_mean%array,3)
  bounds(1,2) = lbound(mri_inst%mri%velocity_mean%array,4)
  bounds(2,2) = ubound(mri_inst%mri%velocity_mean%array,4)
  bounds(1,3) = lbound(mri_inst%mri%velocity_mean%array,5)
  bounds(2,3) = ubound(mri_inst%mri%velocity_mean%array,5)

  !===========================================================================================================
  !=== integral energy for the Kalman window =================================================================
  !===========================================================================================================
  TKE = 0.
  if (sum(bounds(2,1:3)).ne.0) then
  DO k = (bounds(1,3)-1)*kalman_num_spatial_refinements(3)+1,bounds(2,3)*kalman_num_spatial_refinements(3)+1
     DO j = (bounds(1,2)-1)*kalman_num_spatial_refinements(2)+1,bounds(2,2)*kalman_num_spatial_refinements(2)+1
        DO i = (bounds(1,1)-1)*kalman_num_spatial_refinements(1)+1,bounds(2,1)*kalman_num_spatial_refinements(1)+1
           if ( k.le.N3p .and. j.le.N2p .and. i.le.N1p) then
              TKE(1) = TKE(1) + dx1p(i)*dx2p(j)*dx3p(k)*(mean_f(phase,i,j,k,1)**2+mean_f(phase,i,j,k,2)**2+&
                                                         mean_f(phase,i,j,k,3)**2)
        !--------------------------------------------------------------------------------------------------------
        !--- u'_i u'_i(t_k) -------------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------------------------------
              TKE(2) = TKE(2) + dx1p(i)*dx2p(j)*dx3p(k)*( (work1(i,j,k) - mean_f(phase,i,j,k,1))**2 + &
                                                          (work2(i,j,k) - mean_f(phase,i,j,k,2))**2 + &
                                                          (work3(i,j,k) - mean_f(phase,i,j,k,3))**2 )
        !--------------------------------------------------------------------------------------------------------
        !--- <u'_i u'_i>(t_k) -----------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------------------------------
              TKE(3) = TKE(3) + dx1p(i)*dx2p(j)*dx3p(k)*( covar_f(phase,i,j,k,1) + covar_f(phase,i,j,k,2) + &
                                                          covar_f(phase,i,j,k,3) )
           end if
        END DO
     END DO
  END DO
  end if
  CALL MPI_ALLREDUCE(TKE,TKE_global,3,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  TKE_global = TKE_global/2.

  IF (rank == 0) THEN
    WRITE(63,'(4E25.17)') time, TKE_global(1:3)
    CALL flush(63)
  END IF
  !===========================================================================================================


  !===========================================================================================================
  !=== apply kalman filtering ================================================================================
  !===========================================================================================================
  work11(:,:,:) = work1(:,:,:)
  work22(:,:,:) = work2(:,:,:)
  work33(:,:,:) = work3(:,:,:)
  write_gain = 0.
  klmn => kalman_first
  do while(associated(klmn))
    DO k = bounds(1,3), bounds(2,3)
      DO j = bounds(1,2), bounds(2,2)
        DO i = bounds(1,1), bounds(2,1)
          n = 0
          DO kk = (k-1)*kalman_num_spatial_refinements(3)+1,k*kalman_num_spatial_refinements(3)+1
            DO jj = (j-1)*kalman_num_spatial_refinements(2)+1,j*kalman_num_spatial_refinements(2)+1
              DO ii = (i-1)*kalman_num_spatial_refinements(1)+1,i*kalman_num_spatial_refinements(1)+1
                n = n + 1
    !===========================================================================================================
    !=== Fill the matrices and the vectors involved in the kalman filter algorithm =============================
    !===========================================================================================================
    !--------------------------------------------------------------------------------------------------------
    !--- u(x,t) ---------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------
                klmn%muf(3*(n-1)+1) = work1(ii,jj,kk)
                klmn%muf(3*(n-1)+2) = work2(ii,jj,kk)
                klmn%muf(3*(n-1)+3) = work3(ii,jj,kk)
    !--------------------------------------------------------------------------------------------------------
    !--- fill the state vector muf and the covariance matrix pf ---------------------------------------------
    !--- <u'_i u'_j + init_covar>(t_k) ----------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------
                klmn%pf(3*(n-1)+1,3*(n-1)+1) = covar_f(phase,ii,jj,kk,1)
                klmn%pf(3*(n-1)+2,3*(n-1)+2) = covar_f(phase,ii,jj,kk,2)
                klmn%pf(3*(n-1)+3,3*(n-1)+3) = covar_f(phase,ii,jj,kk,3)
                klmn%pf(3*(n-1)+1,3*(n-1)+2) = covar_f(phase,ii,jj,kk,4)
                klmn%pf(3*(n-1)+1,3*(n-1)+3) = covar_f(phase,ii,jj,kk,5)
                klmn%pf(3*(n-1)+2,3*(n-1)+3) = covar_f(phase,ii,jj,kk,6)
                klmn%pf(3*(n-1)+2,3*(n-1)+1) = klmn%pf(3*(n-1)+1,3*(n-1)+2)
                klmn%pf(3*(n-1)+3,3*(n-1)+1) = klmn%pf(3*(n-1)+1,3*(n-1)+3)
                klmn%pf(3*(n-1)+3,3*(n-1)+2) = klmn%pf(3*(n-1)+2,3*(n-1)+3)
              END DO
            END DO
          END DO
    !--------------------------------------------------------------------------------------------------------
    !--- d_k and R_k ----------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------
          klmn%obs_data(1:3) = mri_inst%mri%velocity_mean%array(1:3,phase,i,j,k)
          klmn%obs_covar(1:3,1:3) = mri_inst%mri%velocity_cov%array(1:3,1:3,phase,i,j,k)
    !--------------------------------------------------------------------------------------------------------
    !--- apply kalman filtering- ----------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------
          CALL kalman_filter(klmn%muf,klmn%pf,klmn%obs_oper,klmn%obs_data,klmn%obs_covar, &
                             klmn%K,size(klmn%muf),size(klmn%obs_data),INFO)
  
          n = 0
          DO kk = (k-1)*kalman_num_spatial_refinements(3)+1,k*kalman_num_spatial_refinements(3)+1
            DO jj = (j-1)*kalman_num_spatial_refinements(2)+1,j*kalman_num_spatial_refinements(2)+1
              DO ii = (i-1)*kalman_num_spatial_refinements(1)+1,i*kalman_num_spatial_refinements(1)+1
                n = n + 1
    !--------------------------------------------------------------------------------------------------------
    !--- store u_i(:,:,:) of analysis step ------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------
                  work11(ii,jj,kk) = klmn%muf(3*(n-1)+1)
                  work22(ii,jj,kk) = klmn%muf(3*(n-1)+2)
                  work33(ii,jj,kk) = klmn%muf(3*(n-1)+3)
    !--------------------------------------------------------------------------------------------------------
    !--- write Kalman gain field ----------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------
                  write_gain(ii,jj,kk,1) = klmn%K(3*(n-1)+1,1)
                  write_gain(ii,jj,kk,2) = klmn%K(3*(n-1)+2,2)
                  write_gain(ii,jj,kk,3) = klmn%K(3*(n-1)+3,3)
                  write_gain(ii,jj,kk,4) = klmn%K(3*(n-1)+1,2)
                  write_gain(ii,jj,kk,5) = klmn%K(3*(n-1)+1,3)
                  write_gain(ii,jj,kk,6) = klmn%K(3*(n-1)+2,3)
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
    klmn => klmn%next
  end do
  !===========================================================================================================

  !===========================================================================================================
  !=== then update u_i(:,:,:) ================================================================================
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  !--- worki(:,:,:) --> vel(:,:,:,i) -------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  work1(:,:,:) = work11(:,:,:)
  work2(:,:,:) = work22(:,:,:)
  work3(:,:,:) = work33(:,:,:)
  CALL interpolate2_pre_vel(.TRUE.,1,work1,vel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1))
  CALL interpolate2_pre_vel(.TRUE.,2,work2,vel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),2))
  CALL interpolate2_pre_vel(.TRUE.,3,work3,vel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),3))
  ! exchange across boundaries on vel grid
  CALL exchange_all_all(.TRUE.,vel)
  !===========================================================================================================

  !===========================================================================================================
  !=== write pressure-grid mri data-fields ===================================================================
  !===========================================================================================================
  bounds(1,1) = lbound(mri_dest%vector_feature%array,3)
  bounds(2,1) = ubound(mri_dest%vector_feature%array,3)
  bounds(1,2) = lbound(mri_dest%vector_feature%array,4)
  bounds(2,2) = ubound(mri_dest%vector_feature%array,4)
  bounds(1,3) = lbound(mri_dest%vector_feature%array,5)
  bounds(2,3) = ubound(mri_dest%vector_feature%array,5)

  !--------------------------------------------------------------------------------------------------------
  !--- <u_i>(t_k) -----------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
  mri_dest%vector_feature%array (:,phase,:,:,:) = mri_dest%vector_feature%array (:,phase,:,:,:)*(write_kalm_count/intervals)
  DO k = bounds(1,3), bounds(2,3)
    DO j = bounds(1,2), bounds(2,2)
      DO i = bounds(1,1), bounds(2,1)
        DO kk = k,k+1
          DO jj = j,j+1
            DO ii = i,i+1
              mri_dest%vector_feature%array(1,phase,i,j,k) = mri_dest%vector_feature%array(1,phase,i,j,k) + 0.125*work1(ii,jj,kk)
              mri_dest%vector_feature%array(2,phase,i,j,k) = mri_dest%vector_feature%array(2,phase,i,j,k) + 0.125*work2(ii,jj,kk)
              mri_dest%vector_feature%array(3,phase,i,j,k) = mri_dest%vector_feature%array(3,phase,i,j,k) + 0.125*work3(ii,jj,kk)
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
  mri_dest%vector_feature%array (:,phase,:,:,:) = mri_dest%vector_feature%array (:,phase,:,:,:)/(write_kalm_count/intervals+1)
  !===========================================================================================================

  write_kalm_count = write_kalm_count + 1
  dtime_out_kalm = dtime_kalm_phases(phase)
  time_out_kalm = time_out_kalm + dtime_out_kalm
  !dtime_out_vect = dtime_out_kalm/4.
  write_out_kalm = .FALSE.

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
  USE mr_io_protocol, only: DistFlowMRIPadded
  USE mr_io_parallel_spacetime 

  ! (basic subroutine)

  IMPLICIT NONE

  !--- write time-dependent data after time integration for next restart ---

  CHARACTER(LEN=3) :: next_restart_char
  CHARACTER(LEN=2) :: phase_char
  REAL :: write_mean_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  REAL :: write_covar_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)

  IF (write_restart_yes) THEN

     IF (rank == 0) WRITE(*,'(a,i3,a)') 'writing kalman data for restart',restart,' ...'
     CALL num_to_string(3,restart,next_restart_char)
     !===========================================================================================================
     !=== write means_ and covariances_fields for visualization =================================================
     !===========================================================================================================
     do phase = 1,intervals 

        CALL num_to_string(2,phase,phase_char)

        write_mean_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3) = mean_f(phase,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
        CALL write_hdf('meanf_X_phase'//phase_char//'_restart.'//next_restart_char,'meanfX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean_f(b1L,b2L,b3L,1))
        CALL write_hdf('meanf_Y_phase'//phase_char//'_restart.'//next_restart_char,'meanfY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean_f(b1L,b2L,b3L,2))
        CALL write_hdf('meanf_Z_phase'//phase_char//'_restart.'//next_restart_char,'meanfZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean_f(b1L,b2L,b3L,3))

        write_covar_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6) = covar_f(phase,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)
        CALL write_hdf('covarf_XX_phase'//phase_char//'_restart.'//next_restart_char,'covarfXX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_f(b1L,b2L,b3L,1))
        CALL write_hdf('covarf_YY_phase'//phase_char//'_restart.'//next_restart_char,'covarfYY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_f(b1L,b2L,b3L,2))
        CALL write_hdf('covarf_ZZ_phase'//phase_char//'_restart.'//next_restart_char,'covarfZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_f(b1L,b2L,b3L,3))
        CALL write_hdf('covarf_XY_phase'//phase_char//'_restart.'//next_restart_char,'covarfXY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_f(b1L,b2L,b3L,4))
        CALL write_hdf('covarf_XZ_phase'//phase_char//'_restart.'//next_restart_char,'covarfXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_f(b1L,b2L,b3L,5))
        CALL write_hdf('covarf_YZ_phase'//phase_char//'_restart.'//next_restart_char,'covarfYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_f(b1L,b2L,b3L,6))
       !===========================================================================================================

     end do 

     CALL mr_io_write_parallel_spacetime(MPI_COMM_WORLD, MPI_INFO_NULL, kalman_mri_output_file_path, mri_dest)
     CALL h5open_f(herror)

  END IF

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
  USE mr_io_protocol, only: DistFlowMRIPadded
  USE mr_io_parallel_spacetime 

  !--- read time-dependent data from previous restart before time integration starts ---

  IMPLICIT NONE


  CHARACTER(LEN=2) :: phase_char
  INTEGER :: i
  REAL :: read_mean_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  REAL :: read_covar_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)


  IF (rank == 0) WRITE(*,'(a,i3,a)') 'reading kalman data for restart',restart,' ...'

  CALL num_to_string(3,restart,restart_char)

  do i = 1,intervals

     CALL num_to_string(2,phase,phase_char)

     CALL read2_hdf('meanf_X_phase'//phase_char//'_restart.'//restart_char,'meanfX',S1p,S2p,S3p,N1p,N2p,N3p,0,read_mean_f(b1L,b2L,b3L,1))
     CALL read2_hdf('meanf_Y_phase'//phase_char//'_restart.'//restart_char,'meanfY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_mean_f(b1L,b2L,b3L,2))
     CALL read2_hdf('meanf_Z_phase'//phase_char//'_restart.'//restart_char,'meanfZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_mean_f(b1L,b2L,b3L,3))
     
     mean_f(i,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3) = read_mean_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)

     CALL read2_hdf('covarf_XX_phase'//phase_char//'_restart.'//restart_char,'covarfXX',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_f(b1L,b2L,b3L,1))
     CALL read2_hdf('covarf_YY_phase'//phase_char//'_restart.'//restart_char,'covarfYY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_f(b1L,b2L,b3L,2))
     CALL read2_hdf('covarf_ZZ_phase'//phase_char//'_restart.'//restart_char,'covarfZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_f(b1L,b2L,b3L,3))
     CALL read2_hdf('covarf_XY_phase'//phase_char//'_restart.'//restart_char,'covarfXY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_f(b1L,b2L,b3L,4))
     CALL read2_hdf('covarf_XZ_phase'//phase_char//'_restart.'//restart_char,'covarfXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_f(b1L,b2L,b3L,5))
     CALL read2_hdf('covarf_YZ_phase'//phase_char//'_restart.'//restart_char,'covarfYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_f(b1L,b2L,b3L,6))

     covar_f(i,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6) = read_covar_f(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)

  end do

  CALL mr_io_read_parallel_spacetime(MPI_COMM_WORLD, MPI_INFO_NULL, (/NB1,NB2,NB3/), kalman_mri_output_file_path, mri_dest)
  CALL h5open_f(herror)

  RETURN
 
  END SUBROUTINE read_restart_kalman


  SUBROUTINE kalman_filter(mu, Sigma, H, d, R, K, m, n, INFO)

  implicit none

  integer, intent(in) :: m
  integer, intent(in) :: n
  real, intent(in)  :: H(n, m)
  real, intent(in)  :: d(n)
  real, intent(in)  :: R(n, n)
  real, intent(inout) :: mu(m)
  real, intent(inout) :: Sigma(m, m)
  real, intent(out) :: K(m, n)
  integer :: IPIV(n)
  integer, intent(out) :: INFO
  real    :: work(n)
  real    :: H_trasp(m, n)
  real    :: mu_2(m,1)
  real    :: dvar_2(n,1)
  real    :: Rvar_2(n, n)
  real    :: var_17(m, n)
!  real    :: var_19(m, m)

  H_trasp = transpose(H)
  dvar_2(1:n,1) = d(1:n)
  mu_2(1:m,1) = mu(1:m)
  Rvar_2 = R

  call dgemm('N', 'N', n, 1, m, -1.0, H, n, mu_2, m, 1.0, dvar_2, n)    ! dvar_2 = - H*mu_f + data
  call dsymm('L', 'U', m, n, 1.0, Sigma, m, H_trasp, m, 0.0, var_17, m) ! var_17 = Sigma*H'
  call dgemm('N', 'N', n, n, m, 1.0, H, n, var_17, m, 1.0, Rvar_2, n)   ! Rvar_2 = H*Sigma*H'+R

  call dgetrf(n, n, Rvar_2, n, IPIV, INFO)
  if (INFO.ne.0) then
     write(*,*) 'matrix is numerically singular'
  end if
  call dgetri(n, Rvar_2, n, IPIV, work, n, INFO)                        ! Rvar_2 = (H*Sigma*H'+R)^-1
  if (INFO.ne.0) then
    write(*,*) 'matrix inversion failed'
  end if

  call dgemm('N', 'N', m, n, n, 1.0, var_17, m, Rvar_2, n, 0.0, K, m)  ! K = Sigma*H'*(H*Sigma*H'+R)^-1
  call dgemm('N', 'N', m, 1, n, 1.0, K, m, dvar_2, n, 1.0, mu_2, m)    ! mu_a = mu_f + K*(data - H*mu_f)

  mu(1:m) = mu_2(1:m,1)

!  kalman filter: calculation of p_a
!  call dgemm('N', 'N', m, m, n, 1.0, K, m, H, n, 0.0, var_19, m)       ! var_19  = Sigma*H'*(H*Sigma*H'+R)^-1*H
!  call dgemm('N', 'N', m, m, n, -1.0, var, m, Sigma, n, 0.0, Sigma, m) ! Sigma_a = - Sigma*H'*(H*Sigma*H'+R)^-1*H*Sigma + Sigma
 
  RETURN

  END SUBROUTINE kalman_filter
