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
  USE usr_func
  USE usr_vars
  USE mod_lib
  USE mod_inout
  USE MPI
  USE HDF5
  USE mr_io_protocol, only: DistFlowMRIPadded
  USE mr_io_parallel_spacetime 

  IMPLICIT NONE

  TYPE(kalman_t), pointer :: klmn
  INTEGER :: i,j,k,ii,jj,kk
  INTEGER, DIMENSION(2,5) :: bounds_kalman
  INTEGER :: m,n,l
  INTEGER :: n_procs
  INTEGER, ALLOCATABLE :: n_data_count(:)
  INTEGER, ALLOCATABLE :: i_data(:,:,:)
  INTEGER, ALLOCATABLE :: flag_data(:,:,:)
  REAL :: d_vol
  !added for writing turb_statxz_ with M1 M2 M3 ================================================================
  CHARACTER(LEN=2) ::  phs
  CHARACTER(LEN=3) ::  M1_char
  CHARACTER(LEN=3) ::  M2_char
  CHARACTER(LEN=3) ::  M3_char
  CHARACTER(LEN=50) :: read_dir,write_dir

  call validate_input_kalman() 

  ALLOCATE(mean_gbl(1:intervals,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)); mean_gbl = 0. 
  ALLOCATE(covar_gbl(1:intervals,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)); covar_gbl = 0. 
  ALLOCATE(write_gain (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)); write_gain = 0.

  ALLOCATE(dtime_kalm_phases(1:intervals)); dtime_kalm_phases = 0.

  ALLOCATE(i_data(S1p:N1p,S2p:N2p,S3p:N3p)); i_data = 0
  ALLOCATE(flag_data(S1p:N1p,S2p:N2p,S3p:N3p)); flag_data = 0
  ALLOCATE(klm_flag(S1p:N1p,S2p:N2p,S3p:N3p)); klm_flag = 0

  if (associated(mri_inst).eqv..false.) then
     nullify(mri_inst); allocate(mri_inst)
     mri_inst%domain_padding = kalman_domain_padding
     CALL mr_io_read_parallel_flow_padded(MPI_COMM_WORLD, MPI_INFO_NULL, (/NB1,NB2,NB3/), &
          kalman_mri_input_file_path, mri_inst)
     CALL h5open_f(herror) ! Required as hpc-predict-io closes HDF5 environment with h5close_f
  end if
  bounds_kalman(1,1:5) = lbound(mri_inst%mri%velocity_mean%array)
  bounds_kalman(2,1:5) = ubound(mri_inst%mri%velocity_mean%array)

  dtime_kalm_phases(1:intervals-1) = mri_inst%mri%t_coordinates(2:intervals) - mri_inst%mri%t_coordinates(1:intervals-1)
  dtime_kalm_phases(intervals) = kalman_mri_input_attr_t_heart_cycle_period - &
                                 mri_inst%mri%t_coordinates(intervals) + mri_inst%mri%t_coordinates(1)

  !===========================================================================================================
  !=== read the observed data (index, x, y, z) ===============================================================
  !===========================================================================================================
  l = 0
  DO k = bounds_kalman(1,5), bounds_kalman(2,5)
     DO j = bounds_kalman(1,4), bounds_kalman(2,4)
        DO i = bounds_kalman(1,3), bounds_kalman(2,3)
           klm_flag = 1
           l = l + 1
           i_data((i-1)*kalman_num_spatial_refinements(1)+1:i*kalman_num_spatial_refinements(1)+1, &
                  (j-1)*kalman_num_spatial_refinements(2)+1:j*kalman_num_spatial_refinements(2)+1, &
                  (k-1)*kalman_num_spatial_refinements(3)+1:k*kalman_num_spatial_refinements(3)+1) = l
        END DO 
     END DO
  END DO

  !===========================================================================================================
  !=== construct matrices x_f, H, R, P_f =====================================================================
  !===========================================================================================================
  NULLIFY(kalman_first)
  n = 0
  DO k = bounds_kalman(1,5), bounds_kalman(2,5)
     DO j = bounds_kalman(1,4), bounds_kalman(2,4)
        DO i = bounds_kalman(1,3), bounds_kalman(2,3)
           n = n + 1
           write(*,*) n
           l = count(i_data(S1p:N1p,S2p:N2p,S3p:N3p).eq.n)
           if (.NOT.associated(kalman_first)) then
              allocate(kalman_first)
              klmn => kalman_first
           else
              allocate(klmn%next)
              klmn => klmn%next
           end if
           klmn%i_data = i
           klmn%j_data = j
           klmn%k_data = k
           klmn%m = l
           NULLIFY(klmn%mean);      ALLOCATE(klmn%mean (1:intervals,1:klmn%m,1:3)); klmn%mean = 0.
           NULLIFY(klmn%covar);     ALLOCATE(klmn%covar(1:intervals,1:klmn%m,1:6)); klmn%covar = 0.
           NULLIFY(klmn%muf);       ALLOCATE(klmn%muf(1:3*klmn%m)); klmn%muf = 0.
           NULLIFY(klmn%pf);        ALLOCATE(klmn%pf(1:3*klmn%m,1:3*klmn%m)); klmn%pf = 0.
           NULLIFY(klmn%obs_data);  ALLOCATE(klmn%obs_data(1:3))
           NULLIFY(klmn%obs_covar); ALLOCATE(klmn%obs_covar(1:3,1:3))
           NULLIFY(klmn%obs_oper);  ALLOCATE(klmn%obs_oper(1:3,1:3*klmn%m)); klmn%obs_oper = 0.
           NULLIFY(klmn%K);         ALLOCATE(klmn%K(1:3*klmn%m,1:3)); klmn%K = 0.
           NULLIFY(klmn%x);         ALLOCATE(klmn%x(1:klmn%m))
           NULLIFY(klmn%y);         ALLOCATE(klmn%y(1:klmn%m))
           NULLIFY(klmn%z);         ALLOCATE(klmn%z(1:klmn%m))
           l = 0
           DO kk = S3p, N3p
              DO jj = S2p, N2p
                 DO ii = S1p, N1p
                    if (n.eq.i_data(ii,jj,kk)) then
                       l = l + 1
                       klmn%x(l) = ii
                       klmn%y(l) = jj
                       klmn%z(l) = kk
                    end if
                 END DO
              END DO
           END DO
           d_vol = 0.
           do l = 1,klmn%m
              do ii = 1,3
                 klmn%obs_oper(ii,3*(l-1)+ii) = dx1p(klmn%x(l))*dx2p(klmn%y(l))*dx3p(klmn%z(l))
              end do
              d_vol = d_vol + dx1p(klmn%x(l))*dx2p(klmn%y(l))*dx3p(klmn%z(l)) ! (dx3p = 1. for 2D)
           end do
           klmn%covar(1:intervals,1:klmn%m,1:3) = 1.0e1
           klmn%covar(1:intervals,1:klmn%m,4:6) = 1.0
           klmn%obs_oper = klmn%obs_oper/d_vol
           NULLIFY(klmn%next)
        END DO
     END DO
  END DO
  !===========================================================================================================

  DEALLOCATE(i_data,flag_data)

  write_kalm_count = 0

  phase = 1
  ALLOCATE(repetition(1:intervals)); repetition = 0
 
  !added for writing turb_statxz_ with M1 M2 M3===============================================================
  CALL num_to_string(3,M1,M1_char)
  CALL num_to_string(3,M2,M2_char)
  CALL num_to_string(3,M3,M3_char)
  CALL num_to_string(3,restart,restart_char)

  IF (rank == 0 .AND. dtime_out_kalm /= 0.) THEN
     CALL system('mkdir -p kf_result')
     do n = 1,intervals
        CALL num_to_string(2,n,phs)
        CALL system('mkdir -p kf_result/phase_'//phs )
     end do
     OPEN(43,FILE='tke_domain_kalm_'//restart_char//'.txt',STATUS='UNKNOWN')
     OPEN(44,FILE='tke_window_kalm_'//restart_char//'.txt',STATUS='UNKNOWN')
     OPEN(53,FILE='turb_statxz_kalm'//restart_char//'_'//M1_char//'x'//M2_char//'x'//M3_char//'.txt',STATUS='UNKNOWN') !added M1 M2 M3
  END IF

  RETURN

  END SUBROUTINE open_kalman
 


  SUBROUTINE close_kalman
  ! (basic subroutine)

  USE mod_dims
  USE mod_vars
  USE usr_func
  USE usr_vars

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'

  TYPE(kalman_t), pointer :: klmn

  DEALLOCATE(repetition,klm_flag)
  DEALLOCATE(mean_gbl,covar_gbl,write_gain)
  DEALLOCATE(dtime_kalm_phases)

  if (associated(kalman_first)) then 
     klmn => kalman_first
     CALL destroy_kalman(klmn)
     DEALLOCATE(kalman_first)
  end if

  IF (rank == 0 .AND. dtime_out_kalm /= 0.) THEN
     CLOSE(43)
     CLOSE(44)
     CLOSE(53)
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

     DEALLOCATE(klmn%x,klmn%y,klmn%z)
     DEALLOCATE(klmn%mean,klmn%covar)
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

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
  TYPE(kalman_t), pointer :: klmn
  INTEGER :: i, j, k
  INTEGER :: INFO
  REAL    :: TKE(1:2), TKE_global(1:2), FKE(1:2), FKE_global(1:2)
  REAL    :: mean_xz_km_write(1:M2), mean_xz_km(S2p:N2p), mean_xz_km_global(S2p:N2p)
  CHARACTER(LEN=2) :: phs
  CHARACTER(LEN=50) :: read_dir,write_dir

  phase = mod(write_kalm_count,intervals) + 1

  if (time.ge.time_out_kalm) then
     write_kalm_count = write_kalm_count + 1
     dtime_out_kalm = dtime_kalm_phases(phase)
     time_out_kalm = time_out_kalm + dtime_out_kalm
     dtime_out_vect = dtime_out_kalm/3.
  end if

  repetition(phase) = repetition(phase) + 1

  IF (rank == 0) WRITE(*,'(a,i8,a,i8)') 'kalman repetition', repetition(phase), '   for phase', phase,' ...'

  INFO = 0
  !===========================================================================================================
  !=== interpolation of velocity to pressure grid points =====================================================
  !===========================================================================================================
  ! vel(:,:,:,i) --> worki(:,:,:)
  IF (task == 1) CALL interpolate_vel(.TRUE.)
  !===========================================================================================================

  !===========================================================================================================
  !m_stats = m_stats/U_ref
  !c_stats = c_stats/(U_ref*U_ref)

  !===========================================================================================================
  !=== turbulence statistics =================================================================================
  !===========================================================================================================
  klmn => kalman_first
  do while(associated(klmn))
  !--------------------------------------------------------------------------------------------------------
  !--- <u_i>(t_k) -------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
     klmn%mean(phase,:,:) = klmn%mean(phase,:,:)*(repetition(phase)-1)
     do i = 1, klmn%m
        klmn%mean(phase,i,1) = klmn%mean(phase,i,1) + work1(klmn%x(i),klmn%y(i),klmn%z(i))
        klmn%mean(phase,i,2) = klmn%mean(phase,i,2) + work2(klmn%x(i),klmn%y(i),klmn%z(i))
        klmn%mean(phase,i,3) = klmn%mean(phase,i,3) + work3(klmn%x(i),klmn%y(i),klmn%z(i))
     end do
     klmn%mean(phase,:,:) = klmn%mean(phase,:,:)/repetition(phase)
  !--------------------------------------------------------------------------------------------------------
  !--- <u_i u_j>(t_k) ---------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
     klmn%covar(phase,:,:) = klmn%covar(phase,:,:)*repetition(phase)
     do i = 1, klmn%m
        klmn%covar(phase,i,1) = klmn%covar(phase,i,1) + &
                          (work1(klmn%x(i),klmn%y(i),klmn%z(i))*work1(klmn%x(i),klmn%y(i),klmn%z(i))) 
        klmn%covar(phase,i,2) = klmn%covar(phase,i,2) + &
                          (work2(klmn%x(i),klmn%y(i),klmn%z(i))*work2(klmn%x(i),klmn%y(i),klmn%z(i)))
        klmn%covar(phase,i,3) = klmn%covar(phase,i,3) + &
                          (work3(klmn%x(i),klmn%y(i),klmn%z(i))*work3(klmn%x(i),klmn%y(i),klmn%z(i)))
        klmn%covar(phase,i,4) = klmn%covar(phase,i,4) + &
                          (work1(klmn%x(i),klmn%y(i),klmn%z(i))*work2(klmn%x(i),klmn%y(i),klmn%z(i)))
        klmn%covar(phase,i,5) = klmn%covar(phase,i,5) + &
                          (work1(klmn%x(i),klmn%y(i),klmn%z(i))*work3(klmn%x(i),klmn%y(i),klmn%z(i)))
        klmn%covar(phase,i,6) = klmn%covar(phase,i,6) + &
                          (work2(klmn%x(i),klmn%y(i),klmn%z(i))*work3(klmn%x(i),klmn%y(i),klmn%z(i)))
     end do
     klmn%covar(phase,:,:) = klmn%covar(phase,:,:)/(repetition(phase)+1)
  !--------------------------------------------------------------------------------------------------------
     klmn => klmn%next
  end do
  !===========================================================================================================


  !===========================================================================================================
  !=== Fill the matrices and the vectors involved in the kalman filter algorithm =============================
  !=== apply kalman filtering ================================================================================
  !===========================================================================================================
  klmn => kalman_first
  do while(associated(klmn))
     do i = 1, klmn%m
  !--------------------------------------------------------------------------------------------------------
  !--- u(x,t) --------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
        klmn%muf(3*(i-1)+1) = work1(klmn%x(i),klmn%y(i),klmn%z(i))
        klmn%muf(3*(i-1)+2) = work2(klmn%x(i),klmn%y(i),klmn%z(i))
        klmn%muf(3*(i-1)+3) = work3(klmn%x(i),klmn%y(i),klmn%z(i))
  !--------------------------------------------------------------------------------------------------------
  !--- fill the state vector muf and the covariance matrix pf ---------------------------------------------
  !--- <u'_i u'_j>(t_k) = <u_i u_j>(t_k) - <u_i>(t_k) <u_j>(t_k) ------------------------------------------
  !--------------------------------------------------------------------------------------------------------
        klmn%pf(3*(i-1)+1,3*(i-1)+1) = klmn%covar(phase,i,1) - klmn%mean(phase,i,1)*klmn%mean(phase,i,1)
        klmn%pf(3*(i-1)+2,3*(i-1)+2) = klmn%covar(phase,i,2) - klmn%mean(phase,i,2)*klmn%mean(phase,i,2)
        klmn%pf(3*(i-1)+3,3*(i-1)+3) = klmn%covar(phase,i,3) - klmn%mean(phase,i,3)*klmn%mean(phase,i,3)
        klmn%pf(3*(i-1)+1,3*(i-1)+2) = klmn%covar(phase,i,4) - klmn%mean(phase,i,1)*klmn%mean(phase,i,2)
        klmn%pf(3*(i-1)+1,3*(i-1)+3) = klmn%covar(phase,i,5) - klmn%mean(phase,i,1)*klmn%mean(phase,i,3)
        klmn%pf(3*(i-1)+2,3*(i-1)+3) = klmn%covar(phase,i,6) - klmn%mean(phase,i,2)*klmn%mean(phase,i,3)
        klmn%pf(3*(i-1)+2,3*(i-1)+1) = klmn%pf(3*(i-1)+1,3*(i-1)+2)
        klmn%pf(3*(i-1)+3,3*(i-1)+1) = klmn%pf(3*(i-1)+1,3*(i-1)+3)
        klmn%pf(3*(i-1)+3,3*(i-1)+2) = klmn%pf(3*(i-1)+2,3*(i-1)+3)
     end do
  !--------------------------------------------------------------------------------------------------------
  !--- d_k and R_k ----------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
     klmn%obs_data(1:3) = mri_inst%mri%velocity_mean%array(1:3,phase,klmn%i_data,klmn%j_data,klmn%k_data)
     klmn%obs_covar(1:3,1:3) = mri_inst%mri%velocity_cov%array(1:3,1:3,phase,klmn%i_data,klmn%j_data,klmn%k_data)
  !--------------------------------------------------------------------------------------------------------
  !----apply kalman filtering- ----------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
     CALL kalman_filter(klmn%muf,klmn%pf,klmn%obs_oper,klmn%obs_data,klmn%obs_covar, &
                        klmn%K,size(klmn%muf),size(klmn%obs_data),INFO)
  !--------------------------------------------------------------------------------------------------------
  !--- then update u_i(:,:,:) -----------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
     do i = 1, klmn%m
        work1(klmn%x(i),klmn%y(i),klmn%z(i)) = klmn%muf(3*(i-1)+1)
        work2(klmn%x(i),klmn%y(i),klmn%z(i)) = klmn%muf(3*(i-1)+2)
        work3(klmn%x(i),klmn%y(i),klmn%z(i)) = klmn%muf(3*(i-1)+3)
     end do
     klmn => klmn%next
  end do
  !===========================================================================================================

  !===========================================================================================================
  ! === worki(:,:,:) --> vel(:,:,:,i) ========================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.TRUE.,1,work1,vel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1))
  CALL interpolate2_pre_vel(.TRUE.,2,work2,vel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),2))
  CALL interpolate2_pre_vel(.TRUE.,3,work3,vel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),3))
  ! exchange across boundaries on vel grid
  CALL exchange_all_all(.TRUE.,vel)
  !===========================================================================================================

  !===========================================================================================================
  !=== reynolds_tau ==========================================================================================
  !===========================================================================================================
  IF (wallnorm_dir == 1) THEN
     mean_xz_km = 0.
     DO j = S2p, N2p
        DO k = S3p, N3p
           DO i = S1p, N1p
              mean_xz_km(j) = mean_xz_km(j) + work1(i,j,k)*dx1p(i)*dx3p(k) ! (dx3p = 1. for 2D)
           END DO
        END DO
     END DO
     CALL MPI_ALLREDUCE(mean_xz_km,mean_xz_km_global,(N2p-S2p+1),MPI_REAL8,MPI_SUM,COMM_SLICE2,merror)
     mean_xz_km = mean_xz_km_global / (L1*L3)
     CALL MPI_ALLGATHERv(mean_xz_km,(N2p-S2p+1),MPI_REAL8,mean_xz_km_write,bar2_size,bar2_offset,MPI_REAL8,COMM_BAR2,merror)
     IF (rank == 0) THEN
        WRITE(53,'(7E25.17)') time,y2p(2),mean_xz_km_write(2),sqrt(abs(mean_xz_km_write(2)/(y2p(2)-y2p(1)))/Re)*Re,&
                                   y2p(M2-1),mean_xz_km_write(M2-1),sqrt(abs(mean_xz_km_write(M2-1)/(y2p(M2-1)-y2p(M2)))/Re)*Re
        CALL flush(53)
     END IF
  END IF
  !===========================================================================================================

  !--------------------------------------------------------------------------------------------------------
  !--- <u_i>(t_k) -----------------------------------------------------------------------------------------
  !--- <u_i u_j>(t_k) -------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
  mean_gbl (phase,:,:,:,1:3) = mean_gbl (phase,:,:,:,1:3)*(repetition(phase)-1)
  covar_gbl(phase,:,:,:,1:6) = covar_gbl(phase,:,:,:,1:6)*(repetition(phase)-1)
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           mean_gbl(phase,i,j,k,1) = mean_gbl(phase,i,j,k,1) + work1(i,j,k)
           mean_gbl(phase,i,j,k,2) = mean_gbl(phase,i,j,k,2) + work2(i,j,k)
           mean_gbl(phase,i,j,k,3) = mean_gbl(phase,i,j,k,3) + work3(i,j,k)
           covar_gbl(phase,i,j,k,1) = covar_gbl(phase,i,j,k,1) + work1(i,j,k)*work1(i,j,k)
           covar_gbl(phase,i,j,k,2) = covar_gbl(phase,i,j,k,2) + work2(i,j,k)*work2(i,j,k)
           covar_gbl(phase,i,j,k,3) = covar_gbl(phase,i,j,k,3) + work3(i,j,k)*work3(i,j,k)
           covar_gbl(phase,i,j,k,4) = covar_gbl(phase,i,j,k,4) + work1(i,j,k)*work2(i,j,k)
           covar_gbl(phase,i,j,k,5) = covar_gbl(phase,i,j,k,5) + work1(i,j,k)*work3(i,j,k)
           covar_gbl(phase,i,j,k,6) = covar_gbl(phase,i,j,k,6) + work2(i,j,k)*work3(i,j,k)
        END DO
     END DO
  END DO
  mean_gbl(phase,:,:,:,1:3) = mean_gbl(phase,:,:,:,1:3)/repetition(phase)
  covar_gbl(phase,:,:,:,1:3) = covar_gbl(phase,:,:,:,1:3)/repetition(phase)
  !--------------------------------------------------------------------------------------------------------

  !===========================================================================================================
  !=== integral energy for the whole domain  =================================================================
  !===========================================================================================================
  FKE = 0.
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           FKE(1) = FKE(1) + dx1p(i)*dx2p(j)*dx3p(k)*(mean_gbl(phase,i,j,k,1)**2+mean_gbl(phase,i,j,k,2)**2+&
                                                      mean_gbl(phase,i,j,k,3)**2)
           !--------------------------------------------------------------------------------------------------------
           !--- u'_i u'_i(t_k) -------------------------------------------------------------------------------------
           !--------------------------------------------------------------------------------------------------------
           FKE(2) = FKE(2) + dx1p(i)*dx2p(j)*dx3p(k)*( (work1(i,j,k) - mean_gbl(phase,i,j,k,1))**2 + &
                                                       (work2(i,j,k) - mean_gbl(phase,i,j,k,2))**2 + &
                                                       (work3(i,j,k) - mean_gbl(phase,i,j,k,3))**2 )
        END DO
     END DO
  END DO
  CALL MPI_ALLREDUCE(FKE,FKE_global,2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  FKE_global = FKE_global/2.
  
  TKE = 0.
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           TKE(1) = TKE(1) + dx1p(i)*dx2p(j)*dx3p(k)*(mean_gbl(phase,i,j,k,1)**2+mean_gbl(phase,i,j,k,2)**2+&
                                                      mean_gbl(phase,i,j,k,3)**2)
           !--------------------------------------------------------------------------------------------------------
           !--- <u'_i u'_i>(t_k) -----------------------------------------------------------------------------------
           !--------------------------------------------------------------------------------------------------------
           TKE(2) = TKE(2) + dx1p(i)*dx2p(j)*dx3p(k)*( &
                             (covar_gbl(phase,i,j,k,1) - mean_gbl(phase,i,j,k,1)*mean_gbl(phase,i,j,k,1) ) + &
                             (covar_gbl(phase,i,j,k,2) - mean_gbl(phase,i,j,k,2)*mean_gbl(phase,i,j,k,2) ) + &
                             (covar_gbl(phase,i,j,k,3) - mean_gbl(phase,i,j,k,3)*mean_gbl(phase,i,j,k,3) ) )
        END DO
     END DO
  END DO
  CALL MPI_ALLREDUCE(TKE,TKE_global,2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  TKE_global = TKE_global/2.

  !--- Output ---
  IF (rank == 0) THEN
    WRITE(43,'(5E25.17)') time, FKE_global(1:2),TKE_global(1:2)
    CALL flush(43)
  END IF
 
  !===========================================================================================================
  !=== integral energy for the Kalman window =================================================================
  !===========================================================================================================
  FKE = 0.
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           FKE(1) = FKE(1) + klm_flag(i,j,k)*dx1p(i)*dx2p(j)*dx3p(k)*(mean_gbl(phase,i,j,k,1)**2+mean_gbl(phase,i,j,k,2)**2+&
                                                      mean_gbl(phase,i,j,k,3)**2)
           !--------------------------------------------------------------------------------------------------------
           !--- u'_i u'_i(t_k) -------------------------------------------------------------------------------------
           !--------------------------------------------------------------------------------------------------------
           FKE(2) = FKE(2) + klm_flag(i,j,k)*dx1p(i)*dx2p(j)*dx3p(k)*( (work1(i,j,k) - mean_gbl(phase,i,j,k,1))**2 + &
                                                                       (work2(i,j,k) - mean_gbl(phase,i,j,k,2))**2 + &
                                                                       (work3(i,j,k) - mean_gbl(phase,i,j,k,3))**2 )
        END DO
     END DO
  END DO
  CALL MPI_ALLREDUCE(FKE,FKE_global,2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  FKE_global = FKE_global/2.
  
  TKE = 0.
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           TKE(1) = TKE(1) + klm_flag(i,j,k)*dx1p(i)*dx2p(j)*dx3p(k)*(mean_gbl(phase,i,j,k,1)**2+mean_gbl(phase,i,j,k,2)**2+&
                                                      mean_gbl(phase,i,j,k,3)**2)
           !--------------------------------------------------------------------------------------------------------
           !--- <u'_i u'_i>(t_k) -----------------------------------------------------------------------------------
           !--------------------------------------------------------------------------------------------------------
           TKE(2) = TKE(2) + klm_flag(i,j,k)*dx1p(i)*dx2p(j)*dx3p(k)*( &
                             (covar_gbl(phase,i,j,k,1) - mean_gbl(phase,i,j,k,1)*mean_gbl(phase,i,j,k,1) ) + &
                             (covar_gbl(phase,i,j,k,2) - mean_gbl(phase,i,j,k,2)*mean_gbl(phase,i,j,k,2) ) + &
                             (covar_gbl(phase,i,j,k,3) - mean_gbl(phase,i,j,k,3)*mean_gbl(phase,i,j,k,3) ) )
        END DO
     END DO
  END DO
  CALL MPI_ALLREDUCE(TKE,TKE_global,2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  TKE_global = TKE_global/2.

  !--- Output ---
  IF (rank == 0) THEN
    WRITE(44,'(5E25.17)') time, FKE_global(1:2), TKE_global(1:2)
    CALL flush(44)
  END IF


  write_gain = 0.
  klmn => kalman_first
  do while(associated(klmn))
     do i = 1, klmn%m
        write_gain (klmn%x(i),klmn%y(i),klmn%z(i),1) = klmn%K(3*(i-1)+1,1)
        write_gain (klmn%x(i),klmn%y(i),klmn%z(i),2) = klmn%K(3*(i-1)+2,2)
        write_gain (klmn%x(i),klmn%y(i),klmn%z(i),3) = klmn%K(3*(i-1)+3,3)
        write_gain (klmn%x(i),klmn%y(i),klmn%z(i),4) = klmn%K(3*(i-1)+1,2)
        write_gain (klmn%x(i),klmn%y(i),klmn%z(i),5) = klmn%K(3*(i-1)+1,3)
        write_gain (klmn%x(i),klmn%y(i),klmn%z(i),6) = klmn%K(3*(i-1)+2,3)
     end do
     klmn => klmn%next
  end do

  !===========================================================================================================

  RETURN

  END SUBROUTINE compute_kalman
 
 
  SUBROUTINE write_restart_kalman

  USE mod_dims
  USE mod_vars
  USE mod_lib
  USE mod_inout
  USE usr_func
  USE usr_vars

  ! (basic subroutine)

  IMPLICIT NONE

  !--- write time-dependent data after time integration for next restart ---

  CHARACTER(LEN=3) :: next_restart_char
  CHARACTER(LEN=2) :: phase_char
  TYPE(kalman_t), pointer :: klmn
  INTEGER :: i
  REAL :: write_mean_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  REAL :: write_covar_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)

  IF (write_restart_yes) THEN

     IF (rank == 0) WRITE(*,'(a,i3,a)') 'writing kalman data for restart',restart,' ...'
     CALL num_to_string(3,restart,next_restart_char)
     !===========================================================================================================
     !=== write means_ and covariances_fields for visualization =================================================
     !===========================================================================================================
     do phase = 1,intervals 

        CALL num_to_string(2,phase,phase_char)


        write_mean_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3) = mean_gbl(phase,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
        write_covar_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6) = covar_gbl(phase,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)

        CALL write_hdf('mean_globalX_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'meanX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean_gbl(b1L,b2L,b3L,1))
        CALL write_hdf('mean_globalY_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'meanY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean_gbl(b1L,b2L,b3L,2))
        CALL write_hdf('mean_globalZ_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'meanZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean_gbl(b1L,b2L,b3L,3))
        CALL write_hdf('covar_globalXX_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'covarXX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,1))
        CALL write_hdf('covar_globalYY_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'covarYY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,2))
        CALL write_hdf('covar_globalZZ_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'covarZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,3))
        CALL write_hdf('covar_globalXY_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'covarXY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,4))
        CALL write_hdf('covar_globalXZ_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'covarXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,5))
        CALL write_hdf('covar_globalYZ_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'covarYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,6))


        write_mean_gbl = 0.
        write_covar_gbl = 0.
        klmn => kalman_first
        do while(associated(klmn))
           do i = 1, klmn%m
              write_mean_gbl (klmn%x(i),klmn%y(i),klmn%z(i),1) = klmn%mean (phase,i,1)
              write_mean_gbl (klmn%x(i),klmn%y(i),klmn%z(i),2) = klmn%mean (phase,i,2)
              write_mean_gbl (klmn%x(i),klmn%y(i),klmn%z(i),3) = klmn%mean (phase,i,3)
              write_covar_gbl(klmn%x(i),klmn%y(i),klmn%z(i),1) = klmn%covar(phase,i,1)
              write_covar_gbl(klmn%x(i),klmn%y(i),klmn%z(i),2) = klmn%covar(phase,i,2)
              write_covar_gbl(klmn%x(i),klmn%y(i),klmn%z(i),3) = klmn%covar(phase,i,3)
              write_covar_gbl(klmn%x(i),klmn%y(i),klmn%z(i),4) = klmn%covar(phase,i,4)
              write_covar_gbl(klmn%x(i),klmn%y(i),klmn%z(i),5) = klmn%covar(phase,i,5)
              write_covar_gbl(klmn%x(i),klmn%y(i),klmn%z(i),6) = klmn%covar(phase,i,6)
           end do
           klmn => klmn%next
        end do

        CALL write_hdf('kalmanX_'//phase_char//'_restart.'//next_restart_char,'meanX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean_gbl(b1L,b2L,b3L,1))
        CALL write_hdf('kalmanY_'//phase_char//'_restart.'//next_restart_char,'meanY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean_gbl(b1L,b2L,b3L,2))
        CALL write_hdf('kalmanZ_'//phase_char//'_restart.'//next_restart_char,'meanZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean_gbl(b1L,b2L,b3L,3))
        CALL write_hdf('kalmanXX_'//phase_char//'_restart.'//next_restart_char,'covarXX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,1))
        CALL write_hdf('kalmanYY_'//phase_char//'_restart.'//next_restart_char,'covarYY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,2))
        CALL write_hdf('kalmanZZ_'//phase_char//'_restart.'//next_restart_char,'covarZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,3))
        CALL write_hdf('kalmanXY_'//phase_char//'_restart.'//next_restart_char,'covarXY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,4))
        CALL write_hdf('kalmanXZ_'//phase_char//'_restart.'//next_restart_char,'covarXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,5))
        CALL write_hdf('kalmanYZ_'//phase_char//'_restart.'//next_restart_char,'covarYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,6))
        !===========================================================================================================

     end do 

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

  !--- read time-dependent data from previous restart before time integration starts ---

  IMPLICIT NONE


  CHARACTER(len=80) :: text
  CHARACTER(len=80) :: dummy
  CHARACTER(LEN=2) :: phase_char
  INTEGER           :: ios
  TYPE(kalman_t), pointer :: klmn
  INTEGER :: i
  REAL    :: read_mean_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  REAL    :: read_covar_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)

  IF (rank == 0) WRITE(*,'(a,i3,a)') 'reading kalman data for restart',restart,' ...'

  !========================================================================================================
  !=== neue Restart-Nr. als String fuer File-Namen ========================================================
  !========================================================================================================
  CALL num_to_string(3,restart,restart_char)

  do phase = 1,intervals

     CALL num_to_string(2,phase,phase_char)

     CALL read2_hdf('mean_globalX_phase_'//phase_char//'_kalman_restart.'//restart_char,'meanX',S1p,S2p,S3p,N1p,N2p,N3p,0,read_mean_gbl(b1L,b2L,b3L,1))
     CALL read2_hdf('mean_globalY_phase_'//phase_char//'_kalman_restart.'//restart_char,'meanY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_mean_gbl(b1L,b2L,b3L,2))
     CALL read2_hdf('mean_globalZ_phase_'//phase_char//'_kalman_restart.'//restart_char,'meanZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_mean_gbl(b1L,b2L,b3L,3))
     CALL read2_hdf('covar_globalXX_phase_'//phase_char//'_kalman_restart.'//restart_char,'covarXX',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,1))
     CALL read2_hdf('covar_globalYY_phase_'//phase_char//'_kalman_restart.'//restart_char,'covarYY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,2))
     CALL read2_hdf('covar_globalZZ_phase_'//phase_char//'_kalman_restart.'//restart_char,'covarZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,3))
     CALL read2_hdf('covar_globalXY_phase_'//phase_char//'_kalman_restart.'//restart_char,'covarXY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,4))
     CALL read2_hdf('covar_globalXZ_phase_'//phase_char//'_kalman_restart.'//restart_char,'covarXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,5))
     CALL read2_hdf('covar_globalYZ_phase_'//phase_char//'_kalman_restart.'//restart_char,'covarYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,6))

     mean_gbl(phase,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3) = read_mean_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
     covar_gbl(phase,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6) = read_covar_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)


     !===========================================================================================================
     !=== read means_ and covariances_fields for visualization ==================================================
     !===========================================================================================================
     CALL read2_hdf('kalmanX_'//phase_char//'_restart.'//restart_char,'meanX',S1p,S2p,S3p,N1p,N2p,N3p,0,read_mean_gbl(b1L,b2L,b3L,1))
     CALL read2_hdf('kalmanY_'//phase_char//'_restart.'//restart_char,'meanY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_mean_gbl(b1L,b2L,b3L,2))
     CALL read2_hdf('kalmanZ_'//phase_char//'_restart.'//restart_char,'meanZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_mean_gbl(b1L,b2L,b3L,3))
     CALL read2_hdf('kalmanXX_'//phase_char//'_restart.'//restart_char,'covarXX',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,1))
     CALL read2_hdf('kalmanYY_'//phase_char//'_restart.'//restart_char,'covarYY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,2))
     CALL read2_hdf('kalmanZZ_'//phase_char//'_restart.'//restart_char,'covarZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,3))
     CALL read2_hdf('kalmanXY_'//phase_char//'_restart.'//restart_char,'covarXY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,4))
     CALL read2_hdf('kalmanXZ_'//phase_char//'_restart.'//restart_char,'covarXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,5))
     CALL read2_hdf('kalmanYZ_'//phase_char//'_restart.'//restart_char,'covarYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,6))
     !===========================================================================================================

     klmn => kalman_first
     do while(associated(klmn))
        do i = 1, klmn%m
           klmn%mean(phase,i,1) = read_mean_gbl  (klmn%x(i),klmn%y(i),klmn%z(i),1)
           klmn%mean(phase,i,2) = read_mean_gbl  (klmn%x(i),klmn%y(i),klmn%z(i),2)
           klmn%mean(phase,i,3) = read_mean_gbl  (klmn%x(i),klmn%y(i),klmn%z(i),3)
           klmn%covar(phase,i,1) = read_covar_gbl(klmn%x(i),klmn%y(i),klmn%z(i),1)
           klmn%covar(phase,i,2) = read_covar_gbl(klmn%x(i),klmn%y(i),klmn%z(i),2)
           klmn%covar(phase,i,3) = read_covar_gbl(klmn%x(i),klmn%y(i),klmn%z(i),3)
           klmn%covar(phase,i,4) = read_covar_gbl(klmn%x(i),klmn%y(i),klmn%z(i),4)
           klmn%covar(phase,i,5) = read_covar_gbl(klmn%x(i),klmn%y(i),klmn%z(i),5)
           klmn%covar(phase,i,6) = read_covar_gbl(klmn%x(i),klmn%y(i),klmn%z(i),6)
        end do
        klmn => klmn%next
     end do
     !========================================================================================================
  end do

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
!  real    :: var_20(m, m)

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
