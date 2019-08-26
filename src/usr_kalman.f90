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

  SUBROUTINE open_kalman
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE usr_func
  USE usr_vars
  USE mod_lib
  USE mod_inout

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  !TYPE(DistHPCPredictMRI), pointer :: mri
  TYPE(HPCPredictMRI), pointer :: mri
  TYPE(kalman_t), pointer :: klmn
  INTEGER :: i,j,k,ii,jj,kk
  INTEGER :: m,n,l
  INTEGER :: n_procs
  INTEGER, ALLOCATABLE :: n_data_count(:)
  INTEGER, ALLOCATABLE :: i_data(:,:,:)
  INTEGER, ALLOCATABLE :: flag_data(:,:,:)
  REAL :: d_vol
  REAL :: dist,dist_min
  !added for writing turb_statxz_ with M1 M2 M3 ================================================================
  CHARACTER(LEN=2) ::  phs
  CHARACTER(LEN=3) ::  M1_char
  CHARACTER(LEN=3) ::  M2_char
  CHARACTER(LEN=3) ::  M3_char
  CHARACTER(LEN=50) :: read_dir,write_dir

  ALLOCATE(mean_gbl(1:intervals,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)); mean_gbl = 0. 
  ALLOCATE(covar_gbl(1:intervals,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)); covar_gbl = 0. 
  ALLOCATE(write_mean (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)); write_mean = 0.
  ALLOCATE(write_fluct (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)); write_fluct = 0.
  ALLOCATE(write_covar(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)); write_covar = 0.
  ALLOCATE(write_gain (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)); write_gain = 0.

  ALLOCATE(dtime_kalm_phases(1:intervals)); dtime_kalm_phases = 0.

  dtime_kalm_phases(1:12) = 0.01
  dtime_kalm_phases(13) = 0.03
  dtime_kalm_phases(14:19) = 0.05
  dtime_kalm_phases(20:21) = 0.1
  dtime_kalm_phases(22) = 6./7. - 0.7 + 0.05

  n_data_tot = 906780
  !===========================================================================================================
  !=== read the observed data (index, x, y, z) ===============================================================
  !===========================================================================================================
  read_dir = './data/Results_HDF5/'
  !--------------------------------------------------------------------------------------------------------
  NULLIFY(mri_first); allocate(mri_first)
  mri => mri_first
  !CALL mr_io_read_parallel_hpcpredict(MPI_COMM_WORLD, MPI_INFO_NULL,trim(read_dir)//'bern_experimental_dataset_hpc_predict_mri.h5',mri)
  CALL mr_io_read_hpcpredict(trim(read_dir)//'bern_experimental_dataset_hpc_predict_mri.h5',mri)
  mri%x_coordinates = mri%x_coordinates/1000.
  mri%y_coordinates = mri%y_coordinates/1000.
  mri%z_coordinates = mri%z_coordinates/1000.
  mri%velocity_mean = mri%velocity_mean/U_ref
  mri%velocity_cov  = mri%velocity_cov/(U_ref*U_ref)

  ALLOCATE(i_data(S1p:N1p,S2p:N2p,S3p:N3p)); i_data = 0
  ALLOCATE(flag_data(S1p:N1p,S2p:N2p,S3p:N3p)); flag_data = 0

  DO k = S3p, N3p
     if ( x3p(k).gt.-17.778965e-3 .and. x3p(k).lt.17.601346e-3 ) then
     DO j = S2p, N2p
        if ( x2p(j).gt.-17.923005e-3 .and. x2p(j).lt.35.786866e-3 ) then
        DO i = S1p, N1p
           if ( x1p(i).gt.-17.811665e-3 .and. x1p(i).lt.17.994916e-3 ) then
           dist_min = 10000.0
           DO kk = 1,mri%z_dim
              DO jj = 1,mri%y_dim
                 DO ii = 1,mri%x_dim
                    dist = (mri%x_coordinates(ii)-x1p(i))**2 + &
                           (mri%y_coordinates(jj)-x2p(j))**2 + (mri%z_coordinates(kk)-x3p(k))**2
                    if (dist.lt.dist_min) then
                       i_data(i,j,k,1) = ii
                       i_data(i,j,k,2) = jj
                       i_data(i,j,k,3) = kk
                       dist_min = min(dist_min,dist)
                    end if
                 END DO
              END DO
           END DO
           end if
        END DO 
        end if
     END DO
     end if
  END DO

  !===========================================================================================================
  !=== construct matrices x_f, H, R, P_f =====================================================================
  !===========================================================================================================
  NULLIFY(kalman_first)
  DO kk = 1,mri%z_dim
     DO jj = 1,mri%y_dim
        DO ii = 1,mri%x_dim
           l = count(i_data(S1p:N1p,S2p:N2p,S3p:N3p,1).eq.ii .and. &
                     i_data(S1p:N1p,S2p:N2p,S3p:N3p,2).eq.jj .and. &
                     i_data(S1p:N1p,S2p:N2p,S3p:N3p,3).eq.kk)
           if (l.ne.0) then
              if (.NOT.associated(kalman_first)) then
                 allocate(kalman_first)
                 klmn => kalman_first
              else
                 allocate(klmn%next)
                 klmn => klmn%next
              end if
              klmn%i_data = n
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
              DO k = S3p, N3p
                 DO j = S2p, N2p
                    DO i = S1p, N1p
                       if (klmn%i_data.eq.i_data(i,j,k,1) .and. klmn%j_data.eq.i_data(i,j,k,2) .and. &
                           klmn%k_data.eq.i_data(i,j,k,3)) then
                          l = l + 1
                          klmn%x(l) = i
                          klmn%y(l) = j
                          klmn%z(l) = k
                          klmn%flg = flag_data(i,j,k)
                       end if
                    END DO
                 END DO
              END DO
              d_vol = 0.
              do l = 1,klmn%m
                 do i = 1,3
                    klmn%obs_oper(i,3*(l-1)+i) = dx1p(klmn%x(l))*dx2p(klmn%y(l))*dx3p(klmn%z(l))
                 end do
                 d_vol = d_vol + dx1p(klmn%x(l))*dx2p(klmn%y(l))*dx3p(klmn%z(l)) ! (dx3p = 1. for 2D)
              end do
              klmn%covar(1:intervals,1:klmn%m,1:3) = 1.0e2
              klmn%covar(1:intervals,1:klmn%m,4:6) = 1.0e1
              klmn%obs_oper = klmn%obs_oper/d_vol
              NULLIFY(klmn%next)
           end if
        END DO
     END DO
  END DO
  !===========================================================================================================

  DEALLOCATE(i_data,flag_data)

  write_kalm_count = -1

  phase = intervals
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

  DEALLOCATE(repetition)
  DEALLOCATE(mean_gbl,covar_gbl,write_mean,write_fluct,write_covar,write_gain)
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
  CHARACTER(LEN=8) :: count_char
  CHARACTER(LEN=2) :: phs
  CHARACTER(LEN=50) :: read_dir,write_dir

  if (time.ge.time_out_kalm) then
     write_kalm_count = write_kalm_count + 1
     phase = mod(write_kalm_count,intervals) + 1
     dtime_out_kalm = dtime_kalm_phases(phase)
     time_out_kalm = time_out_kalm + dtime_out_kalm
     dtime_out_vect = dtime_out_kalm
  end if

  IF (rank == 0) WRITE(*,'(a,i8,a)') 'writing kalman data fields',write_kalm_count,' ...'
  CALL num_to_string(8,write_kalm_count,count_char)

  phase = mod(write_kalm_count,intervals) + 1
  repetition(phase) = repetition(phase) + 1

  INFO = 0
  !===========================================================================================================
  !=== interpolation of velocity to pressure grid points =====================================================
  !===========================================================================================================
  ! vel(:,:,:,i) --> worki(:,:,:)
  IF (task == 1) CALL interpolate_vel(.TRUE.)
  !===========================================================================================================

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
     klmn%obs_data(1:3)      = mri%velocity_mean(1:3,    phase,klmn%i_data,klmn%j_data,klmn%k_data)
     klmn%obs_covar(1:3,1:3) = mri%velocity_cov (1:3,1:3,phase,klmn%i_data,klmn%j_data,klmn%k_data)
     klmn%obs_covar(1:3,1:3) = (1-klmn%flg)*klmn%obs_covar(1:3,1:3) + klmn%flg*100.*klmn%obs_covar(1:3,1:3)
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
  !--------------------------------------------------------------------------------------------------------
  !--- <u'_i u'_i>(t_k) -----------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
  write_fluct(:,:,:,1) = covar_gbl(phase,:,:,:,1) - mean_gbl(phase,:,:,:,1)*mean_gbl(phase,:,:,:,1)
  write_fluct(:,:,:,2) = covar_gbl(phase,:,:,:,2) - mean_gbl(phase,:,:,:,2)*mean_gbl(phase,:,:,:,2)
  write_fluct(:,:,:,3) = covar_gbl(phase,:,:,:,3) - mean_gbl(phase,:,:,:,3)*mean_gbl(phase,:,:,:,3)
  TKE = 0.
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           TKE(1) = TKE(1) + dx1p(i)*dx2p(j)*dx3p(k)*(mean_gbl(phase,i,j,k,1)**2+mean_gbl(phase,i,j,k,2)**2+&
                                                      mean_gbl(phase,i,j,k,3)**2)
           TKE(2) = TKE(2) + dx1p(i)*dx2p(j)*dx3p(k)*(write_fluct(i,j,k,1)+write_fluct(i,j,k,2)+&
                                                      write_fluct(i,j,k,3))
        END DO
     END DO
  END DO
  !--------------------------------------------------------------------------------------------------------
  !--- u'_i(t_k) ------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
  write_fluct(:,:,:,1) = work1(:,:,:) - mean_gbl(phase,:,:,:,1)
  write_fluct(:,:,:,2) = work2(:,:,:) - mean_gbl(phase,:,:,:,2)
  write_fluct(:,:,:,3) = work3(:,:,:) - mean_gbl(phase,:,:,:,3)
  FKE = 0.
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           FKE(1) = FKE(1) + dx1p(i)*dx2p(j)*dx3p(k)*(mean_gbl(phase,i,j,k,1)**2+mean_gbl(phase,i,j,k,2)**2+&
                                                      mean_gbl(phase,i,j,k,3)**2)
           FKE(2) = FKE(2) + dx1p(i)*dx2p(j)*dx3p(k)*(write_fluct(i,j,k,1)**2+write_fluct(i,j,k,2)**2+&
                                                      write_fluct(i,j,k,3)**2)
        END DO
     END DO
  END DO
  
  CALL MPI_ALLREDUCE(FKE,FKE_global,2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  CALL MPI_ALLREDUCE(TKE,TKE_global,2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  !--- Absolute kinetic energy ---
  TKE_global = TKE_global/2.
  FKE_global = FKE_global/2.
  !--- Output ---
  IF (rank == 0) THEN
    WRITE(43,'(5E25.17)') time, FKE_global(1:2),TKE_global(1:2)
    CALL flush(43)
  END IF
 
  !===========================================================================================================
  !=== integral energy for the Kalman window =================================================================
  !===========================================================================================================
  write_mean = 0.
  write_fluct = 0.
  write_gain = 0.
  klmn => kalman_first
  do while(associated(klmn))
     do i = 1, klmn%m
        write_mean (klmn%x(i),klmn%y(i),klmn%z(i),1) = mean_gbl(phase,klmn%x(i),klmn%y(i),klmn%z(i),1)
        write_mean (klmn%x(i),klmn%y(i),klmn%z(i),2) = mean_gbl(phase,klmn%x(i),klmn%y(i),klmn%z(i),2)
        write_mean (klmn%x(i),klmn%y(i),klmn%z(i),3) = mean_gbl(phase,klmn%x(i),klmn%y(i),klmn%z(i),3)
        write_fluct(klmn%x(i),klmn%y(i),klmn%z(i),1) = work1(klmn%x(i),klmn%y(i),klmn%z(i)) - write_mean(klmn%x(i),klmn%y(i),klmn%z(i),1)
        write_fluct(klmn%x(i),klmn%y(i),klmn%z(i),2) = work2(klmn%x(i),klmn%y(i),klmn%z(i)) - write_mean(klmn%x(i),klmn%y(i),klmn%z(i),2)
        write_fluct(klmn%x(i),klmn%y(i),klmn%z(i),3) = work3(klmn%x(i),klmn%y(i),klmn%z(i)) - write_mean(klmn%x(i),klmn%y(i),klmn%z(i),3)
        write_gain (klmn%x(i),klmn%y(i),klmn%z(i),1) = klmn%K(3*(i-1)+1,1)
        write_gain (klmn%x(i),klmn%y(i),klmn%z(i),2) = klmn%K(3*(i-1)+2,2)
        write_gain (klmn%x(i),klmn%y(i),klmn%z(i),3) = klmn%K(3*(i-1)+3,3)
        write_gain (klmn%x(i),klmn%y(i),klmn%z(i),4) = klmn%K(3*(i-1)+1,2)
        write_gain (klmn%x(i),klmn%y(i),klmn%z(i),5) = klmn%K(3*(i-1)+1,3)
        write_gain (klmn%x(i),klmn%y(i),klmn%z(i),6) = klmn%K(3*(i-1)+2,3)
     end do
     klmn => klmn%next
  end do

  FKE = 0.
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           FKE(1) = FKE(1) + dx1p(i)*dx2p(j)*dx3p(k)*(write_mean(i,j,k,1)**2  + write_mean(i,j,k,2)**2  + write_mean(i,j,k,3)**2)
           FKE(2) = FKE(2) + dx1p(i)*dx2p(j)*dx3p(k)*(write_fluct(i,j,k,1)**2 + write_fluct(i,j,k,2)**2 + write_fluct(i,j,k,3)**2)
        END DO
     END DO
  END DO
  CALL MPI_ALLREDUCE(FKE,FKE_global,2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  !--- Absolute kinetic energy ---
  FKE_global = FKE_global /2.

  write_fluct = 0.
  klmn => kalman_first
  do while(associated(klmn))
     do i = 1, klmn%m
        write_fluct(klmn%x(i),klmn%y(i),klmn%z(i),1) = covar_gbl(phase,klmn%x(i),klmn%y(i),klmn%z(i),1) - &
                                                       mean_gbl(phase,klmn%x(i),klmn%y(i),klmn%z(i),1)*mean_gbl(phase,klmn%x(i),klmn%y(i),klmn%z(i),1)
        write_fluct(klmn%x(i),klmn%y(i),klmn%z(i),2) = covar_gbl(phase,klmn%x(i),klmn%y(i),klmn%z(i),2) - &
                                                       mean_gbl(phase,klmn%x(i),klmn%y(i),klmn%z(i),2)*mean_gbl(phase,klmn%x(i),klmn%y(i),klmn%z(i),2)
        write_fluct(klmn%x(i),klmn%y(i),klmn%z(i),3) = covar_gbl(phase,klmn%x(i),klmn%y(i),klmn%z(i),3) - &
                                                       mean_gbl(phase,klmn%x(i),klmn%y(i),klmn%z(i),3)*mean_gbl(phase,klmn%x(i),klmn%y(i),klmn%z(i),3)
     end do
     klmn => klmn%next
  end do

  TKE = 0.
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           TKE(1) = TKE(1) + dx1p(i)*dx2p(j)*dx3p(k)*(write_mean(i,j,k,1)**2  + write_mean(i,j,k,2)**2  + write_mean(i,j,k,3)**2)
           TKE(2) = TKE(2) + dx1p(i)*dx2p(j)*dx3p(k)*(write_fluct(i,j,k,1) + write_fluct(i,j,k,2) + write_fluct(i,j,k,3))
        END DO
     END DO
  END DO
  CALL MPI_ALLREDUCE(TKE,TKE_global,2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  !--- Absolute kinetic energy ---
  TKE_global = TKE_global /2.


  !--- Output ---
  IF (rank == 0) THEN
    WRITE(44,'(5E25.17)') time, FKE_global(1:2), TKE_global(1:2)
    CALL flush(44)
  END IF
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
        write_mean = 0.
        write_covar = 0.
        klmn => kalman_first
        do while(associated(klmn))
           do i = 1, klmn%m
              write_mean (klmn%x(i),klmn%y(i),klmn%z(i),1) = klmn%mean (phase,i,1)
              write_mean (klmn%x(i),klmn%y(i),klmn%z(i),2) = klmn%mean (phase,i,2)
              write_mean (klmn%x(i),klmn%y(i),klmn%z(i),3) = klmn%mean (phase,i,3)
              write_covar(klmn%x(i),klmn%y(i),klmn%z(i),1) = klmn%covar(phase,i,1)
              write_covar(klmn%x(i),klmn%y(i),klmn%z(i),2) = klmn%covar(phase,i,2)
              write_covar(klmn%x(i),klmn%y(i),klmn%z(i),3) = klmn%covar(phase,i,3)
              write_covar(klmn%x(i),klmn%y(i),klmn%z(i),4) = klmn%covar(phase,i,4)
              write_covar(klmn%x(i),klmn%y(i),klmn%z(i),5) = klmn%covar(phase,i,5)
              write_covar(klmn%x(i),klmn%y(i),klmn%z(i),6) = klmn%covar(phase,i,6)
           end do
           klmn => klmn%next
        end do

        write_mean_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3) = mean_gbl(phase,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
        write_covar_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6) = covar_gbl(phase,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)

        CALL num_to_string(2,phase,phase_char)

        CALL write_hdf('mean_globalX_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'meanX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean_gbl(b1L,b2L,b3L,1))
        CALL write_hdf('mean_globalY_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'meanY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean_gbl(b1L,b2L,b3L,2))
        CALL write_hdf('mean_globalZ_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'meanZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean_gbl(b1L,b2L,b3L,3))
        CALL write_hdf('covar_globalXX_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'covarXX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,1))
        CALL write_hdf('covar_globalYY_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'covarYY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,2))
        CALL write_hdf('covar_globalZZ_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'covarZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,3))
        CALL write_hdf('covar_globalXY_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'covarXY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,4))
        CALL write_hdf('covar_globalXZ_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'covarXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,5))
        CALL write_hdf('covar_globalYZ_phase_'//phase_char//'_kalman_restart.'//next_restart_char,'covarYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,6))

        CALL write_hdf('kalmanX_'//phase_char//'_restart.'//next_restart_char,'meanX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean(b1L,b2L,b3L,1))
        CALL write_hdf('kalmanY_'//phase_char//'_restart.'//next_restart_char,'meanY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean(b1L,b2L,b3L,2))
        CALL write_hdf('kalmanZ_'//phase_char//'_restart.'//next_restart_char,'meanZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean(b1L,b2L,b3L,3))
        CALL write_hdf('kalmanXX_'//phase_char//'_restart.'//next_restart_char,'covarXX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,1))
        CALL write_hdf('kalmanYY_'//phase_char//'_restart.'//next_restart_char,'covarYY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,2))
        CALL write_hdf('kalmanZZ_'//phase_char//'_restart.'//next_restart_char,'covarZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,3))
        CALL write_hdf('kalmanXY_'//phase_char//'_restart.'//next_restart_char,'covarXY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,4))
        CALL write_hdf('kalmanXZ_'//phase_char//'_restart.'//next_restart_char,'covarXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,5))
        CALL write_hdf('kalmanYZ_'//phase_char//'_restart.'//next_restart_char,'covarYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,6))
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
  REAL    :: read_mean(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  REAL    :: read_covar(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)

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
     CALL read2_hdf('kalmanX_'//phase_char//'_restart.'//restart_char,'meanX',S1p,S2p,S3p,N1p,N2p,N3p,0,read_mean(b1L,b2L,b3L,1))
     CALL read2_hdf('kalmanY_'//phase_char//'_restart.'//restart_char,'meanY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_mean(b1L,b2L,b3L,2))
     CALL read2_hdf('kalmanZ_'//phase_char//'_restart.'//restart_char,'meanZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_mean(b1L,b2L,b3L,3))
     CALL read2_hdf('kalmanXX_'//phase_char//'_restart.'//restart_char,'covarXX',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar(b1L,b2L,b3L,1))
     CALL read2_hdf('kalmanYY_'//phase_char//'_restart.'//restart_char,'covarYY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar(b1L,b2L,b3L,2))
     CALL read2_hdf('kalmanZZ_'//phase_char//'_restart.'//restart_char,'covarZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar(b1L,b2L,b3L,3))
     CALL read2_hdf('kalmanXY_'//phase_char//'_restart.'//restart_char,'covarXY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar(b1L,b2L,b3L,4))
     CALL read2_hdf('kalmanXZ_'//phase_char//'_restart.'//restart_char,'covarXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar(b1L,b2L,b3L,5))
     CALL read2_hdf('kalmanYZ_'//phase_char//'_restart.'//restart_char,'covarYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar(b1L,b2L,b3L,6))
     !===========================================================================================================

     klmn => kalman_first
     do while(associated(klmn))
        do i = 1, klmn%m
           klmn%mean(phase,i,1) = read_mean (klmn%x(i),klmn%y(i),klmn%z(i),1)
           klmn%mean(phase,i,2) = read_mean (klmn%x(i),klmn%y(i),klmn%z(i),2)
           klmn%mean(phase,i,3) = read_mean (klmn%x(i),klmn%y(i),klmn%z(i),3)
           klmn%covar(phase,i,1) = read_covar(klmn%x(i),klmn%y(i),klmn%z(i),1)
           klmn%covar(phase,i,2) = read_covar(klmn%x(i),klmn%y(i),klmn%z(i),2)
           klmn%covar(phase,i,3) = read_covar(klmn%x(i),klmn%y(i),klmn%z(i),3)
           klmn%covar(phase,i,4) = read_covar(klmn%x(i),klmn%y(i),klmn%z(i),4)
           klmn%covar(phase,i,5) = read_covar(klmn%x(i),klmn%y(i),klmn%z(i),5)
           klmn%covar(phase,i,6) = read_covar(klmn%x(i),klmn%y(i),klmn%z(i),6)
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
