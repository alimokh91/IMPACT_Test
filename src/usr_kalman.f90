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

  IMPLICIT NONE

  INCLUDE 'mpif.h'
 
  TYPE(kalman_t), pointer :: klmn
  INTEGER :: i,j,k,ii,jj,kk
  INTEGER :: m,n,l
  INTEGER :: n_procs,strd
  INTEGER, ALLOCATABLE :: n_data_count(:)
  INTEGER, ALLOCATABLE :: i_data(:,:,:)
  REAL :: d_vol

  strd = 2
  n_data = 0
  ALLOCATE(i_data(S1p:N1p,S2p:N2p,S3p:N3p)); i_data = 0
  DO k = S3p, N3p, strd
     DO j = S2p, N2p, strd
        if (x2p(j).gt.0.4.and.x2p(j).le.1.6) then
        DO i = S1p, N1p, strd
           n_data = n_data + 1
           DO ii = i, i+strd-1
              DO jj = j, j+strd-1
                 DO kk = k, k+strd-1
                    IF (ii.le.N1p .and. jj.le.N2p .and. kk.le.N3p) i_data(ii,jj,kk) = n_data
                 END DO
              END DO
           END DO
        END DO
        end if
     END DO
  END DO

  n_procs = NB1*NB2*NB3
  ALLOCATE(n_data_count(1:n_procs))
  call MPI_ALLGATHER(n_data,1,MPI_INTEGER,n_data_count,1,MPI_INTEGER,COMM_CART,merror)
  data_shift = sum(n_data_count(1:rank+1)) - n_data_count(rank+1)
  n_data_tot = sum(n_data_count(1:n_procs))
  DEALLOCATE(n_data_count)

  !===========================================================================================================
  !=== construct matrices x_f, H, R, P_f =====================================================================
  !===========================================================================================================
  NULLIFY(kalman_first)
  DO n = 1,n_data
     if (.NOT.associated(kalman_first)) then
        allocate(kalman_first)
        klmn => kalman_first
     else
        allocate(klmn%next)
        klmn => klmn%next
     end if
     klmn%i_data = n
     klmn%m = count(i_data(S1p:N1p,S2p:N2p,S3p:N3p).eq.n)
     NULLIFY(klmn%mean);      ALLOCATE(klmn%mean(1:klmn%m,1:3)); klmn%mean = 0.
     NULLIFY(klmn%covar);     ALLOCATE(klmn%covar(1:klmn%m,1:6)); klmn%covar = 0.
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
              if (klmn%i_data.eq.i_data(i,j,k)) then
                 l = l + 1
                 klmn%x(l) = i
                 klmn%y(l) = j
                 klmn%z(l) = k
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
        klmn%covar(l,1:3) = 1.0e6*dtime_out_kalm
        klmn%covar(l,4:6) = 1.0e3*dtime_out_kalm
     end do 
     klmn%obs_oper = klmn%obs_oper/d_vol
     NULLIFY(klmn%next)
  END DO
  !===========================================================================================================

  DEALLOCATE(i_data)

  write_kalm_count = 0
 
  IF (rank == 0 .AND. dtime_out_kalm /= 0.) THEN
     OPEN(43,FILE='tke_kalm_'//restart_char//'.txt',STATUS='UNKNOWN')
     OPEN(53,FILE='turb_statxz_kalm_'//restart_char//'.txt',STATUS='UNKNOWN')
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

  if (associated(kalman_first)) then 
     klmn => kalman_first
     CALL destroy_kalman(klmn)
     DEALLOCATE(kalman_first)
  end if

  IF (rank == 0 .AND. dtime_out_kalm /= 0.) THEN
     CLOSE(43)
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
  INTEGER :: i, ii
  INTEGER :: j, jj
  INTEGER :: k, kk
  INTEGER :: INFO
  REAL    :: mean_xz_km_write(1:M2), mean_xz_km(S2p:N2p), mean_xz_km_global(S2p:N2p)
  REAL    :: m_stats(1:n_data,3),c_stats(1:n_data,6)
  REAL    :: write_mean(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  REAL    :: write_covar(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)
  REAL    :: write_gain(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)
  CHARACTER(LEN=8) :: count_char
  CHARACTER(LEN=50) :: read_dir,write_dir

  write_out_kalm = .FALSE.
  time_out_kalm = time_out_kalm + dtime_out_kalm

  CALL num_to_string(8,write_kalm_count,count_char)

  IF (rank == 0) WRITE(*,'(a,i8,a)') 'writing kalman data fields',write_kalm_count,' ...'

  INFO = 0
  !===========================================================================================================
  !=== interpolation of velocity to pressure grid points =====================================================
  !===========================================================================================================
  ! vel(:,:,:,i) --> worki(:,:,:)
  IF (task == 1) CALL interpolate_vel(.TRUE.)
  !===========================================================================================================

  !===========================================================================================================
  !=== read the observed data ================================================================================
  !===========================================================================================================
  read_dir = './data/'
  !--------------------------------------------------------------------------------------------------------
  CALL read_stats_hdf_2D(trim(read_dir)//'d_k.'//count_char,'d_k',n_data_tot,3,1,1,n_data,3,data_shift,0,m_stats)
  !--------------------------------------------------------------------------------------------------------
  CALL read_stats_hdf_2D(trim(read_dir)//'R_k.'//count_char,'R_k',n_data_tot,6,1,1,n_data,6,data_shift,0,c_stats)
  !CALL read_stats_hdf_2D(trim(read_dir)//'R_t.'//count_char,'R_t',n_data_tot,6,1,1,n_data,6,data_shift,0,c_stats)
  !===========================================================================================================

  !===========================================================================================================
  !=== turbulence statistics =================================================================================
  !===========================================================================================================
  klmn => kalman_first
  do while(associated(klmn))
  !--------------------------------------------------------------------------------------------------------
  !--- <u_i>(t) -------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
     if(write_kalm_count.gt.0) klmn%mean(:,:) = klmn%mean(:,:)*(time-dtime_out_kalm-time_start)
     do i = 1, klmn%m
        klmn%mean(i,1) = klmn%mean(i,1) + work1(klmn%x(i),klmn%y(i),klmn%z(i))*dtime_out_kalm
        klmn%mean(i,2) = klmn%mean(i,2) + work2(klmn%x(i),klmn%y(i),klmn%z(i))*dtime_out_kalm
        klmn%mean(i,3) = klmn%mean(i,3) + work3(klmn%x(i),klmn%y(i),klmn%z(i))*dtime_out_kalm
     end do
     klmn%mean(:,:) = klmn%mean(:,:)/(time-time_start)
  !--------------------------------------------------------------------------------------------------------
  !--- <u_i u_j>(t) ---------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
     if(write_kalm_count.gt.0) klmn%covar(:,:) = klmn%covar(:,:)*(time-dtime_out_kalm-time_start)
     do i = 1, klmn%m
        klmn%covar(i,1) = klmn%covar(i,1) + dtime_out_kalm* &
                          (work1(klmn%x(i),klmn%y(i),klmn%z(i))*work1(klmn%x(i),klmn%y(i),klmn%z(i))) 
        klmn%covar(i,2) = klmn%covar(i,2) + dtime_out_kalm* &
                          (work2(klmn%x(i),klmn%y(i),klmn%z(i))*work2(klmn%x(i),klmn%y(i),klmn%z(i)))
        klmn%covar(i,3) = klmn%covar(i,3) + dtime_out_kalm* &
                          (work3(klmn%x(i),klmn%y(i),klmn%z(i))*work3(klmn%x(i),klmn%y(i),klmn%z(i)))
        klmn%covar(i,4) = klmn%covar(i,4) + dtime_out_kalm* &
                          (work1(klmn%x(i),klmn%y(i),klmn%z(i))*work2(klmn%x(i),klmn%y(i),klmn%z(i)))
        klmn%covar(i,5) = klmn%covar(i,5) + dtime_out_kalm* &
                          (work1(klmn%x(i),klmn%y(i),klmn%z(i))*work3(klmn%x(i),klmn%y(i),klmn%z(i)))
        klmn%covar(i,6) = klmn%covar(i,6) + dtime_out_kalm* &
                          (work2(klmn%x(i),klmn%y(i),klmn%z(i))*work3(klmn%x(i),klmn%y(i),klmn%z(i)))
     end do
     klmn%covar(:,:) = klmn%covar(:,:)/(time-time_start)
  !--------------------------------------------------------------------------------------------------------
     klmn => klmn%next
  end do
  !===========================================================================================================

  if (timestep.ne.timestep_old ) then

  !===========================================================================================================
  !=== Fill the matrices and the vectors involved in the kalman filter algorithm =============================
  !=== apply kalman filtering ================================================================================
  !===========================================================================================================
  klmn => kalman_first
  do while(associated(klmn))
     do i = 1, klmn%m
  !--------------------------------------------------------------------------------------------------------
  !--- u(x,t) -------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
        klmn%muf(3*(i-1)+1) = work1(klmn%x(i),klmn%y(i),klmn%z(i))
        klmn%muf(3*(i-1)+2) = work2(klmn%x(i),klmn%y(i),klmn%z(i))
        klmn%muf(3*(i-1)+3) = work3(klmn%x(i),klmn%y(i),klmn%z(i))
  !--------------------------------------------------------------------------------------------------------
  !--- fill the state vector muf and the covariance matrix pf ---------------------------------------------
  !--- <u'_i u'_j>(t) = <u_i u_j>(t) - <u_i>(t) <u_j>(t) --------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
        klmn%pf(3*(i-1)+1,3*(i-1)+1) = klmn%covar(i,1) - klmn%mean(i,1)*klmn%mean(i,1)
        klmn%pf(3*(i-1)+2,3*(i-1)+2) = klmn%covar(i,2) - klmn%mean(i,2)*klmn%mean(i,2)
        klmn%pf(3*(i-1)+3,3*(i-1)+3) = klmn%covar(i,3) - klmn%mean(i,3)*klmn%mean(i,3)
        klmn%pf(3*(i-1)+1,3*(i-1)+2) = klmn%covar(i,4) - klmn%mean(i,1)*klmn%mean(i,2)
        klmn%pf(3*(i-1)+1,3*(i-1)+3) = klmn%covar(i,5) - klmn%mean(i,1)*klmn%mean(i,3)
        klmn%pf(3*(i-1)+2,3*(i-1)+3) = klmn%covar(i,6) - klmn%mean(i,2)*klmn%mean(i,3)
        klmn%pf(3*(i-1)+2,3*(i-1)+1) = klmn%pf(3*(i-1)+1,3*(i-1)+2)
        klmn%pf(3*(i-1)+3,3*(i-1)+1) = klmn%pf(3*(i-1)+1,3*(i-1)+3)
        klmn%pf(3*(i-1)+3,3*(i-1)+2) = klmn%pf(3*(i-1)+2,3*(i-1)+3)
     end do
  !--------------------------------------------------------------------------------------------------------
  !--- d_k and R_k ----------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
     klmn%obs_data(1) = m_stats(klmn%i_data,1)
     klmn%obs_data(2) = m_stats(klmn%i_data,2)
     klmn%obs_data(3) = m_stats(klmn%i_data,3)
     klmn%obs_covar(1,1) = c_stats(klmn%i_data,1)
     klmn%obs_covar(2,2) = c_stats(klmn%i_data,2)
     klmn%obs_covar(3,3) = c_stats(klmn%i_data,3)
     klmn%obs_covar(1,2) = c_stats(klmn%i_data,4)
     klmn%obs_covar(1,3) = c_stats(klmn%i_data,5)
     klmn%obs_covar(2,3) = c_stats(klmn%i_data,6)
     klmn%obs_covar(2,1) = klmn%obs_covar(1,2)
     klmn%obs_covar(3,1) = klmn%obs_covar(1,3)
     klmn%obs_covar(3,2) = klmn%obs_covar(2,3)
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

  end if

  !===========================================================================================================
  !=== write means_ and covariances_fields for visualization =================================================
  !===========================================================================================================
  IF (write_out_vect) THEN
     write_gain = 0.
     write_mean = 0.
     write_covar = 0.
     klmn => kalman_first
     do while(associated(klmn))
        do i = 1, klmn%m
           write_mean (klmn%x(i),klmn%y(i),klmn%z(i),1) = klmn%mean(i,1)
           write_mean (klmn%x(i),klmn%y(i),klmn%z(i),2) = klmn%mean(i,2)
           write_mean (klmn%x(i),klmn%y(i),klmn%z(i),3) = klmn%mean(i,3)
           write_covar(klmn%x(i),klmn%y(i),klmn%z(i),1) = klmn%covar(i,1)
           write_covar(klmn%x(i),klmn%y(i),klmn%z(i),2) = klmn%covar(i,2)
           write_covar(klmn%x(i),klmn%y(i),klmn%z(i),3) = klmn%covar(i,3)
           write_covar(klmn%x(i),klmn%y(i),klmn%z(i),4) = klmn%covar(i,4)
           write_covar(klmn%x(i),klmn%y(i),klmn%z(i),5) = klmn%covar(i,5)
           write_covar(klmn%x(i),klmn%y(i),klmn%z(i),6) = klmn%covar(i,6)
           write_gain (klmn%x(i),klmn%y(i),klmn%z(i),1) = klmn%K(3*(i-1)+1,1)
           write_gain (klmn%x(i),klmn%y(i),klmn%z(i),2) = klmn%K(3*(i-1)+2,2)
           write_gain (klmn%x(i),klmn%y(i),klmn%z(i),3) = klmn%K(3*(i-1)+3,3)
           write_gain (klmn%x(i),klmn%y(i),klmn%z(i),4) = klmn%K(3*(i-1)+1,2)
           write_gain (klmn%x(i),klmn%y(i),klmn%z(i),5) = klmn%K(3*(i-1)+1,3)
           write_gain (klmn%x(i),klmn%y(i),klmn%z(i),6) = klmn%K(3*(i-1)+2,3)
        end do
        klmn => klmn%next
     end do

     CALL num_to_string(8,write_count,count_char)
     write_dir = './kf_result/'

     CALL write_hdf(trim(write_dir)//'kalmanX_'//count_char,'meanX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean(b1L,b2L,b3L,1))
     CALL write_hdf(trim(write_dir)//'kalmanY_'//count_char,'meanY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean(b1L,b2L,b3L,2))
     CALL write_hdf(trim(write_dir)//'kalmanZ_'//count_char,'meanZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean(b1L,b2L,b3L,3))
     CALL write_hdf(trim(write_dir)//'kalmanXX_'//count_char,'covarXX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,1))
     CALL write_hdf(trim(write_dir)//'kalmanYY_'//count_char,'covarYY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,2))
     CALL write_hdf(trim(write_dir)//'kalmanZZ_'//count_char,'covarZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,3))
     !CALL write_hdf(trim(write_dir)//'kalmanXY_'//count_char,'covarXY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,4))
     !CALL write_hdf(trim(write_dir)//'kalmanXZ_'//count_char,'covarXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,5))
     !CALL write_hdf(trim(write_dir)//'kalmanYZ_'//count_char,'covarYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,6))
     CALL write_hdf(trim(write_dir)//'kalman_gain_XX_'//count_char,'gainXX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_gain(b1L,b2L,b3L,1))
     CALL write_hdf(trim(write_dir)//'kalman_gain_YY_'//count_char,'gainYY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_gain(b1L,b2L,b3L,2))
     CALL write_hdf(trim(write_dir)//'kalman_gain_ZZ_'//count_char,'gainZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_gain(b1L,b2L,b3L,3))
     !CALL write_hdf(trim(write_dir)//'kalman_gain_XY_'//count_char,'gainXY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_gain(b1L,b2L,b3L,4))
     !CALL write_hdf(trim(write_dir)//'kalman_gain_XZ_'//count_char,'gainXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_gain(b1L,b2L,b3L,5))
     !CALL write_hdf(trim(write_dir)//'kalman_gain_YZ_'//count_char,'gainYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_gain(b1L,b2L,b3L,6))
  END IF
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

  write_kalm_count = write_kalm_count + 1

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
  TYPE(kalman_t), pointer :: klmn
  INTEGER :: i
  REAL    :: write_mean(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  REAL    :: write_covar(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)

  IF (write_restart_yes) THEN

   !===========================================================================================================
   !=== write means_ and covariances_fields for visualization =================================================
   !===========================================================================================================
   write_mean = 0.
   write_covar = 0.
   klmn => kalman_first
   do while(associated(klmn))
      do i = 1, klmn%m
         write_mean (klmn%x(i),klmn%y(i),klmn%z(i),1) = klmn%mean (i,1)
         write_mean (klmn%x(i),klmn%y(i),klmn%z(i),2) = klmn%mean (i,2)
         write_mean (klmn%x(i),klmn%y(i),klmn%z(i),3) = klmn%mean (i,3)
         write_covar(klmn%x(i),klmn%y(i),klmn%z(i),1) = klmn%covar(i,1)
         write_covar(klmn%x(i),klmn%y(i),klmn%z(i),2) = klmn%covar(i,2)
         write_covar(klmn%x(i),klmn%y(i),klmn%z(i),3) = klmn%covar(i,3)
         write_covar(klmn%x(i),klmn%y(i),klmn%z(i),4) = klmn%covar(i,4)
         write_covar(klmn%x(i),klmn%y(i),klmn%z(i),5) = klmn%covar(i,5)
         write_covar(klmn%x(i),klmn%y(i),klmn%z(i),6) = klmn%covar(i,6)
      end do
      klmn => klmn%next
   end do

   IF (rank == 0) WRITE(*,'(a,i3,a)') 'writing kalman data for restart',restart,' ...'
   CALL num_to_string(3,restart,next_restart_char)

   CALL write_hdf('kalman_restartX_'//next_restart_char,'meanX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean(b1L,b2L,b3L,1))
   CALL write_hdf('kalman_restartY_'//next_restart_char,'meanY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean(b1L,b2L,b3L,2))
   CALL write_hdf('kalman_restartZ_'//next_restart_char,'meanZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean(b1L,b2L,b3L,3))
   CALL write_hdf('kalman_restartXX_'//next_restart_char,'covarXX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,1))
   CALL write_hdf('kalman_restartYY_'//next_restart_char,'covarYY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,2))
   CALL write_hdf('kalman_restartZZ_'//next_restart_char,'covarZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,3))
   CALL write_hdf('kalman_restartXY_'//next_restart_char,'covarXY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,4))
   CALL write_hdf('kalman_restartXZ_'//next_restart_char,'covarXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,5))
   CALL write_hdf('kalman_restartYZ_'//next_restart_char,'covarYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,6))
   !===========================================================================================================

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
  INTEGER           :: ios
  TYPE(kalman_t), pointer :: klmn
  INTEGER :: i
  REAL    :: write_mean(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  REAL    :: write_covar(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)

  IF (rank == 0) WRITE(*,'(a,i3,a)') 'reading kalman data for restart',restart,' ...'

  !========================================================================================================
  !=== neue Restart-Nr. als String fuer File-Namen ========================================================
  !========================================================================================================
  CALL num_to_string(3,restart,restart_char)

  !===========================================================================================================
  !=== read means_ and covariances_fields for visualization ==================================================
  !===========================================================================================================
  CALL read2_hdf('kalman_restartX_'//restart_char,'meanX',S1p,S2p,S3p,N1p,N2p,N3p,0,write_mean(b1L,b2L,b3L,1))
  CALL read2_hdf('kalman_restartY_'//restart_char,'meanY',S1p,S2p,S3p,N1p,N2p,N3p,0,write_mean(b1L,b2L,b3L,2))
  CALL read2_hdf('kalman_restartZ_'//restart_char,'meanZ',S1p,S2p,S3p,N1p,N2p,N3p,0,write_mean(b1L,b2L,b3L,3))
  CALL read2_hdf('kalman_restartXX_'//restart_char,'covarXX',S1p,S2p,S3p,N1p,N2p,N3p,0,write_covar(b1L,b2L,b3L,1))
  CALL read2_hdf('kalman_restartYY_'//restart_char,'covarYY',S1p,S2p,S3p,N1p,N2p,N3p,0,write_covar(b1L,b2L,b3L,2))
  CALL read2_hdf('kalman_restartZZ_'//restart_char,'covarZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,write_covar(b1L,b2L,b3L,3))
  CALL read2_hdf('kalman_restartXY_'//restart_char,'covarXY',S1p,S2p,S3p,N1p,N2p,N3p,0,write_covar(b1L,b2L,b3L,4))
  CALL read2_hdf('kalman_restartXZ_'//restart_char,'covarXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,write_covar(b1L,b2L,b3L,5))
  CALL read2_hdf('kalman_restartYZ_'//restart_char,'covarYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,write_covar(b1L,b2L,b3L,6))
  !===========================================================================================================

  klmn => kalman_first
  do while(associated(klmn))
     do i = 1, klmn%m
        klmn%mean(i,1) = write_mean (klmn%x(i),klmn%y(i),klmn%z(i),1)
        klmn%mean(i,2) = write_mean (klmn%x(i),klmn%y(i),klmn%z(i),2)
        klmn%mean(i,3) = write_mean (klmn%x(i),klmn%y(i),klmn%z(i),3)
        klmn%covar(i,1) = write_covar(klmn%x(i),klmn%y(i),klmn%z(i),1)
        klmn%covar(i,2) = write_covar(klmn%x(i),klmn%y(i),klmn%z(i),2)
        klmn%covar(i,3) = write_covar(klmn%x(i),klmn%y(i),klmn%z(i),3)
        klmn%covar(i,4) = write_covar(klmn%x(i),klmn%y(i),klmn%z(i),4)
        klmn%covar(i,5) = write_covar(klmn%x(i),klmn%y(i),klmn%z(i),5)
        klmn%covar(i,6) = write_covar(klmn%x(i),klmn%y(i),klmn%z(i),6)
     end do
     klmn => klmn%next
  end do
  !========================================================================================================

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
