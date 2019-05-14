!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Dario De Marinis, ARTORG CVE, University of Bern (dario.demarinis@artorg.unibe.ch)            *
!* November 2017                                                                                             *
!*************************************************************************************************************

! file stats.f90
! file containing subroutines used for writing and reading statistical data to and from files.

!!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8

  SUBROUTINE open_stats
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE usr_func
  USE usr_vars
  USE mod_lib !added for writing turb_statxz_ with M1 M2 M3

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
  TYPE(stats_group_t), pointer :: stats_group
  TYPE(stats_t), pointer :: stats
  INTEGER :: l,n,i,j,k,ii,jj,kk,id
  INTEGER :: n_procs,strd
  INTEGER, ALLOCATABLE :: n_data_count(:)
  INTEGER, ALLOCATABLE :: i_data(:,:,:)
  REAL :: d_vol

  !added for writing turb_statxz_ with M1 M2 M3===============================================================
  CHARACTER(LEN=3) ::  M1_char
  CHARACTER(LEN=3) ::  M2_char
  CHARACTER(LEN=3) ::  M3_char
  !===========================================================================================================

  NULLIFY(stats_group_first)

  strd = 1
  id = 0

  do while (strd.le.2) 

     id = id + 1
     if (.not.associated(stats_group_first)) then
        allocate(stats_group_first)
        stats_group => stats_group_first
     else
        allocate(stats_group%next)
        stats_group => stats_group%next
     end if
     stats_group%group_id = id
     NULLIFY(stats_group%next)

     !===========================================================================================================
     !=== evaluate n_data for this stats_group ==================================================================
     !===========================================================================================================
     stats_group%n_data = 0
     ALLOCATE(i_data(S1p:N1p,S2p:N2p,S3p:N3p)); i_data = 0
     DO k = S3p, N3p, strd
        DO j = S2p, N2p, strd
           DO i = S1p, N1p, strd
              if (x1p(i).gt.1.0.and.x1p(i).le.2.0) then
              stats_group%n_data = stats_group%n_data + 1
              DO ii = i, i+strd-1
                 DO jj = j, j+strd-1
                    DO kk = k, k+strd-1
                       IF (ii.le.N1p .and. jj.le.N2p .and. kk.le.N3p) i_data(ii,jj,kk) = stats_group%n_data
                    END DO
                 END DO
              END DO
              end if
           END DO
        END DO
     END DO
     n_procs = NB1*NB2*NB3
     ALLOCATE(n_data_count(1:n_procs))
     call MPI_ALLGATHER(stats_group%n_data,1,MPI_INTEGER,n_data_count,1,MPI_INTEGER,COMM_CART,merror)
     stats_group%data_shift = sum(n_data_count(1:rank+1)) - n_data_count(rank+1)
     stats_group%n_data_tot = sum(n_data_count(1:n_procs))
     DEALLOCATE(n_data_count)
     !===========================================================================================================
     
     !===========================================================================================================
     !=== construct matrices d_k, H, R ==========================================================================
     !===========================================================================================================
     NULLIFY(stats_group%stats_first)
     DO n = 1,stats_group%n_data
        if (.NOT.associated(stats_group%stats_first)) then
           allocate(stats_group%stats_first)
           stats => stats_group%stats_first
        else
           allocate(stats%next)
           stats => stats%next
        end if
        stats%i_data = n
        stats%m = count(i_data(S1p:N1p,S2p:N2p,S3p:N3p).eq.n)
        NULLIFY(stats%mean_xyz);   ALLOCATE(stats%mean_xyz(1:3)); stats%mean_xyz = 0.
        NULLIFY(stats%covar_xyz);  ALLOCATE(stats%covar_xyz(1:6)); stats%covar_xyz = 0.
        NULLIFY(stats%mean_xyzt);  ALLOCATE(stats%mean_xyzt(1:3)); stats%mean_xyzt = 0.
        NULLIFY(stats%covar_xyzt); ALLOCATE(stats%covar_xyzt(1:6)); stats%covar_xyzt = 0.
        NULLIFY(stats%x);          ALLOCATE(stats%x(1:stats%m))
        NULLIFY(stats%y);          ALLOCATE(stats%y(1:stats%m))
        NULLIFY(stats%z);          ALLOCATE(stats%z(1:stats%m))
        NULLIFY(stats%wgt);        ALLOCATE(stats%wgt(1:stats%m))
        l = 0
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 if (stats%i_data.eq.i_data(i,j,k)) then
                    l = l + 1
                    stats%x(l) = i
                    stats%y(l) = j
                    stats%z(l) = k
                 end if
              END DO
           END DO
        END DO
        d_vol = 0.
        DO l = 1,stats%m    
           stats%wgt(l) = dx1p(stats%x(l))*dx2p(stats%y(l))*dx3p(stats%z(l))
           d_vol = d_vol + stats%wgt(l)
        END DO
        stats%wgt(1:stats%m) = stats%wgt(1:stats%m)/d_vol
        NULLIFY(stats%next)
     END DO
     DEALLOCATE(i_data)
     !===========================================================================================================







     id = id + 1
     if (.not.associated(stats_group_first)) then
        allocate(stats_group_first)
        stats_group => stats_group_first
     else
        allocate(stats_group%next)
        stats_group => stats_group%next
     end if
     stats_group%group_id = id
     NULLIFY(stats_group%next)

     !===========================================================================================================
     !=== evaluate n_data for this stats_group ==================================================================
     !===========================================================================================================
     stats_group%n_data = 0
     ALLOCATE(i_data(S1p:N1p,S2p:N2p,S3p:N3p)); i_data = 0
     DO k = S3p, N3p, strd
        DO j = S2p, N2p, strd
           if (x2p(j).gt.0.4.and.x2p(j).le.1.6) then
           DO i = S1p, N1p, strd
              stats_group%n_data = stats_group%n_data + 1
              DO ii = i, i+strd-1
                 DO jj = j, j+strd-1
                    DO kk = k, k+strd-1
                       IF (ii.le.N1p .and. jj.le.N2p .and. kk.le.N3p) i_data(ii,jj,kk) = stats_group%n_data
                    END DO
                 END DO
              END DO
           END DO
           end if
        END DO
     END DO
     n_procs = NB1*NB2*NB3
     ALLOCATE(n_data_count(1:n_procs))
     call MPI_ALLGATHER(stats_group%n_data,1,MPI_INTEGER,n_data_count,1,MPI_INTEGER,COMM_CART,merror)
     stats_group%data_shift = sum(n_data_count(1:rank+1)) - n_data_count(rank+1)
     stats_group%n_data_tot = sum(n_data_count(1:n_procs))
     DEALLOCATE(n_data_count)
     !===========================================================================================================
     
     !===========================================================================================================
     !=== construct matrices d_k, H, R ==========================================================================
     !===========================================================================================================
     NULLIFY(stats_group%stats_first)
     DO n = 1,stats_group%n_data
        if (.NOT.associated(stats_group%stats_first)) then
           allocate(stats_group%stats_first)
           stats => stats_group%stats_first
        else
           allocate(stats%next)
           stats => stats%next
        end if
        stats%i_data = n
        stats%m = count(i_data(S1p:N1p,S2p:N2p,S3p:N3p).eq.n)
        NULLIFY(stats%mean_xyz);   ALLOCATE(stats%mean_xyz(1:3)); stats%mean_xyz = 0.
        NULLIFY(stats%covar_xyz);  ALLOCATE(stats%covar_xyz(1:6)); stats%covar_xyz = 0.
        NULLIFY(stats%mean_xyzt);  ALLOCATE(stats%mean_xyzt(1:3)); stats%mean_xyzt = 0.
        NULLIFY(stats%covar_xyzt); ALLOCATE(stats%covar_xyzt(1:6)); stats%covar_xyzt = 0.
        NULLIFY(stats%x);          ALLOCATE(stats%x(1:stats%m))
        NULLIFY(stats%y);          ALLOCATE(stats%y(1:stats%m))
        NULLIFY(stats%z);          ALLOCATE(stats%z(1:stats%m))
        NULLIFY(stats%wgt);        ALLOCATE(stats%wgt(1:stats%m))
        l = 0
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 if (stats%i_data.eq.i_data(i,j,k)) then
                    l = l + 1
                    stats%x(l) = i
                    stats%y(l) = j
                    stats%z(l) = k
                 end if
              END DO
           END DO
        END DO
        d_vol = 0.
        DO l = 1,stats%m    
           stats%wgt(l) = dx1p(stats%x(l))*dx2p(stats%y(l))*dx3p(stats%z(l))
           d_vol = d_vol + stats%wgt(l)
        END DO
        stats%wgt(1:stats%m) = stats%wgt(1:stats%m)/d_vol
        NULLIFY(stats%next)
     END DO
     DEALLOCATE(i_data)
     !===========================================================================================================

     strd = strd + 1
  end do

  num_windows = id !define num_windows in mod_vars for write covariance in xdmf
  !===========================================================================================================
  write_stats_count = 0
  
  !added for writing turb_statxz_ with M1 M2 M3===============================================================
  CALL num_to_string(3,M1,M1_char)
  CALL num_to_string(3,M2,M2_char)
  CALL num_to_string(3,M3,M3_char)
  !============================================================================================================

  IF (rank == 0 .AND. dtime_out_scal /= 0.) THEN
     OPEN(33,FILE='tke_'//restart_char//'.txt',STATUS='UNKNOWN')
     OPEN(23,FILE='turb_statxz_'//restart_char//'_'//M1_char//'x'//M2_char//'x'//M3_char//'.txt',STATUS='UNKNOWN') !added M1 M2 M3 info
  END IF

  END SUBROUTINE open_stats
 
  SUBROUTINE close_stats
  ! (basic subroutine)

  USE mod_dims
  USE mod_vars
  USE usr_func
  USE usr_vars

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
  TYPE(stats_group_t), pointer :: stats_group
  TYPE(stats_t), pointer :: stats

  IF (rank == 0 .AND. dtime_out_scal /= 0.) THEN
     CLOSE(33)
     CLOSE(23)
  END IF

  if(associated(stats_group_first)) then
     stats_group => stats_group_first
     CALL destroy_stats_group(stats_group)
     DEALLOCATE(stats_group_first)
  end if

  IF (rank == 0 .AND. dtime_out_kalm /= 0.) THEN
     CLOSE(43)
     CLOSE(53)
  END IF

  contains

  RECURSIVE SUBROUTINE destroy_stats_group(stats_group)

     implicit none

     TYPE(stats_group_t), pointer, intent(inout) :: stats_group
     TYPE(stats_t), pointer :: stats

     ! ---------------------------------

     if (associated(stats_group%next)) then
        CALL destroy_stats_group(stats_group%next)
        DEALLOCATE(stats_group%next)
     end if

     if(associated(stats_group%stats_first)) then
        stats => stats_group%stats_first
        CALL destroy_stats(stats)
        DEALLOCATE(stats_group%stats_first)
     end if

  END SUBROUTINE destroy_stats_group


  RECURSIVE SUBROUTINE destroy_stats(stats)

     implicit none

     TYPE(stats_t), pointer, intent(inout) :: stats

     ! ---------------------------------

     if (associated(stats%next)) then
        call destroy_stats(stats%next)
        DEALLOCATE(stats%next)
     end if

     DEALLOCATE(stats%x,stats%y,stats%z)
     DEALLOCATE(stats%wgt)
     DEALLOCATE(stats%mean_xyz,stats%covar_xyz)
     DEALLOCATE(stats%mean_xyzt,stats%covar_xyzt)

  END SUBROUTINE destroy_stats


  END SUBROUTINE close_stats


  SUBROUTINE compute_stats
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE mod_lib
  USE mod_diff
  USE mod_inout
  USE usr_func
  USE usr_vars

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  TYPE(stats_group_t), pointer :: stats_group
  TYPE(stats_t), pointer :: stats
  INTEGER :: i, ii
  INTEGER :: j, jj
  INTEGER :: k, kk
  REAL    :: TKE(1:2), TKE_global(1:2), mean_xz_write(1:M2), mean_xz(S2p:N2p), mean_xz_global(S2p:N2p)
  REAL    :: write_mean(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  REAL    :: write_covar(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)
  REAL, ALLOCATABLE :: d_k(:,:), R_k(:,:), R_t(:,:)
  CHARACTER(LEN=8) :: count_char
  CHARACTER(LEN=1) :: id
  CHARACTER(LEN=50) :: write_dir
 
  write_out_scal = .FALSE.
  time_out_scal = time_out_scal + dtime_out_scal

  !===========================================================================================================
  !=== interpolation of velocity to pressure grid points =====================================================
  !===========================================================================================================
  ! vel(:,:,:,i) --> worki(:,:,:)
  IF (task == 1) CALL interpolate_vel(.TRUE.)
  !===========================================================================================================

  !-----------------------------------------------------------------------------------------------------------
  !--- <u>_xz(y) ---------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (wallnorm_dir == 1) THEN
     mean_xz = 0.
     DO j = S2p, N2p
        DO k = S3p, N3p
           DO i = S1p, N1p
              mean_xz(j) = mean_xz(j) + work1(i,j,k)*dx1p(i)*dx3p(k) ! (dx3p = 1. for 2D)
           END DO
        END DO
     END DO
     CALL MPI_ALLREDUCE(mean_xz,mean_xz_global,(N2p-S2p+1),MPI_REAL8,MPI_SUM,COMM_SLICE2,merror)
     mean_xz = mean_xz_global / (L1*L3)
     CALL MPI_ALLGATHERv(mean_xz,(N2p-S2p+1),MPI_REAL8,mean_xz_write,bar2_size,bar2_offset,MPI_REAL8,COMM_BAR2,merror)
     IF (rank == 0) THEN
        ! write into file 'turb_statxz_'//restart_char//'.txt'
        WRITE(23,'(7E25.17)') time,y2p(2),      mean_xz_write(2),   sqrt(abs(mean_xz_write(2)/(y2p(2)-y2p(1)))/Re)*Re,& ! Re_tau = Re * <u_tau(y=0)>_xz
                                   y2p(M2-1),   mean_xz_write(M2-1),sqrt(abs(mean_xz_write(M2-1)/(y2p(M2-1)-y2p(M2)))/Re)*Re !Re_tau = Re * <u_tau(y=M2)>_xz
        CALL flush(23)
     END IF
  END IF
  !===========================================================================================================


  stats_group => stats_group_first
  do while(associated(stats_group))

     !===========================================================================================================
     !=== turbulence statistics (space averaging) ===============================================================
     !===========================================================================================================
     stats => stats_group%stats_first
     do while(associated(stats))
     !-----------------------------------------------------------------------------------------------------------
     !--- <u_i>_xyz(voxel) --------------------------------------------------------------------------------------
     !-----------------------------------------------------------------------------------------------------------
        stats%mean_xyz = 0.
        do i = 1, stats%m
           stats%mean_xyz(1) = stats%mean_xyz(1) + work1(stats%x(i),stats%y(i),stats%z(i))*stats%wgt(i)
           stats%mean_xyz(2) = stats%mean_xyz(2) + work2(stats%x(i),stats%y(i),stats%z(i))*stats%wgt(i)
           stats%mean_xyz(3) = stats%mean_xyz(3) + work3(stats%x(i),stats%y(i),stats%z(i))*stats%wgt(i)
        end do
     !-----------------------------------------------------------------------------------------------------------
     !--- <u_i u_j>_xyz -----------------------------------------------------------------------------------------
     !-----------------------------------------------------------------------------------------------------------
        stats%covar_xyz = 0.
        do i = 1, stats%m
           stats%covar_xyz(1) = stats%covar_xyz(1) + stats%wgt(i)* &
                                (work1(stats%x(i),stats%y(i),stats%z(i))*work1(stats%x(i),stats%y(i),stats%z(i)))
           stats%covar_xyz(2) = stats%covar_xyz(2) + stats%wgt(i)* &
                                (work2(stats%x(i),stats%y(i),stats%z(i))*work2(stats%x(i),stats%y(i),stats%z(i)))
           stats%covar_xyz(3) = stats%covar_xyz(3) + stats%wgt(i)* &
                                (work3(stats%x(i),stats%y(i),stats%z(i))*work3(stats%x(i),stats%y(i),stats%z(i)))
           stats%covar_xyz(4) = stats%covar_xyz(4) + stats%wgt(i)* &
                                (work1(stats%x(i),stats%y(i),stats%z(i))*work2(stats%x(i),stats%y(i),stats%z(i)))
           stats%covar_xyz(5) = stats%covar_xyz(5) + stats%wgt(i)* &
                                (work1(stats%x(i),stats%y(i),stats%z(i))*work3(stats%x(i),stats%y(i),stats%z(i)))
           stats%covar_xyz(6) = stats%covar_xyz(6) + stats%wgt(i)* &
                                (work2(stats%x(i),stats%y(i),stats%z(i))*work3(stats%x(i),stats%y(i),stats%z(i)))
        end do
        stats => stats%next
     end do
     !===========================================================================================================
   
     if (timestep.ne.timestep_old) then
        !===========================================================================================================
        !=== turbulence statistics (space+time averaging) ==========================================================
        !===========================================================================================================
        stats => stats_group%stats_first
        do while(associated(stats))
        !-----------------------------------------------------------------------------------------------------------
        !--- <u_i>_xyzt --------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------------
           stats%mean_xyzt(1:3) = stats%mean_xyzt(1:3)*(time-dtime_out_scal-time_start)
           stats%mean_xyzt(1:3) = stats%mean_xyzt(1:3) + stats%mean_xyz(1:3)*dtime_out_scal
           stats%mean_xyzt(1:3) = stats%mean_xyzt(1:3)/(time-time_start)
        !-----------------------------------------------------------------------------------------------------------
        !--- <u_i u_j>_xyzt ----------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------------
           stats%covar_xyzt(1:6) = stats%covar_xyzt(1:6)*(time-dtime_out_scal-time_start)
           stats%covar_xyzt(1:6) = stats%covar_xyzt(1:6) + stats%covar_xyz(1:6)*dtime_out_scal
           stats%covar_xyzt(1:6) = stats%covar_xyzt(1:6)/(time-time_start)
           stats => stats%next
        end do
        !===========================================================================================================
     end if

     !===========================================================================================================
     !=== write d_k, R_k and R(t) at the current time t=k =======================================================
     !===========================================================================================================
     allocate(d_k(1:stats_group%n_data,1:3)); d_k = 0.
     allocate(R_k(1:stats_group%n_data,1:6)); R_k = 0.
     allocate(R_t(1:stats_group%n_data,1:6)); R_t = 0.
     stats => stats_group%stats_first
     do while(associated(stats))
     !--- <u_i>_xyz ---------------------------------------------------------------------------------------------
        d_k(stats%i_data,1) = stats%mean_xyz(1)
        d_k(stats%i_data,2) = stats%mean_xyz(2)
        d_k(stats%i_data,3) = stats%mean_xyz(3)
     !--- <u'_i u'_j>_xyz ---------------------------------------------------------------------------------------
        R_k(stats%i_data,1) = stats%covar_xyz(1) - stats%mean_xyz(1)*stats%mean_xyz(1) 
        R_k(stats%i_data,2) = stats%covar_xyz(2) - stats%mean_xyz(2)*stats%mean_xyz(2) 
        R_k(stats%i_data,3) = stats%covar_xyz(3) - stats%mean_xyz(3)*stats%mean_xyz(3) 
        R_k(stats%i_data,4) = stats%covar_xyz(4) - stats%mean_xyz(1)*stats%mean_xyz(2) 
        R_k(stats%i_data,5) = stats%covar_xyz(5) - stats%mean_xyz(1)*stats%mean_xyz(3) 
        R_k(stats%i_data,6) = stats%covar_xyz(6) - stats%mean_xyz(2)*stats%mean_xyz(3) 
     !--- <u'_i u'_j>_xyzt --------------------------------------------------------------------------------------
        R_t(stats%i_data,1) = stats%covar_xyzt(1) - stats%mean_xyzt(1)*stats%mean_xyzt(1) 
        R_t(stats%i_data,2) = stats%covar_xyzt(2) - stats%mean_xyzt(2)*stats%mean_xyzt(2) 
        R_t(stats%i_data,3) = stats%covar_xyzt(3) - stats%mean_xyzt(3)*stats%mean_xyzt(3) 
        R_t(stats%i_data,4) = stats%covar_xyzt(4) - stats%mean_xyzt(1)*stats%mean_xyzt(2) 
        R_t(stats%i_data,5) = stats%covar_xyzt(5) - stats%mean_xyzt(1)*stats%mean_xyzt(3) 
        R_t(stats%i_data,6) = stats%covar_xyzt(6) - stats%mean_xyzt(2)*stats%mean_xyzt(3) 
        stats => stats%next
     end do
   
     CALL num_to_string(8,write_stats_count,count_char)
     IF (rank == 0) WRITE(*,'(a,i8,a)') 'writing data stats',write_stats_count,' ...'

     CALL num_to_string(1,stats_group%group_id,id)
     write_dir = './data_'//id//'/'
   
     CALL write_stats_hdf_2D(trim(write_dir)//'d_k.'//count_char,'d_k', &
                           stats_group%n_data_tot,3,1,1,stats_group%n_data,3,stats_group%data_shift,0,d_k) 
     CALL write_stats_hdf_2D(trim(write_dir)//'R_k.'//count_char,'R_k', &
                            stats_group%n_data_tot,6,1,1,stats_group%n_data,6,stats_group%data_shift,0,R_k)
     !CALL write_stats_hdf_2D(trim(write_dir)//'R_t.'//count_char,'R_t', &
     !                        stats_group%n_data_tot,6,1,1,stats_group%n_data,6,stats_group%data_shift,0,R_t)
     !===========================================================================================================

     !===========================================================================================================
     !write fluctuations_fields for visualization
     !===========================================================================================================
     IF (write_out_vect) THEN
        write_mean = 0.
        write_covar = 0.
        stats => stats_group%stats_first
        do while(associated(stats))
           do i = 1,stats%m
              write_mean(stats%x(i),stats%y(i),stats%z(i),1)  = d_k(stats%i_data,1)
              write_mean(stats%x(i),stats%y(i),stats%z(i),2)  = d_k(stats%i_data,2)
              write_mean(stats%x(i),stats%y(i),stats%z(i),3)  = d_k(stats%i_data,3)
              write_covar(stats%x(i),stats%y(i),stats%z(i),1) = R_k(stats%i_data,1)
              write_covar(stats%x(i),stats%y(i),stats%z(i),2) = R_k(stats%i_data,2)
              write_covar(stats%x(i),stats%y(i),stats%z(i),3) = R_k(stats%i_data,3)
              write_covar(stats%x(i),stats%y(i),stats%z(i),4) = R_k(stats%i_data,4)
              write_covar(stats%x(i),stats%y(i),stats%z(i),5) = R_k(stats%i_data,5)
              write_covar(stats%x(i),stats%y(i),stats%z(i),6) = R_k(stats%i_data,6)
           end do 
           stats => stats%next
        end do

        CALL num_to_string(8,write_count,count_char)
        IF (rank == 0) WRITE(*,'(a,i8,a)') 'writing data fields',write_count,' ...'

        CALL num_to_string(1,stats_group%group_id,id)
        write_dir = './data_'//id//'/'
        !    write_hdf(filename,dsetname,
                      !SS1,SS2,SS3,NN1,NN2,NN3,vel_dir,stride,phi)

        CALL write_hdf(trim(write_dir)//'dataX_'//count_char,'meanX', &
                       S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean(b1L,b2L,b3L,1))
        CALL write_hdf(trim(write_dir)//'dataY_'//count_char,'meanY', &
                       S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean(b1L,b2L,b3L,2))
        CALL write_hdf(trim(write_dir)//'dataZ_'//count_char,'meanZ', &
                       S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean(b1L,b2L,b3L,3))
        CALL write_hdf(trim(write_dir)//'dataXX_'//count_char,'covarXX', &
                       S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,1))
        CALL write_hdf(trim(write_dir)//'dataYY_'//count_char,'covarYY', &
                       S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,2))
        CALL write_hdf(trim(write_dir)//'dataZZ_'//count_char,'covarZZ', &
                       S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,3))

        CALL write_hdf(trim(write_dir)//'dataXY_'//count_char,'covarXY', &
                      S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,4))
        CALL write_hdf(trim(write_dir)//'dataXZ_'//count_char,'covarXZ', &
                      S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,5))
        CALL write_hdf(trim(write_dir)//'dataYZ_'//count_char,'covarYZ', &
                      S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar(b1L,b2L,b3L,6))
     END IF

     deallocate(d_k,R_k,R_t)

     stats_group => stats_group%next
  end do
   
  write_stats_count = write_stats_count + 1

  !===========================================================================================================

!  !===========================================================================================================
!  !=== integral energy =======================================================================================
!  !===========================================================================================================
!  TKE = 0.
!  DO k = S3p, N3p
!     DO j = S2p, N2p
!        DO i = S1p, N1p
!           TKE(1) = TKE(1) + dx1p(i)*dx2p(j)*dx3p(k)*(work1(i,j,k)**2 + work2(i,j,k)**2 + work3(i,j,k)**2)
!           TKE(2) = TKE(2) + dx1p(i)*dx2p(j)*dx3p(k)*(write_mean(i,j,k,1)**2+write_mean(i,j,k,2)**2+&
!                                                      write_mean(i,j,k,3)**2)
!        END DO
!     END DO
!  END DO
!  CALL MPI_ALLREDUCE(TKE,TKE_global,2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
!  !--- Kinetic energy density ---
!  TKE = TKE_global / (2.*L1*L2*L3)
!  !--- Absolute kinetic energy ---
!  TKE_global = TKE_global / 2.
!  !--- Output ---
!  IF (rank == 0) THEN
!    WRITE(33,'(5E25.17)') time, TKE_global(1:2), TKE(1:2)
!    CALL flush(33)
!  END IF
!  !===========================================================================================================

  RETURN

  END SUBROUTINE compute_stats
 

  SUBROUTINE write_restart_stats

  USE mod_dims
  USE mod_vars
  USE mod_lib
  USE mod_diff
  USE mod_inout
  USE usr_func
  USE usr_vars

  IMPLICIT NONE

  !--- write time-dependent data after time integration for next restart ---
  ! this is just an example/template

  CHARACTER(LEN=3) :: next_restart_char
  CHARACTER(LEN=1) :: id
  TYPE(stats_group_t), pointer :: stats_group
  TYPE(stats_t), pointer :: stats
  INTEGER :: i
  REAL, ALLOCATABLE :: mean_xyzt (:,:), covar_xyzt(:,:)

  IF (rank == 0) WRITE(*,'(a,i3,a)') 'writing stats data for restart',restart,' ...'

  IF (write_restart_yes) THEN

     stats_group => stats_group_first
     do while(associated(stats_group))

        ALLOCATE(mean_xyzt(1:stats_group%n_data,1:3)); mean_xyzt = 0.
        ALLOCATE(covar_xyzt(1:stats_group%n_data,1:6)); covar_xyzt = 0.

        stats => stats_group%stats_first
        do while(associated(stats))
           mean_xyzt(stats%i_data,1) = stats%mean_xyzt(1)
           mean_xyzt(stats%i_data,2) = stats%mean_xyzt(2)
           mean_xyzt(stats%i_data,3) = stats%mean_xyzt(3)
           covar_xyzt(stats%i_data,1) = stats%covar_xyzt(1)
           covar_xyzt(stats%i_data,2) = stats%covar_xyzt(2)
           covar_xyzt(stats%i_data,3) = stats%covar_xyzt(3)
           covar_xyzt(stats%i_data,4) = stats%covar_xyzt(4)
           covar_xyzt(stats%i_data,5) = stats%covar_xyzt(5)
           covar_xyzt(stats%i_data,6) = stats%covar_xyzt(6)
           stats => stats%next
        end do

        !========================================================================================================
        !=== neue Restart-Nr. als String fuer File-Namen ========================================================
        !========================================================================================================
        CALL num_to_string(3,restart,next_restart_char)
        !========================================================================================================

        CALL num_to_string(1,stats_group%group_id,id)

        !========================================================================================================
        !=== Schreiben ==========================================================================================
        !========================================================================================================
        CALL write_stats_hdf_2D('mean_stats'//id//'_restart.'//next_restart_char, 'mean_stats' ,&
                                stats_group%n_data_tot,3,1,1,stats_group%n_data,3,stats_group%data_shift,0,mean_xyzt)
        CALL write_stats_hdf_2D('covar_stats'//id//'_restart.'//next_restart_char,'covar_stats',&
                                stats_group%n_data_tot,6,1,1,stats_group%n_data,6,stats_group%data_shift,0,covar_xyzt)
     
        DEALLOCATE(mean_xyzt,covar_xyzt)
 
        stats_group => stats_group%next
     end do

  END IF

  RETURN

  END SUBROUTINE write_restart_stats


  ! subroutine that reads data from previous restart before starting time integration 
  SUBROUTINE read_restart_stats

  USE mod_dims
  USE mod_vars
  USE mod_lib
  USE mod_diff
  USE mod_inout
  USE usr_func
  USE usr_vars

  !--- read time-dependent data from previous restart before time integration starts ---

  IMPLICIT NONE
  
  CHARACTER(LEN=1) :: id
  TYPE(stats_group_t), pointer :: stats_group
  TYPE(stats_t), pointer :: stats
  INTEGER :: i,ios
  REAL, ALLOCATABLE :: mean_xyzt (:,:), covar_xyzt(:,:)

  IF (rank == 0) WRITE(*,'(a,i3,a)') 'reading stats data for restart',restart,' ...'

  !========================================================================================================
  !=== neue Restart-Nr. als String fuer File-Namen ========================================================
  !========================================================================================================
  CALL num_to_string(3,restart,restart_char)
  !========================================================================================================

  stats_group => stats_group_first
  do while(associated(stats_group))

     ALLOCATE(mean_xyzt(1:stats_group%n_data,1:3))
     ALLOCATE(covar_xyzt(1:stats_group%n_data,1:6))

     CALL num_to_string(1,stats_group%group_id,id)
     !read the observed data 
     CALL read_stats_hdf_2D('mean_stats'//id//'_restart.'//restart_char, 'mean_stats' , &
                            stats_group%n_data_tot,3,1,1,stats_group%n_data,3,stats_group%data_shift,0,mean_xyzt)
     CALL read_stats_hdf_2D('covar_stats'//id//'_restart.'//restart_char,'covar_stats', &
                            stats_group%n_data_tot,6,1,1,stats_group%n_data,6,stats_group%data_shift,0,covar_xyzt)

     stats => stats_group%stats_first
     do while(associated(stats))
        stats%mean_xyzt(1) = mean_xyzt(stats%i_data,1)
        stats%mean_xyzt(2) = mean_xyzt(stats%i_data,2)
        stats%mean_xyzt(3) = mean_xyzt(stats%i_data,3)
        stats%covar_xyzt(1) = covar_xyzt(stats%i_data,1)
        stats%covar_xyzt(2) = covar_xyzt(stats%i_data,2)
        stats%covar_xyzt(3) = covar_xyzt(stats%i_data,3)
        stats%covar_xyzt(4) = covar_xyzt(stats%i_data,4)
        stats%covar_xyzt(5) = covar_xyzt(stats%i_data,5)
        stats%covar_xyzt(6) = covar_xyzt(stats%i_data,6)
        stats => stats%next
     end do

     DEALLOCATE(mean_xyzt,covar_xyzt)

     stats_group => stats_group%next
  end do

  RETURN
 
  END SUBROUTINE read_restart_stats
