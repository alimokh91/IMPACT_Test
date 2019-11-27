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
  USE MPI
  USE HDF5

  IMPLICIT NONE
 
  INTEGER :: n
  CHARACTER(LEN=2) ::  phs
  CHARACTER(LEN=3) ::  M1_char, M2_char, M3_char
 
  ALLOCATE(mean_gbl(1:intervals,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)); mean_gbl = 0.
  ALLOCATE(covar_gbl(1:intervals,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)); covar_gbl = 0.
 
  ALLOCATE(dtime_phases(1:intervals)); dtime_phases(1:intervals) = dtime_out_scal 
  if (allocated(dtime_kalm_phases)) dtime_phases(1:intervals) = dtime_kalm_phases(1:intervals)
  ALLOCATE(repetition(1:intervals)); repetition = 0

  write_stats_count = 0
  phase = 1
 
  !===========================================================================================================
  !added for writing turb_statxz_ with M1 M2 M3===============================================================
  !===========================================================================================================
  IF ( rank==0 ) THEN
     CALL num_to_string(3,M1,M1_char)
     CALL num_to_string(3,M2,M2_char)
     CALL num_to_string(3,M3,M3_char)
     CALL num_to_string(3,restart,restart_char)
     CALL system('mkdir -p stats_result')
     do n = 1,intervals
        CALL num_to_string(2,n,phs)
        CALL system('mkdir -p stats_result/phase_'//phs )
     end do
     OPEN(23,FILE='./stats_result/turb_statxz_'//restart_char//'_'//M1_char//'x'//M2_char//'x'//M3_char//'.txt',STATUS='UNKNOWN') !added M1 M2 M3 info
     OPEN(33,FILE='./stats_result/tke_domain_'//restart_char//'_'//M1_char//'x'//M2_char//'x'//M3_char//'.txt',STATUS='UNKNOWN')
  END IF
  !=============================================================================================================
 
  CLOSE(10)
  CLOSE(11)
 
  END SUBROUTINE open_stats
 
  SUBROUTINE close_stats
  ! (basic subroutine)
 
  USE mod_dims
  USE mod_vars
  USE usr_func
  USE usr_vars
  USE MPI
 
  IMPLICIT NONE
 
  IF (rank == 0 .AND. dtime_out_scal /= 0.) THEN
    CLOSE(23)
    CLOSE(33)
  END IF
 
  DEALLOCATE(mean_gbl,covar_gbl)
  DEALLOCATE(dtime_phases)
  DEALLOCATE(repetition)

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
  USE MPI
  USE HDF5
 
  IMPLICIT NONE
 
  INTEGER :: i,j,k
  REAL    :: TKE(1:3), TKE_global(1:3)
  REAL    :: mean_xz_write(1:M2), mean_xz(S2p:N2p), mean_xz_global(S2p:N2p)

  phase = mod(write_stats_count,intervals) + 1
  repetition(phase) = repetition(phase) + 1
 
  IF (rank == 0) WRITE(*,'(a,i8,a,i8,a)') 'data stats repetition', repetition(phase), '   for phase', phase,' ...'
  !===========================================================================================================
  !=== interpolation of velocity to pressure grid points =====================================================
  !===========================================================================================================
  ! vel(:,:,:,i) --> worki(:,:,:)
  IF (task == 1) CALL interpolate_vel(.TRUE.)
  !===========================================================================================================

  !===========================================================================================================
  !=== turbulence statistics =================================================================================
  !===========================================================================================================
  !--------------------------------------------------------------------------------------------------------
  !--- <u_i u_j>(t_k) = <u'_i u'_j>(t_k) + <u_i>(t_k) <u_j>(t_k) ------------------------------------------
  !--------------------------------------------------------------------------------------------------------
  covar_gbl(phase,:,:,:,1) = covar_gbl(phase,:,:,:,1) + mean_gbl(phase,:,:,:,1)*mean_gbl(phase,:,:,:,1)
  covar_gbl(phase,:,:,:,2) = covar_gbl(phase,:,:,:,2) + mean_gbl(phase,:,:,:,2)*mean_gbl(phase,:,:,:,2)
  covar_gbl(phase,:,:,:,3) = covar_gbl(phase,:,:,:,3) + mean_gbl(phase,:,:,:,3)*mean_gbl(phase,:,:,:,3)
  covar_gbl(phase,:,:,:,4) = covar_gbl(phase,:,:,:,4) + mean_gbl(phase,:,:,:,1)*mean_gbl(phase,:,:,:,2)
  covar_gbl(phase,:,:,:,5) = covar_gbl(phase,:,:,:,5) + mean_gbl(phase,:,:,:,1)*mean_gbl(phase,:,:,:,3)
  covar_gbl(phase,:,:,:,6) = covar_gbl(phase,:,:,:,6) + mean_gbl(phase,:,:,:,2)*mean_gbl(phase,:,:,:,3)
  !--------------------------------------------------------------------------------------------------------
  !--- <u_i>(t_k) -----------------------------------------------------------------------------------------
  !--- <u_i u_j>(t_k) -------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------
  mean_gbl (phase,:,:,:,1:3) = mean_gbl (phase,:,:,:,1:3)*(repetition(phase)-1)
  covar_gbl(phase,:,:,:,1:6) = covar_gbl(phase,:,:,:,1:6)*(repetition(phase)-1)
  DO k =  b3L,(N3+b3U)
    DO j = b2L,(N2+b2U)
       DO i = b1L,(N1+b1U)
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
  mean_gbl (phase,:,:,:,1:3) = mean_gbl (phase,:,:,:,1:3)/repetition(phase)
  covar_gbl(phase,:,:,:,1:6) = covar_gbl(phase,:,:,:,1:6)/repetition(phase)
  !-----------------------------------------------------------------------------------------------------------
  !--- <u'_i u'_j>(t_k) --------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  covar_gbl(phase,:,:,:,1) = covar_gbl(phase,:,:,:,1) - mean_gbl(phase,:,:,:,1)*mean_gbl(phase,:,:,:,1)
  covar_gbl(phase,:,:,:,2) = covar_gbl(phase,:,:,:,2) - mean_gbl(phase,:,:,:,2)*mean_gbl(phase,:,:,:,2)
  covar_gbl(phase,:,:,:,3) = covar_gbl(phase,:,:,:,3) - mean_gbl(phase,:,:,:,3)*mean_gbl(phase,:,:,:,3)
  covar_gbl(phase,:,:,:,4) = covar_gbl(phase,:,:,:,4) - mean_gbl(phase,:,:,:,1)*mean_gbl(phase,:,:,:,2)
  covar_gbl(phase,:,:,:,5) = covar_gbl(phase,:,:,:,5) - mean_gbl(phase,:,:,:,1)*mean_gbl(phase,:,:,:,3)
  covar_gbl(phase,:,:,:,6) = covar_gbl(phase,:,:,:,6) - mean_gbl(phase,:,:,:,2)*mean_gbl(phase,:,:,:,3)

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
       WRITE(23,'(5E25.17)') time,mean_xz_write(2),   sqrt(abs(mean_xz_write(2)/(y2p(2)-y2p(1)))/Re)*Re,& !Re_tau=Re<u_tau(y=0)>_xz
                                  mean_xz_write(M2-1),sqrt(abs(mean_xz_write(M2-1)/(y2p(M2-1)-y2p(M2)))/Re)*Re !Re_tau=Re<u_tau(y=M2)>_xz
       CALL flush(23)
    END IF
  END IF
  !===========================================================================================================

  !===========================================================================================================
  !=== integral energy for the whole domain  =================================================================
  !===========================================================================================================
  TKE = 0.
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
     !--------------------------------------------------------------------------------------------------------
     !--- <u_i><u_i>(t_k) -----------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
           TKE(1) = TKE(1) + dx1p(i)*dx2p(j)*dx3p(k)*(mean_gbl(phase,i,j,k,1)**2+mean_gbl(phase,i,j,k,2)**2+&
                                                      mean_gbl(phase,i,j,k,3)**2)
     !--------------------------------------------------------------------------------------------------------
     !--- u'_i u'_i(t_k) -------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
           TKE(2) = TKE(2) + dx1p(i)*dx2p(j)*dx3p(k)*( (work1(i,j,k) - mean_gbl(phase,i,j,k,1))**2 + &
                                                       (work2(i,j,k) - mean_gbl(phase,i,j,k,2))**2 + &
                                                       (work3(i,j,k) - mean_gbl(phase,i,j,k,3))**2 )
     !--------------------------------------------------------------------------------------------------------
     !--- <u'_i u'_i>(t_k) -----------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
           TKE(3) = TKE(3) + dx1p(i)*dx2p(j)*dx3p(k)*( covar_gbl(phase,i,j,k,1) + covar_gbl(phase,i,j,k,2) + &
                                                       covar_gbl(phase,i,j,k,3) )
        END DO
     END DO
  END DO
  CALL MPI_ALLREDUCE(TKE,TKE_global,3,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  TKE_global = TKE_global/2.

  !--- Output ---
  IF (rank == 0) THEN
    WRITE(33,'(4E25.17)') time, TKE_global(1:3)
    CALL flush(33)
  END IF
  !===========================================================================================================

  if (write_out_scal) then
     write_stats_count = write_stats_count + 1
     dtime_out_scal = dtime_phases(phase)
     time_out_scal = time_out_scal + dtime_out_scal
     write_out_scal = .FALSE.
  end if
 
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
  USE MPI
  USE HDF5
 
  IMPLICIT NONE
 
  !--- write time-dependent data after time integration for next restart ---
  ! this is just an example/template

  CHARACTER(LEN=3) :: next_restart_char
  CHARACTER(LEN=2) :: phase_char
  INTEGER :: i
  REAL :: write_mean_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  REAL :: write_covar_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)
 
  IF (rank == 0) WRITE(*,'(a,i3,a)') 'writing stats data for restart',restart,' ...'
 
  !========================================================================================================
  !=== neue Restart-Nr. als String fuer File-Namen ========================================================
  !========================================================================================================
  CALL num_to_string(3,restart,next_restart_char)
  !========================================================================================================

  IF (write_restart_yes) THEN

    do i = 1,intervals
 
      CALL num_to_string(2,i,phase_char)
  
      write_mean_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3) = mean_gbl(i,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3) 
      CALL write_hdf('./stats_result/mean_X_phase'//phase_char//'_restart'//next_restart_char,'meanX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean_gbl(b1L,b2L,b3L,1))
      CALL write_hdf('./stats_result/mean_Y_phase'//phase_char//'_restart'//next_restart_char,'meanY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean_gbl(b1L,b2L,b3L,2))
      CALL write_hdf('./stats_result/mean_Z_phase'//phase_char//'_restart'//next_restart_char,'meanZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_mean_gbl(b1L,b2L,b3L,3))

      write_covar_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6) = covar_gbl(i,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6) 
      CALL write_hdf('./stats_result/covar_XX_phase'//phase_char//'_restart'//next_restart_char,'covarXX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,1))
      CALL write_hdf('./stats_result/covar_YY_phase'//phase_char//'_restart'//next_restart_char,'covarYY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,2))
      CALL write_hdf('./stats_result/covar_ZZ_phase'//phase_char//'_restart'//next_restart_char,'covarZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,3))
      CALL write_hdf('./stats_result/covar_XY_phase'//phase_char//'_restart'//next_restart_char,'covarXY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,4))
      CALL write_hdf('./stats_result/covar_XZ_phase'//phase_char//'_restart'//next_restart_char,'covarXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,5))
      CALL write_hdf('./stats_result/covar_YZ_phase'//phase_char//'_restart'//next_restart_char,'covarYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_covar_gbl(b1L,b2L,b3L,6))
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
  USE MPI
  USE HDF5
 
  !--- read time-dependent data from previous restart before time integration starts ---

  IMPLICIT NONE
 
  CHARACTER(LEN=2) :: id,phase_char
  INTEGER :: i
  REAL :: read_mean_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  REAL :: read_covar_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)
 
  IF (rank == 0) WRITE(*,'(a,i3,a)') 'reading stats data for restart',restart,' ...'
 
  !========================================================================================================
  !=== neue Restart-Nr. als String fuer File-Namen ========================================================
  !========================================================================================================
  CALL num_to_string(3,restart,restart_char)
  !========================================================================================================

  do i = 1,intervals
 
     CALL num_to_string(2,i,phase_char)
   
     CALL read2_hdf('./stats_result/mean_X_phase'//phase_char//'_restart'//restart_char,'meanX',S1p,S2p,S3p,N1p,N2p,N3p,0,read_mean_gbl(b1L,b2L,b3L,1))
     CALL read2_hdf('./stats_result/mean_Y_phase'//phase_char//'_restart'//restart_char,'meanY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_mean_gbl(b1L,b2L,b3L,2))
     CALL read2_hdf('./stats_result/mean_Z_phase'//phase_char//'_restart'//restart_char,'meanZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_mean_gbl(b1L,b2L,b3L,3))
     mean_gbl(i,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3) = read_mean_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3) 
 
     CALL read2_hdf('./stats_result/covar_XX_phase'//phase_char//'_restart'//restart_char,'covarXX',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,1))
     CALL read2_hdf('./stats_result/covar_YY_phase'//phase_char//'_restart'//restart_char,'covarYY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,2))
     CALL read2_hdf('./stats_result/covar_ZZ_phase'//phase_char//'_restart'//restart_char,'covarZZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,3))
     CALL read2_hdf('./stats_result/covar_XY_phase'//phase_char//'_restart'//restart_char,'covarXY',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,4))
     CALL read2_hdf('./stats_result/covar_XZ_phase'//phase_char//'_restart'//restart_char,'covarXZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,5))
     CALL read2_hdf('./stats_result/covar_YZ_phase'//phase_char//'_restart'//restart_char,'covarYZ',S1p,S2p,S3p,N1p,N2p,N3p,0,read_covar_gbl(b1L,b2L,b3L,6))
     covar_gbl(i,b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6) = read_covar_gbl(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)

  end do
 
  RETURN
 
  END SUBROUTINE read_restart_stats
