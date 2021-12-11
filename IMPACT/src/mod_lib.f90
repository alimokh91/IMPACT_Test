!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch                    *
!*************************************************************************************************************

!> @file mod_lib.f90
!! routine library file

!> Subroutine library module. It uses modules mod_dims, mod_vars, mod_exchange, mod_diff and mod_coeffs.
MODULE mod_lib
  
  
  USE mod_dims
  USE mod_vars
  USE mod_exchange
  USE mod_diff
  USE mod_coeffs
  !USE mpi 
  
  PRIVATE
  
  PUBLIC get_dtime, get_beta
  PUBLIC init_BC
  PUBLIC fill_corners, level_pressure
  PUBLIC init_alarm, check_alarm, check_signal
  PUBLIC num_to_string
  PUBLIC iteration_stats
  
  
  INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> Calculates time step.
  !! @todo TEST!!! Noch nicht 100% durchgetested ...
  SUBROUTINE get_dtime
  
  IMPLICIT NONE
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  REAL                   ::  dd1
  REAL                   ::  alpha_re, alpha_im
  REAL                   ::  beta(1:2), beta_global(1:2)
  REAL                   ::  angle, radius
  
  REAL                   ::  kdx_max_conv(1:3)
  REAL                   ::  kdx_max_visc(1:3)
  
  REAL                   ::  array_real(1:2)
  !LOGICAL                ::  array_log (1:5)
  LOGICAL                ::  array_log (1:4)
  
  REAL                   ::  newtime
  REAL                   ::  pi
  
  REAL                   ::  f_max, p_per_side !bbecsek 150510
  REAL                   ::  u_max, h_min, u_max_global, h_min_global
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Die Auswertung bezieht sich generell nur auf den Feldbereich. Bei Ausfluss-Rand-          !
  !                bedingungen ist die CFL-Bedingung naturgemäss irrelevant, Einfluss-Randbedingungen        !
  !                sollten generell implizit in der Zeit diskretisiert werden (Dirichlet-, Neumann- oder     !
  !                Robin-Randbedingungen).                                                                   !
  !              - Bei den Geschwindigkeiten werden der Einfachheit halber nur die Druck-Gitterpunkte be-    !
  !                trachtet, nicht die eigentlich Geschwindigkeitsgitter. Die Abschaetzung sollte aber auf   !
  !                der sicheren Seite liegen.                                                                !
  !              - Explizite Differenzen: k^2_mod_max <= ||L||_inf                                           !
  !              - Implizite Differenzen: k^2_mod_max <= max((kdx_max/dx)^2)  (angenommen, nicht exakt!)     !
  !              - Imaginärteile (z.B. bei schiefen Randstencils) werden der Einfachheit halber vernach-     !
  !                lässigt.                                                                                  !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !--- pi ---
  pi = 2.*ABS(ACOS(0.))
  
  
  dtime_old = dtime
  
  IF ((Int_dtime /= 0 .AND. MOD(timestep,Int_dtime) == 0) .OR. new_dtime) THEN
     
     IF (rank == 0 .AND. write_stout_yes) WRITE(*,'(a)') 'new dtime ...'
     
     new_dtime = .FALSE.
     
     
     kdx_max_conv = pi
     kdx_max_visc = pi
     
     
     ! vel(:,:,:,i) --> worki(:,:,:)
     CALL interpolate_vel(.TRUE.)
     
     beta = 0.
     
     
     !========================================================================================================
     DO k = S31, N31 ! Nur Feldbereich! (Intervallgrenzen von Tangentialgeschwindigkeiten werden verwendet)
        DO j = S21, N21
           DO i = S12, N12
              !-----------------------------------------------------------------------------------------------
                 alpha_re = 0.
                 DO ii = b1L, b1U
                    alpha_re = alpha_re + ABS(cp11(ii,i))
                 END DO
                 DO jj = b2L, b2U
                    alpha_re = alpha_re + ABS(cp22(jj,j))
                 END DO
                 IF (dimens == 3) THEN
                    DO kk = b3L, b3U
                       alpha_re = alpha_re + ABS(cp33(kk,k))
                    END DO
                 END IF
                 alpha_re = alpha_re / Re
              !-----------------------------------------------------------------------------------------------
                 dd1 = 0.
                 DO ii = b1L, b1U
                    dd1 = dd1 + ABS(cp1(ii,i))
                 END DO
                 alpha_im = dd1*ABS(work1(i,j,k))
                 
                 dd1 = 0.
                 DO jj = b2L, b2U
                    dd1 = dd1 + ABS(cp2(jj,j))
                 END DO
                 alpha_im = alpha_im + dd1*ABS(work2(i,j,k))
                 
                 IF (dimens == 3) THEN
                    dd1 = 0.
                    DO kk = b3L, b3U
                       dd1 = dd1 + ABS(cp3(kk,k))
                    END DO
                    alpha_im = alpha_im + dd1*ABS(work3(i,j,k))
                 END IF
              !-----------------------------------------------------------------------------------------------
              IF (ABS(alpha_im) .LT. 10.**(-16)) THEN
                 angle = pi/2.
              ELSE
                 angle = ATAN(alpha_re/alpha_im) ! angle = 0,pi/2
              END IF
              
              radius = SQRT(alpha_re**2 + alpha_im**2)
              
              beta(1) = MAX(beta(1),radius  /stabilitylimit(NINT(2.*REAL(n_stab-1)*angle/pi)))
              beta(2) = MAX(beta(2),alpha_im/stabilitylimit(0                               ))
           END DO
        END DO
     END DO
     !========================================================================================================
      !========================================================================================================
      !--- ADD TIMESTEPPING FOR FSI --------------------------------------------------------------------------
      !========================================================================================================
      !IF (fsi_dtime_yes) THEN
!      h_min = MIN(MINVAL(x1p, MASK = x1p .GT. 0.), MINVAL(x2p, MASK = x2p .GT. 0.), MINVAL(x3p, MASK = x3p .GT. 0.))
!      
!      u_max = 0.
!      DO k = S3p, N3p
!         DO j = S2p, N2p
!            DO i = S1p, N1p
!               u_max = MAX(u_max, SQRT(work1(i,j,k)**2. + work2(i,j,k)**2. + work3(i,j,k)**2.))
!            END DO
!         END DO
!      END DO    
!
!      CALL MPI_REDUCE(h_min,h_min_global,1,MPI_REAL8,MPI_MIN,0,COMM_CART,merror) ! TEST!!! ok?
!      CALL MPI_REDUCE(u_max,u_max_global,1,MPI_REAL8,MPI_MAX,0,COMM_CART,merror) ! TEST!!! ok?

      !END IF
 
     !========================================================================================================
     
     !CALL MPI_ALLREDUCE(beta,beta_global,2,MPI_REAL8,MPI_MAX,COMM_CART,merror)
     CALL MPI_REDUCE(beta,beta_global,2,MPI_REAL8,MPI_MAX,0,COMM_CART,merror) ! TEST!!! ok?
     
     IF (rank == 0) THEN ! TEST!!! ok?
        IF (timeint_mode == 0) THEN
           dtime = beta_global(2) / CFL ! Definition: CFL = 0,1
        ELSE
           dtime = beta_global(1) / CFL ! Definition: CFL = 0,1
        END IF
        

        IF (ABS(dtime) .LT. 1./dtime_max) THEN
           dtime = dtime_max
        ELSE
           dtime = 1./ dtime
        END IF
      
        !========================================================================================================
        !--- ADD TIMESTEPPING FOR FSI --------------------------------------------------------------------------
        !========================================================================================================
        !IF (fsi_dtime_yes) THEN
!           IF (dtime .GE. CFL*h_min_global/u_max_global .AND. timestep .GT. 0) THEN
!              dtime = CFL*h_min_global/u_max_global
!           END IF
        !END IF

        IF (timestep == 0 .AND. dtime .GT. dtime0) dtime = dtime0
     END IF
     
  END IF
  
  !===========================================================================================================
  
  array_real = 0.
  array_log  = .FALSE.
  
  IF (rank == 0) THEN
     ! sicherstellen, dass Zeitschritte nicht zu klein werden:
     newtime = time+dtime + 10.**(-12)
     
     IF (newtime .GE. time_out_vect .AND. dtime_out_vect /= 0.) THEN
        dtime          = time_out_vect - time
        write_out_vect = .TRUE.
        new_dtime      = .TRUE.
     END IF
     
     IF (newtime .GE. time_out_scal .AND. dtime_out_scal /= 0.) THEN
        dtime          = time_out_scal - time
        write_out_scal = .TRUE.
        new_dtime      = .TRUE.
     END IF
     
    ! IF (newtime .GE. time_out_kalm .AND. dtime_out_kalm /= 0.) THEN
    !    dtime          = time_out_kalm - time
    !    write_out_kalm = .TRUE.
    !    new_dtime      = .TRUE.
    ! END IF
 
     IF (newtime .GE. time_end) THEN
        dtime          = time_end - time
        finish_yes     = .TRUE.
     END IF
     
     IF (timestep+1 .GE. n_timesteps) THEN
        finish_yes     = .TRUE.
     END IF
     
     array_real = (/time,dtime/)
     !array_log  = (/finish_yes,new_dtime,write_out_vect,write_out_scal,write_out_kalm/)
     array_log  = (/finish_yes,new_dtime,write_out_vect,write_out_scal/)  
  END IF
  
  CALL MPI_BCAST(array_real,2,MPI_REAL8  ,0,COMM_CART,merror)
  CALL MPI_BCAST(array_log ,5,MPI_LOGICAL,0,COMM_CART,merror)
  
  time  = array_real(1)
  dtime = array_real(2)
  
  finish_yes     = array_log(1)
  new_dtime      = array_log(2)
  write_out_vect = array_log(3)
  write_out_scal = array_log(4)
  !write_out_kalm = array_log(5)
  
  
  IF (timestep == 0) dtime_old = dtime
  
  dtime_average = dtime_average + dtime
  
  
  END SUBROUTINE get_dtime
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! Noch nicht 100% getestet ...
  !> subroutine to calculate the parameter beta that is used for characterizing the condition of the
  !! Helmholtz operator with additional information on the time-step limit.
  SUBROUTINE get_beta
  
  IMPLICIT NONE
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  REAL                   ::  dd1, beta, beta_global
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: -                                                                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
  
     IF (rank == 0) WRITE(*,'(a)') 'determining beta ...'
     
     beta = 0.
     
     !========================================================================================================
     DO k = S31, N31
        DO j = S21, N21
           DO i = S11, N11
              dd1 = 0.
              DO ii = b1L, b1U
                 dd1 = dd1 + ABS(cu11(ii,i))
              END DO
              DO jj = b2L, b2U
                 dd1 = dd1 + ABS(cp22(jj,j))
              END DO
              IF (dimens == 3) THEN
                 DO kk = b3L, b3U
                    dd1 = dd1 + ABS(cp33(kk,k))
                 END DO
              END IF
              beta = MAX(beta,dd1)
           END DO
        END DO
     END DO
     
     DO k = S32, N32
        DO j = S22, N22
           DO i = S12, N12
              dd1 = 0.
              DO ii = b1L, b1U
                 dd1 = dd1 + ABS(cp11(ii,i))
              END DO
              DO jj = b2L, b2U
                 dd1 = dd1 + ABS(cv22(jj,j))
              END DO
              IF (dimens == 3) THEN
                 DO kk = b3L, b3U
                    dd1 = dd1 + ABS(cp33(kk,k))
                 END DO
              END IF
              beta = MAX(beta,dd1)
           END DO
        END DO
     END DO
     
     IF (dimens == 3) THEN
        DO k = S33, N33
           DO j = S23, N23
              DO i = S13, N13
                 dd1 = 0.
                 DO ii = b1L, b1U
                    dd1 = dd1 + ABS(cp11(ii,i))
                 END DO
                 DO jj = b2L, b2U
                    dd1 = dd1 + ABS(cp22(jj,j))
                 END DO
                 DO kk = b3L, b3U
                    dd1 = dd1 + ABS(cw33(kk,k))
                 END DO
                 beta = MAX(beta,dd1)
              END DO
           END DO
        END DO
     END IF
     
     !========================================================================================================
     
     beta = beta / Re
     
     !CALL MPI_ALLREDUCE(beta,beta_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror)
     CALL MPI_REDUCE(beta,beta_global,1,MPI_REAL8,MPI_MAX,0,COMM_CART,merror) ! TEST!!!
     
     IF (rank == 0) THEN
        OPEN(10,FILE='test_beta.txt',STATUS='UNKNOWN')
        WRITE(10,'(a,E25.17)') '||L||_infty  = ', beta_global
        CLOSE(10)
     END IF
     
  
  
  END SUBROUTINE get_beta
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE init_BC ! TEST!!! Name nicht ganz passend ...
  
  IMPLICIT NONE
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  INTEGER                ::  m
  
  
  !bc11 = 0. ! TEST!!! nach usr_initcond.f90 bzw. usr_boundcond.f90 verschoben ...
  !bc12 = 0.
  !bc13 = 0.
  !
  !bc21 = 0.
  !bc22 = 0.
  !bc23 = 0.
  !
  !bc31 = 0.
  !bc32 = 0.
  !bc33 = 0.
  
  
  !===========================================================================================================
  !=== Ausfluss-RB initialisieren ============================================================================
  !===========================================================================================================
  IF (outlet(1,1,2) .AND. (BC_1L == 1 .OR. BC_1L == 3)) bc21(S22B:N22B,S32B:N32B,1) = vel(1, S22B:N22B,S32B:N32B,2)
  IF (outlet(1,2,2) .AND. (BC_1U == 1 .OR. BC_1U == 3)) bc21(S22B:N22B,S32B:N32B,2) = vel(N1,S22B:N22B,S32B:N32B,2)
  
  IF (outlet(1,1,3) .AND. (BC_1L == 1 .OR. BC_1L == 3)) bc31(S23B:N23B,S33B:N33B,1) = vel(1 ,S23B:N23B,S33B:N33B,3)
  IF (outlet(1,2,3) .AND. (BC_1U == 1 .OR. BC_1U == 3)) bc31(S23B:N23B,S33B:N33B,2) = vel(N1,S23B:N23B,S33B:N33B,3)
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(2,1,1) .AND. (BC_2L == 1 .OR. BC_2L == 3)) bc12(S11B:N11B,S31B:N31B,1) = vel(S11B:N11B,1 ,S31B:N31B,1)
  IF (outlet(2,2,1) .AND. (BC_2U == 1 .OR. BC_2U == 3)) bc12(S11B:N11B,S31B:N31B,2) = vel(S11B:N11B,N2,S31B:N31B,1)
  
  IF (outlet(2,1,3) .AND. (BC_2L == 1 .OR. BC_2L == 3)) bc32(S13B:N13B,S33B:N33B,1) = vel(S13B:N13B,1 ,S33B:N33B,3)
  IF (outlet(2,2,3) .AND. (BC_2U == 1 .OR. BC_2U == 3)) bc32(S13B:N13B,S33B:N33B,2) = vel(S13B:N13B,N2,S33B:N33B,3)
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(3,1,1) .AND. (BC_3L == 1 .OR. BC_3L == 3)) bc13(S11B:N11B,S21B:N21B,1) = vel(S11B:N11B,S21B:N21B,1 ,1)
  IF (outlet(3,2,1) .AND. (BC_3U == 1 .OR. BC_3U == 3)) bc13(S11B:N11B,S21B:N21B,2) = vel(S11B:N11B,S21B:N21B,N3,1)
  
  IF (outlet(3,1,2) .AND. (BC_3L == 1 .OR. BC_3L == 3)) bc23(S12B:N12B,S22B:N22B,1) = vel(S12B:N12B,S22B:N22B,1 ,2)
  IF (outlet(3,2,2) .AND. (BC_3U == 1 .OR. BC_3U == 3)) bc23(S12B:N12B,S22B:N22B,2) = vel(S12B:N12B,S22B:N22B,N3,2)
  !===========================================================================================================
  IF (outlet(1,1,1) .AND. (BC_1L == 1 .OR. BC_1L == 2)) THEN
     i = 1
     DO k = S31B, N31B
        DO j = S21B, N21B
           bc11(j,k,1) = cIup(d1L,i)*vel(i+d1L,j,k,1)
!pgi$ unroll = n:8
           DO ii = d1L+1, d1U
              bc11(j,k,1) = bc11(j,k,1) + cIup(ii,i)*vel(i+ii,j,k,1)
           END DO
        END DO
     END DO
  END IF
  IF (outlet(1,2,1) .AND. (BC_1U == 1 .OR. BC_1U == 2)) THEN
     i = N1
     DO k = S31B, N31B
        DO j = S21B, N21B
           bc11(j,k,2) = cIup(d1L,i)*vel(i+d1L,j,k,1)
!pgi$ unroll = n:8
           DO ii = d1L+1, d1U
              bc11(j,k,2) = bc11(j,k,2) + cIup(ii,i)*vel(i+ii,j,k,1)
           END DO
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(2,1,2) .AND. (BC_2L == 1 .OR. BC_2L == 2)) THEN
     j = 1
     DO k = S32B, N32B
        DO i = S12B, N12B
           bc22(i,k,1) = cIvp(d2L,j)*vel(i,j+d2L,k,2)
!pgi$ unroll = n:8
           DO jj = d2L+1, d2U
              bc22(i,k,1) = bc22(i,k,1) + cIvp(jj,j)*vel(i,j+jj,k,2)
           END DO
        END DO
     END DO
  END IF
  IF (outlet(2,2,2) .AND. (BC_2U == 1 .OR. BC_2U == 2)) THEN
     j = N2
     DO k = S32B, N32B
        DO i = S12B, N12B
           bc22(i,k,2) = cIvp(d2L,j)*vel(i,j+d2L,k,2)
!pgi$ unroll = n:8
           DO jj = d2L+1, d2U
              bc22(i,k,2) = bc22(i,k,2) + cIvp(jj,j)*vel(i,j+jj,k,2)
           END DO
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(3,1,3) .AND. (BC_3L == 1 .OR. BC_3L == 2)) THEN
     k = 1
     DO j = S23B, N23B
        DO i = S13B, N13B
           bc33(i,j,1) = cIwp(d3L,k)*vel(i,j,k+d3L,3)
!pgi$ unroll = n:8
           DO kk = d3L+1, d3U
              bc33(i,j,1) = bc33(i,j,1) + cIwp(kk,k)*vel(i,j,k+kk,3)
           END DO
        END DO
     END DO
  END IF
  IF (outlet(3,2,3) .AND. (BC_3U == 1 .OR. BC_3U == 2)) THEN
     k = N3
     DO j = S23B, N23B
        DO i = S13B, N13B
           bc33(i,j,2) = cIwp(d3L,k)*vel(i,j,k+d3L,3)
!pgi$ unroll = n:8
           DO kk = d3L+1, d3U
              bc33(i,j,2) = bc33(i,j,2) + cIwp(kk,k)*vel(i,j,k+kk,3)
           END DO
        END DO
     END DO
  END IF
  !===========================================================================================================
  
  
  
  END SUBROUTINE init_BC
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE level_pressure
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k
  REAL                   ::  pre0, pre0_global
  
  
  IF ((Int_lev_pre /= 0 .AND. MOD(timestep,Int_lev_pre) == 0) .OR. write_out_vect) THEN
     
     pre0 = 0.
     
     DO k = S3p, N3p
        DO j = S2p, N2p
!pgi$ unroll = n:8
           DO i = S1p, N1p
              pre0 = pre0 + pre(i,j,k)
           END DO
        END DO
     END DO
     
     CALL MPI_ALLREDUCE(pre0,pre0_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     
     pre0 = pre0_global/REAL(dim1)/REAL(dim2)/REAL(dim3) ! TEST!!! wegen i4 gefaehrlich!
     
     pre(S1p:N1p,S2p:N2p,S3p:N3p) = pre(S1p:N1p,S2p:N2p,S3p:N3p) - pre0
     
  END IF
  
  
  END SUBROUTINE level_pressure
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE fill_corners(phi)
  
  IMPLICIT NONE
  
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  
  ! Kanten/Ecken sind durch die Randbedingungen immer Divergenz-frei, so dass der Druck hier unbestimmt ist.
  ! Hier Mittelung der Nachbarpunkte
  
  IF (BC_1L .GT. 0 .AND. BC_2L .GT. 0) phi(1 ,1 ,1:N3) = (phi(2   ,1 ,1:N3) + phi(1 ,2   ,1:N3) + phi(2   ,2   ,1:N3))/3.
  IF (BC_1L .GT. 0 .AND. BC_2U .GT. 0) phi(1 ,N2,1:N3) = (phi(2   ,N2,1:N3) + phi(1 ,N2-1,1:N3) + phi(2   ,N2-1,1:N3))/3.
  IF (BC_1U .GT. 0 .AND. BC_2L .GT. 0) phi(N1,1 ,1:N3) = (phi(N1-1,1 ,1:N3) + phi(N1,2   ,1:N3) + phi(N1-1,2   ,1:N3))/3.
  IF (BC_1U .GT. 0 .AND. BC_2U .GT. 0) phi(N1,N2,1:N3) = (phi(N1-1,N2,1:N3) + phi(N1,N2-1,1:N3) + phi(N1-1,N2-1,1:N3))/3.

  IF (BC_1L .GT. 0 .AND. BC_3L .GT. 0) phi(1 ,1:N2,1 ) = (phi(2   ,1:N2,1 ) + phi(1 ,1:N2,2   ) + phi(2   ,1:N2,2   ))/3.
  IF (BC_1L .GT. 0 .AND. BC_3U .GT. 0) phi(1 ,1:N2,N3) = (phi(2   ,1:N2,N3) + phi(1 ,1:N2,N3-1) + phi(2   ,1:N2,N3-1))/3.
  IF (BC_1U .GT. 0 .AND. BC_3L .GT. 0) phi(N1,1:N2,1 ) = (phi(N1-1,1:N2,1 ) + phi(N1,1:N2,2   ) + phi(N1-1,1:N2,2   ))/3.
  IF (BC_1U .GT. 0 .AND. BC_3U .GT. 0) phi(N1,1:N2,N3) = (phi(N1-1,1:N2,N3) + phi(N1,1:N2,N3-1) + phi(N1-1,1:N2,N3-1))/3.

  IF (BC_2L .GT. 0 .AND. BC_3L .GT. 0) phi(1:N1,1 ,1 ) = (phi(1:N1,2   ,1 ) + phi(1:N1,1 ,2   ) + phi(1:N1,2   ,2   ))/3.
  IF (BC_2L .GT. 0 .AND. BC_3U .GT. 0) phi(1:N1,1 ,N3) = (phi(1:N1,2   ,N3) + phi(1:N1,1 ,N3-1) + phi(1:N1,2   ,N3-1))/3.
  IF (BC_2U .GT. 0 .AND. BC_3L .GT. 0) phi(1:N1,N2,1 ) = (phi(1:N1,N2-1,1 ) + phi(1:N1,N2,2   ) + phi(1:N1,N2-1,2   ))/3.
  IF (BC_2U .GT. 0 .AND. BC_3U .GT. 0) phi(1:N1,N2,N3) = (phi(1:N1,N2-1,N3) + phi(1:N1,N2,N3-1) + phi(1:N1,N2-1,N3-1))/3.
  
  
  END SUBROUTINE fill_corners
  
  
  
  
  
  
  
  
  
  
  !> Checks for file 'queue.txt' and initializes an alarm with the custom print routine check_alarm thath also sets
  !! the variable 'finish_yes=.TRUE.'
  SUBROUTINE init_alarm
  
  IMPLICIT NONE
  
  INTEGER                ::  wtime, margin ! TEST!!! margin nach jobscript verlagert!
  INTEGER                ::  ios
  
  
  !*************
  margin = 0 ! [seconds] ! TEST!!! margin nach jobscript verlagert!
  !*************
  
  
  IF (rank == 0) THEN
     OPEN(10,FILE='queue.txt',ACTION='read',STATUS='old',IOSTAT=ios)
     
     IF (ios == 0) THEN
        READ (10,*) wtime
        CLOSE(10)
        
        IF (wtime .GT. margin) wtime = wtime - margin ! TEST!!! margin nach jobscript verlagert!
        
        WRITE(*,*)
        WRITE(*,'(a,i6,a)') '... found queue.txt, exiting after', wtime, ' seconds ...'
        WRITE(*,*)
        
        CALL ALARM(wtime,check_alarm)
     ELSE
        WRITE(*,*)
        WRITE(*,*) '... found no queue.txt ...'
        WRITE(*,*)
     END IF
  END IF
  
  
  END SUBROUTINE init_alarm
  
  
  
  
  
  
  
  
  
  
  !> Routine for displaying that the alarm sign was received and for setting 'finish_yes=.TRUE.'
  SUBROUTINE check_alarm
  
  IMPLICIT NONE
  
  
  IF (rank == 0) THEN
     
     WRITE(*,*) '+++++++++++++++++++++++++++++++++++'
     WRITE(*,*) '+ ... received the alarm sign ... +'
     WRITE(*,*) '+++++++++++++++++++++++++++++++++++'
     
     finish_yes = .TRUE.
     
  END IF
  
  
  END SUBROUTINE check_alarm
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE check_signal
  
  IMPLICIT NONE
  
  INTEGER                ::  signal
  INTEGER                ::  timestep_end
  REAL                   ::  time_end
  
  
  IF (rank == 0) THEN
     
     OPEN (10, FILE='send_signal.txt', STATUS='UNKNOWN')
     READ (10,*) signal
     READ (10,*) timestep_end
     READ (10,*) time_end
     CLOSE(10)
     
     IF (signal == 1) finish_yes = .TRUE.
     IF (signal == 2 .AND. (timestep == timestep_end .OR. time .GE. time_end)) finish_yes = .TRUE.
     
  END IF
  
  CALL MPI_BCAST(finish_yes,1,MPI_LOGICAL,0,COMM_CART,merror)
  
  
  END SUBROUTINE check_signal
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE num_to_string(n_digits,num,num_char)
  
  IMPLICIT NONE
  
  INTEGER     , INTENT(in ) ::  n_digits
  INTEGER     , INTENT(in ) ::  num
  CHARACTER(LEN=n_digits), INTENT(out) ::  num_char
  
  INTEGER                   ::  tenthousands, thousands, hundreds, tens, ones
  
  
  !===========================================================================================================
  IF (num .LT. 0) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot convert negative integers to a string.'
     CALL MPI_FINALIZE(merror)
     STOP
  !===========================================================================================================
  ELSE IF (n_digits .GE. 9) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot convert integers > 99999999 to a string.'
     CALL MPI_FINALIZE(merror)
     STOP
  !===========================================================================================================
  ELSE IF (n_digits == 1) THEN
     WRITE(num_char,'(i1.1)') num
  !===========================================================================================================
  ELSE IF (n_digits == 2) THEN
     WRITE(num_char,'(i2.2)') num
  !===========================================================================================================
  ELSE IF (n_digits == 3) THEN
     WRITE(num_char,'(i3.3)') num
  !===========================================================================================================
  ELSE IF (n_digits == 4) THEN
     WRITE(num_char,'(i4.4)') num
  !===========================================================================================================
  ELSE IF (n_digits == 5) THEN
     !tenthousands =  num / 10000
     !thousands    = (num - 10000 * tenthousands) / 1000
     !hundreds     = (num - 10000 * tenthousands - 1000 * thousands) / 100
     !tens         = (num - 10000 * tenthousands - 1000 * thousands - 100 * hundreds) / 10
     !ones         =  num - 10000 * tenthousands - 1000 * thousands - 100 * hundreds - tens*10
     !
     !num_char(1:1) = CHAR(ICHAR('0')+tenthousands)
     !num_char(2:2) = CHAR(ICHAR('0')+thousands)
     !num_char(3:3) = CHAR(ICHAR('0')+hundreds)
     !num_char(4:4) = CHAR(ICHAR('0')+tens)
     !num_char(5:5) = CHAR(ICHAR('0')+ones)
     WRITE(num_char,'(i5.5)') num
  !===========================================================================================================
  ELSE IF (n_digits == 6) THEN
     WRITE(num_char,'(i6.6)') num
  !===========================================================================================================
  ELSE IF (n_digits == 7) THEN
     WRITE(num_char,'(i7.7)') num
  !===========================================================================================================
  ELSE IF (n_digits == 8) THEN
     WRITE(num_char,'(i8.8)') num
  END IF
  !===========================================================================================================
  
  END SUBROUTINE num_to_string
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE iteration_stats
  
  ! revised: 24.10.07
  
  IMPLICIT NONE
  
  REAL                ::  ratioO_tot
  REAL                ::  ratioH_tot
  REAL                ::  ratioP_tot
  
  INTEGER             ::  countO_tot
  INTEGER             ::  countP_tot
  INTEGER             ::  countH_tot
  
  REAL                ::  ratioP_all
  REAL                ::  ratioH_all
  
  INTEGER             ::  countP_all
  INTEGER             ::  countH_all
  
  REAL                ::  countO_av
  REAL                ::  countP_av(1:2)
  REAL                ::  countH_av(1:3)
  REAL                ::  countP_all_av
  REAL                ::  countH_all_av
  
  INTEGER             ::  tsteps
  
  
  tsteps = timestep - timestep_old
  
  
  IF (rank == 0 .AND. tsteps .GE. 1) THEN
     OPEN(10, FILE='test_iteration_stats_restart'//restart_char//'.txt', STATUS='UNKNOWN')
     
     dtime_average = dtime_average / REAL(tsteps)
     
     countO_tot = 0
     countP_tot = 0
     countH_tot = 0
     
     ratioO_tot = 0.
     ratioP_tot = 0.
     ratioH_tot = 0.
     
     WRITE(10,'(a       )') ''
     WRITE(10,'(a, G13.5)') '                        Re =', Re
     WRITE(10,'(a,3E13.5)') '                Dimensions =', L1 , L2 , L3
     WRITE(10,'(a,3i5   )') '                Resolution =', M1 , M2 , M3
     WRITE(10,'(a,3i5   )') '                     Procs =', NB1, NB2, NB3
     WRITE(10,'(a, E13.5)') '                   <dtime> =', dtime_average
     WRITE(10,'(a, i7   )') '               n_timesteps =', tsteps
     WRITE(10,'(a, E13.5)') '        elapsed time [sec] =', REAL(elatime)/1000.
     WRITE(10,'(a, E13.5)') 'el. time/n_timesteps [sec] =', REAL(elatime)/1000./REAL(tsteps)
     WRITE(10,'(a,2E13.5)') '               div(vel(0)) =', max_div_init(1:2)
     WRITE(10,'(a       )') ''
     WRITE(10,'(a       )') 'number of iterations (number/timestep, <conv. ratio>):'
     WRITE(10,'(a       )') ''
     
     DO substep = 1, RK_steps
        
        countP_all = countP(substep,1) + countP(substep,2)
        countH_all = countH(substep,1) + countH(substep,2) + countH(substep,3)
        
        ratioP_all = ratioP(substep,1) + ratioP(substep,2)
        ratioH_all = ratioH(substep,1) + ratioH(substep,2) + ratioH(substep,3)
        
        !-----------------------------------------------------------------------------------------------------
        
        countO_tot = countO_tot + countO(substep)
        countP_tot = countP_tot + countP_all
        countH_tot = countH_tot + countH_all
        
        ratioO_tot = ratioO_tot + ratioO(substep)
        ratioP_tot = ratioP_tot + ratioP_all
        ratioH_tot = ratioH_tot + ratioH_all
        
        !-----------------------------------------------------------------------------------------------------
        
        countO_av      = REAL(countO(substep))     / REAL(tsteps)
        countP_av(1:2) = REAL(countP(substep,1:2)) / REAL(tsteps)
        countH_av(1:3) = REAL(countH(substep,1:3)) / REAL(tsteps)
        countP_all_av  = REAL(countP_all)          / REAL(tsteps)
        countH_all_av  = REAL(countH_all)          / REAL(tsteps)
        
        
        IF (countO(substep)   /= 0) ratioO(substep)   = ratioO(substep)   / REAL(countO(substep)  )
        
        IF (countP(substep,1) /= 0) ratioP(substep,1) = ratioP(substep,1) / REAL(countP(substep,1))
        IF (countP(substep,2) /= 0) ratioP(substep,2) = ratioP(substep,2) / REAL(countP(substep,2))
        
        IF (countH(substep,1) /= 0) ratioH(substep,1) = ratioH(substep,1) / REAL(countH(substep,1))
        IF (countH(substep,2) /= 0) ratioH(substep,2) = ratioH(substep,2) / REAL(countH(substep,2))
        IF (countH(substep,3) /= 0) ratioH(substep,3) = ratioH(substep,3) / REAL(countH(substep,3))
        
        IF (countP_all        /= 0) ratioP_all        = ratioP_all        / REAL(countP_all       )
        IF (countH_all        /= 0) ratioH_all        = ratioH_all        / REAL(countH_all       )
        
        !-----------------------------------------------------------------------------------------------------
        
        WRITE(10,'(a,i2,a  )') '   substep', substep, ':'
        WRITE(10,'(a       )') '      '
        WRITE(10,'(a,2G12.4)') '      outer Iterations:', countO_av    , ratioO    (substep)
        WRITE(10,'(a       )') '         '
        WRITE(10,'(a,2G12.4)') '         Poisson   (1):', countP_av(1) , ratioP    (substep,1)
        WRITE(10,'(a,2G12.4)') '         Poisson   (2):', countP_av(2) , ratioP    (substep,2)
        WRITE(10,'(a,2G12.4)') '         all          :', countP_all_av, ratioP_all
        WRITE(10,'(a       )') '         '
        WRITE(10,'(a,2G12.4)') '         Helmholtz (1):', countH_av(1) , ratioH    (substep,1)
        WRITE(10,'(a,2G12.4)') '         Helmholtz (2):', countH_av(2) , ratioH    (substep,2)
        WRITE(10,'(a,2G12.4)') '         Helmholtz (3):', countH_av(3) , ratioH    (substep,3)
        WRITE(10,'(a,2G12.4)') '         all          :', countH_all_av, ratioH_all
        WRITE(10,'(a       )') '         '
        
     END DO
     
     IF (countO_tot /= 0) ratioO_tot = ratioO_tot / REAL(countO_tot)
     IF (countP_tot /= 0) ratioP_tot = ratioP_tot / REAL(countP_tot)
     IF (countH_tot /= 0) ratioH_tot = ratioH_tot / REAL(countH_tot)
     
     WRITE(10,'(a       )') '   total:'
     WRITE(10,'(a       )') '      '
     WRITE(10,'(a,2G12.4)') '      outer Iterations:', REAL(countO_tot)/REAL(tsteps), ratioO_tot
     WRITE(10,'(a,2G12.4)') '         Poisson      :', REAL(countP_tot)/REAL(tsteps), ratioP_tot
     WRITE(10,'(a,2G12.4)') '         Helmholtz    :', REAL(countH_tot)/REAL(tsteps), ratioH_tot
     WRITE(10,'(a       )') '         '
     
     CLOSE(10)
     
  END IF
  
  END SUBROUTINE iteration_stats
  
  
  
  
END MODULE mod_lib
