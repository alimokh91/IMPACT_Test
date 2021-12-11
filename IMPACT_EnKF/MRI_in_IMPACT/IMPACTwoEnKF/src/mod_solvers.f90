!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch                    *
!* October 2014                                                                                              *
!*************************************************************************************************************

!> module containing solver routines and preconditioners for solving the NSE. It uses modules mod_dim,
!! mod_vars, mod_exchange, mod_diff, mod_laplace, mod_helmholtz and mod_inout.
!! @note: what is the difference between BiCGstab and BiCGstab2? (one has an extra option for problem_type, 
!!        anything else?)
MODULE mod_solvers
  
  
  USE mod_dims
  USE mod_vars
  USE mod_exchange
  USE mod_diff
  USE mod_laplace
  USE mod_helmholtz
  USE mod_inout
  !USE mpi

  PRIVATE
  
  PUBLIC outer_iteration, explicit, twostep
  PUBLIC force_massflow, apply_nullspace, get_nullspace, solve_nullspace
  PUBLIC solve_Helmholtz
  PUBLIC multigridV, multigridF, restrict, interpolate
  PUBLIC restrict_Helmholtz, interpolate_Helmholtz
  PUBLIC BiCGstab, Richardson
  PUBLIC get_norms, product_scalar, multadd1, multadd2
  PUBLIC status_iteration
  PUBLIC handle_corner_rhs
  
  PUBLIC BiCGstab2, get_norms2, product_scalar2, apply_nullspace2 ! TEST!!!
  PUBLIC relax_restrict, interpolate_relax, plain_restrict, interpolate_mg, relax_bottom
  
  
  INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  SUBROUTINE outer_iteration
  
  IMPLICIT NONE
  
  REAL                   ::  norm, norm_prev
  REAL                   ::  epsP
  
  INTEGER                ::  i, j, k
  INTEGER                ::  m, counter
  LOGICAL                ::  exit_yes
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Beim zweiten Poisson-Problem im Vorkonditionierer hat es sich bislang bewährt, die Lösung !
  !                aus dem ersten Poisson-Problem als Startfeld zu benutzen.                                 !
  !              - Richardson-Iteration scheint langsamer als BiCGstab zu sein (~30%).                       !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== RHS fuer Druck-Gleichungssystem =======================================================================
  !===========================================================================================================
  IF (rank == 0 .AND. log_iteration_yes) WRITE(10,'(a)') 'solving Helmholtz equations for RHS ...'
  
  DO m = 1, dimens
     
     !--- Helmholtz-Gleichung lösen --------------------------------------------------------------------------
     IF (init_pre(substep)) THEN
        CALL solve_Helmholtz(m,epsU,n_it_Helmh_vel,init_vel(substep),rhs(b1L,b2L,b3L,m),vel(b1L,b2L,b3L,m),.TRUE.,.TRUE.)
     ELSE
        CALL gradient(m,pre,gpre)
        
        IF (m == 1) gpre(S11B:N11B,S21B:N21B,S31B:N31B) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) - gpre(S11B:N11B,S21B:N21B,S31B:N31B)
        IF (m == 2) gpre(S12B:N12B,S22B:N22B,S32B:N32B) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) - gpre(S12B:N12B,S22B:N22B,S32B:N32B)
        IF (m == 3) gpre(S13B:N13B,S23B:N23B,S33B:N33B) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) - gpre(S13B:N13B,S23B:N23B,S33B:N33B)
        
        CALL solve_Helmholtz(m,epsU,n_it_Helmh_vel,init_vel(substep),gpre,vel(b1L,b2L,b3L,m),.TRUE.,.TRUE.)
     END IF
     
  END DO
  
  !--- Divergenz bilden --------------------------------------------------------------------------------------
  CALL divergence2(vel,res)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Norm des Residuums ====================================================================================
  !===========================================================================================================
  IF (corner_yes) CALL handle_corner_rhs(1,res)
  CALL get_norms(S1p,S2p,S3p,N1p,N2p,N3p,res,2,.TRUE.,.FALSE.,normInf=norm)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Druckiteration ========================================================================================
  !===========================================================================================================
  counter = 0
  
  LOOP: DO
     
     !========================================================================================================
     !=== Norm an Abbruchkriterium testen ====================================================================
     !========================================================================================================
     IF (rank == 0 .AND. log_iteration_yes) THEN
        WRITE(10,'(a)')
        WRITE(10,'(a)') '-------------------------------------------------------------------------------------'
        WRITE(10,'(a)') 'current residual:'
     END IF
     
     CALL status_iteration(epsU,norm,counter,n_it_outer,exit_yes,.FALSE.,.FALSE.)
     
     ! ACHTUNG!!! Das findet die xt3 nicht so toll:
     IF (rank == 0 .AND. log_iteration_yes) WRITE(10,'(a)') '-------------------------------------------------------------------------------------'
     IF (rank == 0 .AND. log_iteration_yes) CALL flush(10)
     
     IF (exit_yes .OR. norm == 0.) EXIT LOOP
     !========================================================================================================
     
     
     !========================================================================================================
     !=== nächste Iteration ==================================================================================
     !========================================================================================================
     counter = counter + 1
     
     IF (rank == 0 .AND. log_iteration_yes) WRITE(10,'(a)')
     IF (rank == 0 .AND. log_iteration_yes) WRITE(10,'(a)') 'solving Poisson equations ...'
     
     epsP = precRatio(substep,counter)*norm
     
     IF (epsP .LE. epsU) epsP = epsU
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Vorkonditionierer ==================================================================================
     !========================================================================================================
     IF (precond_outer == 1 .OR. precond_outer == 2) THEN
        ! RHS korrigieren (sollte eigentlich nicht mehr notwendig sein): ! TEST!!!
        IF (nullspace_yes) CALL apply_nullspace(res)
        
        number_poisson = 1
        CALL BiCGstab  (precOffset(substep,counter)*epsP,n_it_Poisson,.TRUE.,S1p,S2p,S3p,N1p,N2p,N3p,                &
                   &                       res,work3,2,.TRUE.,.TRUE.,precond_Poisson)
        !CALL Richardson(precOffset(substep,counter)*epsP,n_it_Poisson,.TRUE.,S1p,S2p,S3p,N1p,N2p,N3p,               &
        !           &                       res,work3,2,.TRUE.,.TRUE.,precond_Poisson)
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (precond_outer == 2) THEN
        
        DO m = 1, dimens
           direction = m
           CALL gradient(m,work3,gpre)
           CALL bc_extrapolation(m,gpre)
           CALL product_Helmholtz(gpre,work1)
           CALL bc_extrapolation(m,work1)
           CALL divergence(m,work1,work2)
        END DO
        
        ! Kanten/Ecken sind durch die Randbedingungen immer Divergenz-frei:
        IF (corner_yes   ) CALL handle_corner_rhs(1,work2)
        
        ! RHS korrigieren:
        IF (nullspace_yes) CALL apply_nullspace(work2)
        
        number_poisson = 2
        CALL BiCGstab  (precOffset(substep,counter)*epsP,n_it_Poisson,.FALSE.,S1p,S2p,S3p,N1p,N2p,N3p,               &
                   &                       work2,work3,2,.TRUE.,.TRUE.,precond_Poisson)
        !CALL Richardson(precOffset(substep,counter)*epsP,n_it_Poisson,.FALSE.,S1p,S2p,S3p,N1p,N2p,N3p,              &
        !           &                       work2,work3,2,.TRUE.,.TRUE.,precond_Poisson)
     END IF
     !========================================================================================================
     
     
     !========================================================================================================
     !=== aktuelles Residuum =================================================================================
     !========================================================================================================
     IF (init_pre(substep) .AND. counter == 1) THEN
        pre(S1p:N1p,S2p:N2p,S3p:N3p) = work3(S1p:N1p,S2p:N2p,S3p:N3p)
     ELSE
        pre(S1p:N1p,S2p:N2p,S3p:N3p) = pre(S1p:N1p,S2p:N2p,S3p:N3p) + work3(S1p:N1p,S2p:N2p,S3p:N3p)
        !pre(S1p:N1p,S2p:N2p,S3p:N3p) = 1.25*pre(S1p:N1p,S2p:N2p,S3p:N3p) + work3(S1p:N1p,S2p:N2p,S3p:N3p) ! Relaxation ...
     END IF
     
     IF (rank == 0 .AND. log_iteration_yes) WRITE(10,'(a)')
     IF (rank == 0 .AND. log_iteration_yes) WRITE(10,'(a)') 'solving Helmholtz equations for residual ...'
     
     DO m = 1, dimens
        
        CALL gradient(m,work3,gpre)
        CALL solve_Helmholtz(m,epsU,n_it_Helmh_vel,.TRUE.,gpre,work1,.TRUE.,.TRUE.)
        
        IF (m == 1) vel(S11B:N11B,S21B:N21B,S31B:N31B,1) = vel(S11B:N11B,S21B:N21B,S31B:N31B,1) - work1(S11B:N11B,S21B:N21B,S31B:N31B)
        IF (m == 2) vel(S12B:N12B,S22B:N22B,S32B:N32B,2) = vel(S12B:N12B,S22B:N22B,S32B:N32B,2) - work1(S12B:N12B,S22B:N22B,S32B:N32B)
        IF (m == 3) vel(S13B:N13B,S23B:N23B,S33B:N33B,3) = vel(S13B:N13B,S23B:N23B,S33B:N33B,3) - work1(S13B:N13B,S23B:N23B,S33B:N33B)
        
     END DO
     
     !--- Divergenz bilden -----------------------------------------------------------------------------------
     CALL divergence2(vel,res)
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Norm des Residuums =================================================================================
     !========================================================================================================
     norm_prev = norm
     IF (corner_yes) CALL handle_corner_rhs(1,res)
     CALL get_norms(S1p,S2p,S3p,N1p,N2p,N3p,res,2,.TRUE.,.FALSE.,normInf=norm)
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Abbruchkriterium für innere Iteration ==============================================================
     !========================================================================================================
     precRatio(substep,counter) = norm/norm_prev
     IF (precRatio(substep,counter) .GT. 1.) precRatio(substep,counter) = 1.
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Iterationsstatistik ================================================================================
     !========================================================================================================
     ratioO(substep) = ratioO(substep) + LOG10(norm/norm_prev)
     countO(substep) = countO(substep) + 1
     !========================================================================================================
     
  END DO LOOP
  
  
  !===========================================================================================================
  !=== Massenfluss korrigieren ===============================================================================
  !===========================================================================================================
  IF (forcing_mode == 2) CALL force_massflow(.TRUE.)
  
  
  END SUBROUTINE outer_iteration
  
  
  

  
  
  
  
  
  SUBROUTINE explicit
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k, m
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Beim zweiten Poisson-Problem im Vorkonditionierer hat es sich bislang bewährt, die Lösung !
  !                aus dem ersten Poisson-Problem als Startfeld zu benutzen.                                 !
  !              - Richardson-Iteration scheint langsamer als BiCGstab zu sein (~30%).                       !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== RHS fuer Druck-Gleichungssystem =======================================================================
  !===========================================================================================================
  DO m = 1, dimens
     
     ! TEST!!! weglassen und stattdessen rhs verwenden? Testroutinen lassen sich dann nicht mehr verwenden ...
     IF (m == 1) vel(S11B:N11B,S21B:N21B,S31B:N31B,1) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1)
     IF (m == 2) vel(S12B:N12B,S22B:N22B,S32B:N32B,2) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2)
     IF (m == 3) vel(S13B:N13B,S23B:N23B,S33B:N33B,3) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3)
     
     !--- Randbedingungen extrapolieren ----------------------------------------------------------------------
     CALL bc_extrapolation(m,vel(b1L,b2L,b3L,m))
     
  END DO
  
  !--- Divergenz bilden --------------------------------------------------------------------------------------
  CALL divergence2(vel,res)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Poisson-Lösung ========================================================================================
  !===========================================================================================================
  ! Kanten/Ecken sind durch die Randbedingungen immer Divergenz-frei:
  IF (corner_yes) CALL handle_corner_rhs(1,res)
  
  IF (rank == 0 .AND. log_iteration_yes) WRITE(10,'(a)') 'Poisson problem:'
  
  number_poisson = 1
  CALL BiCGstab  (epsU,n_it_Poisson,init_pre(substep),S1p,S2p,S3p,N1p,N2p,N3p,res,pre,2,.TRUE.,.TRUE.,precond_Poisson)
  !CALL Richardson(epsU,n_it_Poisson,init_pre(substep),S1p,S2p,S3p,N1p,N2p,N3p,res,pre,2,.TRUE.,.TRUE.,precond_Poisson)
  
  IF (rank == 0 .AND. log_iteration_yes) CALL flush(10)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Korrektur der Lösung ==================================================================================
  !===========================================================================================================
  DO m = 1, dimens
     
     CALL gradient(m,pre,gpre)
     
     !--- Randbedingungen extrapolieren ----------------------------------------------------------------------
     CALL bc_extrapolation(m,gpre)
     
     
     IF (m == 1) THEN
        DO k = S31B, N31B
           DO j = S21B, N21B
!pgi$ unroll = n:8
              DO i = S11B, N11B
                 vel(i,j,k,m) = vel(i,j,k,m) - gpre(i,j,k)
              END DO
           END DO
        END DO
     END IF
     IF (m == 2) THEN
        DO k = S32B, N32B
           DO j = S22B, N22B
!pgi$ unroll = n:8
              DO i = S12B, N12B
                 vel(i,j,k,m) = vel(i,j,k,m) - gpre(i,j,k)
              END DO
           END DO
        END DO
     END IF
     IF (m == 3) THEN
        DO k = S33B, N33B
           DO j = S23B, N23B
!pgi$ unroll = n:8
              DO i = S13B, N13B
                 vel(i,j,k,m) = vel(i,j,k,m) - gpre(i,j,k)
              END DO
           END DO
        END DO
     END IF
     
  END DO
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Massenfluss korrigieren ===============================================================================
  !===========================================================================================================
  IF (forcing_mode == 2) CALL force_massflow(.FALSE.)
  !===========================================================================================================
  
  
  END SUBROUTINE explicit
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE twostep
  
  IMPLICIT NONE
  
  INTEGER                ::  m
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Diese Zeitintegration hat trotz Runge-Kutta-Schema max. Konvergenzordnung eins, da hier   !
  !                Gp anstelle von H¯¹Gp gerechnet wird, d.h. es liegt ein Fehler in den Impulsgleichungen   !
  !                vor, abhängig von Delta t/(Re Delta x^2).                                                 !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !=== RHS fuer Druck-Gleichungssystem =======================================================================
  !--- Helmholtz-Gleichung lösen -----------------------------------------------------------------------------
  IF (rank == 0 .AND. log_iteration_yes) WRITE(10,'(a)') 'solving Helmholtz equation for RHS ...'
  
  DO m = 1, dimens
     CALL solve_Helmholtz(m,epsU,n_it_Helmh_vel,.TRUE.,rhs(b1L,b2L,b3L,m),vel(b1L,b2L,b3L,m),.TRUE.,.TRUE.)
  END DO
  
  !--- Divergenz bilden --------------------------------------------------------------------------------------
  CALL divergence2(vel,res)
  !===========================================================================================================
  
  
  !=== H-transformierten Druck bestimmen =====================================================================
  IF (rank == 0 .AND. log_iteration_yes) WRITE(10,'(a)')
  IF (rank == 0 .AND. log_iteration_yes) WRITE(10,'(a)') 'solving Poisson equation ...'
  
  ! Kanten/Ecken sind durch die Randbedingungen immer Divergenz-frei:
  IF (corner_yes   ) CALL handle_corner_rhs(1,res)
  
  CALL BiCGstab(epsU,n_it_Poisson,init_pre(substep),S1p,S2p,S3p,N1p,N2p,N3p,res,pre,2,.TRUE.,.TRUE.,precond_Poisson)
  !===========================================================================================================
  
  
  !=== Geschwindigkeiten bestimmen ===========================================================================
  DO m = 1, dimens
     CALL gradient(m,pre,gpre)
     
     IF (m == 1) vel(S11B:N11B,S21B:N21B,S31B:N31B,m) = vel(S11B:N11B,S21B:N21B,S31B:N31B,m) - gpre(S11B:N11B,S21B:N21B,S31B:N31B)
     IF (m == 2) vel(S12B:N12B,S22B:N22B,S32B:N32B,m) = vel(S12B:N12B,S22B:N22B,S32B:N32B,m) - gpre(S12B:N12B,S22B:N22B,S32B:N32B)
     IF (m == 3) vel(S13B:N13B,S23B:N23B,S33B:N33B,m) = vel(S13B:N13B,S23B:N23B,S33B:N33B,m) - gpre(S13B:N13B,S23B:N23B,S33B:N33B)
  END DO
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Massenfluss korrigieren ===============================================================================
  !===========================================================================================================
  IF (forcing_mode == 2) CALL force_massflow(.FALSE.)
  !===========================================================================================================
  
  
  END SUBROUTINE twostep
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE force_massflow(impl_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in)    ::  impl_yes
  
  INTEGER                ::  i, j, k, m
  REAL                   ::  flux(1:2), flux_global(1:2), dV
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Massenflusskorrektur kann Divergenzfreiheit nicht stoeren, daher ist nachtraegliche       !
  !                Korrektur moeglich:                                                                       !
  !                                                                                                          !
  ! m  = desired massflow                                                                                    !
  !                                                                                                          !
  ! A  = |H  G|                                                                                              !
  !      |D  0|                                                                                              !
  !                                                                                                          !
  ! A^-1 = |H^-1*(I-G*(D*H^-1*G)^-1*D*H^-1)  ... |                                                           !
  !        |          (D*H^-1*G)^-1*D*H^-1   ... |                                                           !
  ! x  = [u ,p]^T                                                                                            !
  ! b  = [f ,0]^T                                                                                            !
  ! j  = [1 ,0]^T                                                                                            !
  ! jt = [1t,0]^T                                                                                            !
  !                                                                                                          !
  ! |A  j|*|x| = |b|  ==>  |I     A^-1*j|*|x| = |   A^-1*b  |                                                !
  ! |jt 0| |k|   |m|       |0  jt*A^-1*j| |k|   |jt*A^-1*b-m|                                                !
  !                                                                                                          !
  ! s := jt*A^-1                                                                                             !
  ! t :=    A^-1*j                                                                                           !
  !                                                                                                          !
  ! ==> x = A^-1*b - t*(s*b-m)/(s*j)                                                                         !
  !                                                                                                          !
  !                                                                                                          !
  ! s = |1t*H^-1 ... |,  because jt*H^-1*G = 0 on equidistant grid in streamwise direction                   !
  !                                                                                                          !
  ! t = |H^-1*1|,        because  D*H^-1*j = 0                                                               !
  !     |  0   |                                                                                             !
  !                                                                                                          !
  ! ==> u = A^-1*f - H^-1*1*(1t*H^-1*f-m)/(1t*H^-1*1)                                                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  m = bulkflow_dir
  direction = m
  
  !===========================================================================================================
  IF (m == 1) THEN
     
     IF (impl_yes) THEN
        gpre = 1.
        
        IF (BC_2L .GT. 0) gpre(S11B:N11B,1 ,S31B:N31B) = 0.
        IF (BC_2U .GT. 0) gpre(S11B:N11B,N2,S31B:N31B) = 0.
        IF (BC_3L .GT. 0) gpre(S11B:N11B,S21B:N21B,1 ) = 0.
        IF (BC_3U .GT. 0) gpre(S11B:N11B,S21B:N21B,N3) = 0.
        
        CALL BiCGstab(epsU,n_it_Helmh_vel,init_vel(substep),S11B,S21B,S31B,N11B,N21B,N31B,gpre,dig,1,.TRUE.,.TRUE.,precond_Helmh_vel)
     ELSE
        dig = 1.
        
        IF (BC_2L .GT. 0) dig (S11B:N11B,1 ,S31B:N31B) = 0.
        IF (BC_2U .GT. 0) dig (S11B:N11B,N2,S31B:N31B) = 0.
        IF (BC_3L .GT. 0) dig (S11B:N11B,S21B:N21B,1 ) = 0.
        IF (BC_3U .GT. 0) dig (S11B:N11B,S21B:N21B,N3) = 0.
     END IF
     
     !--------------------------------------------------------------------------------------------------------
     !--- Massenfluss ----------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     ! Anmerkung: Eigentlich wuerde eine einzige Gitterpunktflaeche reichen, so ist es aber genauer.
     ! flux(1) = 1t*H^-1*1
     ! flux(2) =    H^-1*f
     
     flux = 0.
     
     DO k = S31B, N31B
        DO j = S21B, N21B
!pgi$ unroll = n:8
           DO i = S11B, N11B ! Anmerkung: dx1p(i) darf verwendet werden, weil Gitter ohnehin aequidistant (s. Anmerkungen!)
              dV = dx1p(i)*dx2p(j)*dx3p(k)
              flux(1) = flux(1) + dig(i,j,k  )*dV
              flux(2) = flux(2) + vel(i,j,k,m)*dV
           END DO
        END DO
     END DO
     
     CALL MPI_ALLREDUCE(flux(1:2),flux_global(1:2),2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     
     !--------------------------------------------------------------------------------------------------------
     !--- Forcing (k=(H^-1*f-m)/(1t*H^-1*1)) -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     flux(1) = (flux_global(2) - vel_bulk*L1*L2*L3) / flux_global(1) ! TEST!!! nur auf Kanalstroemung beschraenkt!!
     
     gpre(S11B:N11B,S21B:N21B,S31B:N31B) = flux(1)
     
     IF (BC_2L .GT. 0) gpre(S11B:N11B,1 ,S31B:N31B) = 0.
     IF (BC_2U .GT. 0) gpre(S11B:N11B,N2,S31B:N31B) = 0.
     IF (BC_3L .GT. 0) dig (S11B:N11B,S21B:N21B,1 ) = 0.
     IF (BC_3U .GT. 0) dig (S11B:N11B,S21B:N21B,N3) = 0.
     
     IF (impl_yes) THEN
        CALL BiCGstab(epsU,n_it_Helmh_vel,init_vel(substep),S11B,S21B,S31B,N11B,N21B,N31B,gpre,dig,1,.TRUE.,.TRUE.,precond_Helmh_vel)
        vel(S11B:N11B,S21B:N21B,S31B:N31B,m) = vel(S11B:N11B,S21B:N21B,S31B:N31B,m) - dig (S11B:N11B,S21B:N21B,S31B:N31B)
     ELSE
        vel(S11B:N11B,S21B:N21B,S31B:N31B,m) = vel(S11B:N11B,S21B:N21B,S31B:N31B,m) - gpre(S11B:N11B,S21B:N21B,S31B:N31B)
     END IF
     
  !===========================================================================================================
  ELSE IF (m == 2) THEN
     
     IF (impl_yes) THEN
        gpre = 1.
        
        IF (BC_1L .GT. 0) gpre(1 ,S22B:N22B,S32B:N32B) = 0.
        IF (BC_1U .GT. 0) gpre(N1,S22B:N22B,S32B:N32B) = 0.
        IF (BC_3L .GT. 0) gpre(S12B:N12B,S22B:N22B,1 ) = 0.
        IF (BC_3U .GT. 0) gpre(S12B:N12B,S22B:N22B,N3) = 0.
        
        CALL BiCGstab(epsU,n_it_Helmh_vel,init_vel(substep),S12B,S22B,S32B,N12B,N22B,N32B,gpre,dig,1,.TRUE.,.TRUE.,precond_Helmh_vel)
     ELSE
        dig = 1.
        
        IF (BC_1L .GT. 0) dig (1 ,S22B:N22B,S32B:N32B) = 0.
        IF (BC_1U .GT. 0) dig (N1,S22B:N22B,S32B:N32B) = 0.
        IF (BC_3L .GT. 0) dig (S12B:N12B,S22B:N22B,1 ) = 0.
        IF (BC_3U .GT. 0) dig (S12B:N12B,S22B:N22B,N3) = 0.
     END IF
     
     !--------------------------------------------------------------------------------------------------------
     !--- Massenfluss ----------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     ! Anmerkung: Eigentlich wuerde eine einzige Gitterpunktflaeche reichen, so ist es aber genauer.
     ! flux(1) = 1t*H^-1*1
     ! flux(2) =    H^-1*f
     
     flux = 0.
     
     DO k = S32B, N32B
        DO j = S22B, N22B
!pgi$ unroll = n:8
           DO i = S12B, N12B ! Anmerkung: dx2p(j) darf verwendet werden, weil Gitter ohnehin aequidistant (s. Anmerkungen!)
              dV = dx1p(i)*dx2p(j)*dx3p(k)
              flux(1) = flux(1) + dig(i,j,k  )*dV
              flux(2) = flux(2) + vel(i,j,k,m)*dV
           END DO
        END DO
     END DO
     
     CALL MPI_ALLREDUCE(flux(1:2),flux_global(1:2),2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     
     !--------------------------------------------------------------------------------------------------------
     !--- Forcing (k=(H^-1*f-m)/(1t*H^-1*1)) -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     flux(1) = (flux_global(2) - vel_bulk*L1*L2*L3) / flux_global(1) ! TEST!!! nur auf Kanalstroemung beschraenkt!!
     
     gpre(S12B:N12B,S22B:N22B,S32B:N32B) = flux(1)
     
     IF (BC_1L .GT. 0) gpre(1 ,S22B:N22B,S32B:N32B) = 0.
     IF (BC_1U .GT. 0) gpre(N1,S22B:N22B,S32B:N32B) = 0.
     IF (BC_3L .GT. 0) dig (S12B:N12B,S22B:N22B,1 ) = 0.
     IF (BC_3U .GT. 0) dig (S12B:N12B,S22B:N22B,N3) = 0.
     
     IF (impl_yes) THEN
        CALL BiCGstab(epsU,n_it_Helmh_vel,init_vel(substep),S12B,S22B,S32B,N12B,N22B,N32B,gpre,dig,1,.TRUE.,.TRUE.,precond_Helmh_vel)
        vel(S12B:N12B,S22B:N22B,S32B:N32B,m) = vel(S12B:N12B,S22B:N22B,S32B:N32B,m) - dig (S12B:N12B,S22B:N22B,S32B:N32B)
     ELSE
        vel(S12B:N12B,S22B:N22B,S32B:N32B,m) = vel(S12B:N12B,S22B:N22B,S32B:N32B,m) - gpre(S12B:N12B,S22B:N22B,S32B:N32B)
     END IF
     
  !===========================================================================================================
  ELSE IF (m == 3) THEN
     
     IF (impl_yes) THEN
        gpre = 1.
        
        IF (BC_1L .GT. 0) gpre(1 ,S23B:N23B,S33B:N33B) = 0.
        IF (BC_1U .GT. 0) gpre(N1,S23B:N23B,S33B:N33B) = 0.
        IF (BC_2L .GT. 0) gpre(S13B:N13B,1 ,S33B:N33B) = 0.
        IF (BC_2U .GT. 0) gpre(S13B:N13B,N2,S33B:N33B) = 0.
        
        CALL BiCGstab(epsU,n_it_Helmh_vel,init_vel(substep),S13B,S23B,S33B,N13B,N23B,N33B,gpre,dig,1,.TRUE.,.TRUE.,precond_Helmh_vel)
     ELSE
        dig = 1.
        
        IF (BC_1L .GT. 0) dig (1 ,S23B:N23B,S33B:N33B) = 0.
        IF (BC_1U .GT. 0) dig (N1,S23B:N23B,S33B:N33B) = 0.
        IF (BC_2L .GT. 0) dig (S13B:N13B,1 ,S33B:N33B) = 0.
        IF (BC_2U .GT. 0) dig (S13B:N13B,N2,S33B:N33B) = 0.
     END IF
     
     !--------------------------------------------------------------------------------------------------------
     !--- Massenfluss ----------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     ! Anmerkung: Eigentlich wuerde eine einzige Gitterpunktflaeche reichen, so ist es aber genauer.
     ! flux(1) = 1t*H^-1*1
     ! flux(2) =    H^-1*f
     
     flux = 0.
     
     DO k = S33B, N33B
        DO j = S23B, N23B
!pgi$ unroll = n:8
           DO i = S13B, N13B ! Anmerkung: dx3p(k) darf verwendet werden, weil Gitter ohnehin aequidistant (s. Anmerkungen!)
              dV = dx1p(i)*dx2p(j)*dx3p(k)
              flux(1) = flux(1) + dig(i,j,k  )*dV
              flux(2) = flux(2) + vel(i,j,k,m)*dV
           END DO
        END DO
     END DO
     
     CALL MPI_ALLREDUCE(flux(1:2),flux_global(1:2),2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     
     !--------------------------------------------------------------------------------------------------------
     !--- Forcing (k=(H^-1*f-m)/(1t*H^-1*1)) -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     flux(1) = (flux_global(2) - vel_bulk*L1*L2*L3) / flux_global(1) ! TEST!!! nur auf Kanalstroemung beschraenkt!!
     
     gpre(S13B:N13B,S23B:N23B,S33B:N33B) = flux(1)
     
     IF (BC_1L .GT. 0) gpre(1 ,S23B:N23B,S33B:N33B) = 0.
     IF (BC_1U .GT. 0) gpre(N1,S23B:N23B,S33B:N33B) = 0.
     IF (BC_2L .GT. 0) dig (S13B:N13B,1 ,S33B:N33B) = 0.
     IF (BC_2U .GT. 0) dig (S13B:N13B,N2,S33B:N33B) = 0.
     
     IF (impl_yes) THEN
        CALL BiCGstab(epsU,n_it_Helmh_vel,init_vel(substep),S13B,S23B,S33B,N13B,N23B,N33B,gpre,dig,1,.TRUE.,.TRUE.,precond_Helmh_vel)
        vel(S13B:N13B,S23B:N23B,S33B:N33B,m) = vel(S13B:N13B,S23B:N23B,S33B:N33B,m) - dig (S13B:N13B,S23B:N23B,S33B:N33B)
     ELSE
        vel(S13B:N13B,S23B:N23B,S33B:N33B,m) = vel(S13B:N13B,S23B:N23B,S33B:N33B,m) - gpre(S13B:N13B,S23B:N23B,S33B:N33B)
     END IF
     
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE force_massflow
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE apply_nullspace(res)
  
  IMPLICIT NONE
  
  REAL   , INTENT(inout) ::  res(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  flux, flux_global
  INTEGER                ::  i, j, k
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! note: - null space correction of pressure RHS and velocity RHS is currently not "consistent" (even the   !
  !         orthogonal projection is NOT equivalent!)                                                        !
  !       - null space correction of velocity RHS should make correction of pressure RHS obsolete            !
  !----------------------------------------------------------------------------------------------------------!
  
  
  flux = 0.
  
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           flux = flux + psi(i,j,k)*res(i,j,k)
        END DO
     END DO
  END DO
  
  CALL MPI_ALLREDUCE(flux,flux_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  flux = flux_global
  
  IF (rank == 0) WRITE(10,'(a,E25.17)') 'flux =', flux
  
  
  ! Orthogonale Projektion:
  res(S1p:N1p,S2p:N2p,S3p:N3p) = res(S1p:N1p,S2p:N2p,S3p:N3p) - flux*psi(S1p:N1p,S2p:N2p,S3p:N3p)
  
  
  END SUBROUTINE apply_nullspace
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE apply_nullspace2(g,psi,res,problem_type) ! TEST!!! Original ersetzen!
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  g
  
  REAL   , INTENT(in   ) ::  psi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL   , INTENT(inout) ::  res(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  INTEGER, INTENT(in   ) ::  problem_type
  
  REAL                   ::  flux, flux_global
  
  INTEGER                ::  S1R, N1R, i, dim1
  INTEGER                ::  S2R, N2R, j, dim2
  INTEGER                ::  S3R, N3R, k, dim3
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! note: - null space correction of pressure RHS and velocity RHS is currently not "consistent" (even the   !
  !         orthogonal projection is NOT equivalent!)                                                        !
  !       - also RHS correction of multigrid Poisson problem (finest grid) is normally slightly different    !
  !         from actual Poisson problem                                                                      !
  !----------------------------------------------------------------------------------------------------------!
  
  
  S1R  = SNB(1,1,g)
  S2R  = SNB(1,2,g)
  S3R  = SNB(1,3,g)
  
  N1R  = SNB(2,1,g)
  N2R  = SNB(2,2,g)
  N3R  = SNB(2,3,g)
  
  flux = 0.
  
  IF (problem_type == 2) THEN
     DO k = S3R, N3R
        DO j = S2R, N2R
           DO i = S1R, N1R
              flux = flux + psi(i,j,k)*res(i,j,k)
           END DO
        END DO
     END DO
  END IF
  
  IF (problem_type == 4 .OR. problem_type == 5) THEN
     DO k = S3R, N3R
        DO j = S2R, N2R
           DO i = S1R, N1R
              flux = flux + res(i,j,k) ! The left null space of the transposed matrix is a constant ...
           END DO
        END DO
     END DO
  END IF
  
  CALL MPI_ALLREDUCE(flux,flux_global,1,MPI_REAL8,MPI_SUM,comm1(g),merror)
  flux = flux_global
  
  !IF (rank == 0) WRITE(10,'(a,i3,a,E25.17)') 'grid #', g, ', flux =', flux
  
  
  ! Orthogonale Projektion:
  IF (problem_type == 2) THEN
     res(S1R:N1R,S2R:N2R,S3R:N3R) = res(S1R:N1R,S2R:N2R,S3R:N3R) - flux*psi(S1R:N1R,S2R:N2R,S3R:N3R)
  END IF
  
  IF (problem_type == 4 .OR. problem_type == 5) THEN
     dim1 = (NN(1,g)-1)*NB(1,g)+1
     dim2 = (NN(2,g)-1)*NB(2,g)+1
     dim3 = (NN(3,g)-1)*NB(3,g)+1
     IF (BC_1L_global == -1) dim1 = dim1-1
     IF (BC_2L_global == -1) dim2 = dim2-1
     IF (BC_3L_global == -1) dim3 = dim3-1
     
     flux = flux / REAL(dim1*dim2*dim3) ! hier keine Wurzel, weil oben auch nicht mit einem psi durchmultipliziert wird!
     
     res(S1R:N1R,S2R:N2R,S3R:N3R) = res(S1R:N1R,S2R:N2R,S3R:N3R) - flux
  END IF
  
  IF (corner_yes) CALL handle_corner_rhs(g,res) ! TEST!!! verifizieren ...
  
  
  END SUBROUTINE apply_nullspace2
  
  
  
  
  
  
  
  
  
  
  !> subroutine that first calculates the transpose of the div and grad stencils and then solves for the 
  !! nullspace.
  SUBROUTINE get_nullspace
  
  IMPLICIT NONE
  
  INTEGER                ::  g, m
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  REAL   , ALLOCATABLE   ::  work1(:,:,:)
  REAL   , ALLOCATABLE   ::  work2(:,:,:)
  
  REAL                   ::  eps, eps_global
  REAL                   ::  norm(1:3), norm_global(1:3)
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: -                                                                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! Richtig wäre eine Korrektur direkt von rhs, so dass <PSI,rhs> = 0: PSI = H^-1 D^T psi mit D H^-1 G psi^T = 0
  ! Allerdings: Viel zu aufwändig!
  
  
  !===========================================================================================================
  !=== transponierte Stencils bilden =========================================================================
  !===========================================================================================================
  cDu1T = 0.
  cDv2T = 0.
  cDw3T = 0.
  
  cGp1T = 0.
  cGp2T = 0.
  cGp3T = 0.
  !-----------------------------------------------------------------------------------------------------------
  ALLOCATE(work1(b1L:b1U,b1L:(M1+b1U),1:2))
  ALLOCATE(work2(b1L:b1U,b1L:(M1+b1U),1:2))
  
  work1 = 0.
  DO i = S1p, N1p
     work1(d1L:d1U,i+iShift,1) = cDu1(d1L:d1U,i)
  END DO
  DO i = S11B, N11B
     work1(g1L:g1U,i+iShift,2) = cGp1(g1L:g1U,i)
  END DO
  
  CALL MPI_ALLREDUCE(work1,work2,2*(b1U-b1L+1)*(M1+b1U-b1L+1),MPI_REAL8,MPI_SUM,COMM_BAR1,merror)
  
  IF (BC_1L_global == -1) THEN
     DO i = b1L, -1
        work2(b1L:b1U,2 +ls1+i,1) = work2(b1L:b1U,M1+ls1+1+i,1)
        work2(b1L:b1U,2 +ls1+i,2) = work2(b1L:b1U,M1+ls1+1+i,2)
     END DO
     DO i = 1, b1U
        work2(b1L:b1U,M1+ls1+i,1) = work2(b1L:b1U,1 +ls1  +i,1)
        work2(b1L:b1U,M1+ls1+i,2) = work2(b1L:b1U,1 +ls1  +i,2)
     END DO
  END IF
  
  DO i = S11B, N11B
     DO ii = g1L, g1U
        cDu1T(ii,i) = work2(-ii,i+ii+iShift,1)
     END DO
  END DO
  
  DO i = S1p, N1p
     DO ii = d1L, d1U
        cGp1T(ii,i) = work2(-ii,i+ii+iShift,2)
     END DO
  END DO
  
  DEALLOCATE(work1)
  DEALLOCATE(work2)
  !-----------------------------------------------------------------------------------------------------------
  ALLOCATE(work1(b2L:b2U,b2L:(M2+b2U),1:2))
  ALLOCATE(work2(b2L:b2U,b2L:(M2+b2U),1:2))
  
  work1 = 0.
  DO j = S2p, N2p
     work1(d2L:d2U,j+jShift,1) = cDv2(d2L:d2U,j)
  END DO
  DO j = S22B, N22B
     work1(g2L:g2U,j+jShift,2) = cGp2(g2L:g2U,j)
  END DO
  
  CALL MPI_ALLREDUCE(work1,work2,2*(b2U-b2L+1)*(M2+b2U-b2L+1),MPI_REAL8,MPI_SUM,COMM_BAR2,merror)
  
  IF (BC_2L_global == -1) THEN
     DO j = b2L, -1
        work2(b2L:b2U,2 +ls2+j,1) = work2(b2L:b2U,M2+ls2+1+j,1)
        work2(b2L:b2U,2 +ls2+j,2) = work2(b2L:b2U,M2+ls2+1+j,2)
     END DO
     DO j = 1, b2U
        work2(b2L:b2U,M2+ls2+j,1) = work2(b2L:b2U,1 +ls2  +j,1)
        work2(b2L:b2U,M2+ls2+j,2) = work2(b2L:b2U,1 +ls2  +j,2)
     END DO
  END IF
  
  DO j = S22B, N22B
     DO jj = g2L, g2U
        cDv2T(jj,j) = work2(-jj,j+jj+jShift,1)
     END DO
  END DO
  
  DO j = S2p, N2p
     DO jj = d2L, d2U
        cGp2T(jj,j) = work2(-jj,j+jj+jShift,2)
     END DO
  END DO
  
  DEALLOCATE(work1)
  DEALLOCATE(work2)
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
  
  ALLOCATE(work1(b3L:b3U,b3L:(M3+b3U),1:2))
  ALLOCATE(work2(b3L:b3U,b3L:(M3+b3U),1:2))
  
  work1 = 0.
  DO k = S3p, N3p
     work1(d3L:d3U,k+kShift,1) = cDw3(d3L:d3U,k)
  END DO
  DO k = S33B, N33B
     work1(g3L:g3U,k+kShift,2) = cGp3(g3L:g3U,k)
  END DO
  
  CALL MPI_ALLREDUCE(work1,work2,2*(b3U-b3L+1)*(M3+b3U-b3L+1),MPI_REAL8,MPI_SUM,COMM_BAR3,merror)
  
  IF (BC_3L_global == -1) THEN
     DO k = b3L, -1
        work2(b3L:b3U,2 +ls3+k,1) = work2(b3L:b3U,M3+ls3+1+k,1)
        work2(b3L:b3U,2 +ls3+k,2) = work2(b3L:b3U,M3+ls3+1+k,2)
     END DO
     DO k = 1, b3U
        work2(b3L:b3U,M3+ls3+k,1) = work2(b3L:b3U,1 +ls3  +k,1)
        work2(b3L:b3U,M3+ls3+k,2) = work2(b3L:b3U,1 +ls3  +k,2)
     END DO
  END IF
  
  DO k = S33B, N33B
     DO kk = g3L, g3U
        cDw3T(kk,k) = work2(-kk,k+kk+kShift,1)
     END DO
  END DO
  
  DO k = S3p, N3p
     DO kk = d3L, d3U
        cGp3T(kk,k) = work2(-kk,k+kk+kShift,2)
     END DO
  END DO
  
  DEALLOCATE(work1)
  DEALLOCATE(work2)
  
  END IF
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== Null-Raum lesen / berechnen ===========================================================================
  !===========================================================================================================
  IF (rank == 0 .AND. write_stout_yes) THEN
     IF (read_nullspace_yes) THEN
        WRITE(*,'(a)') 'reading left null-space of preconditioner ...'
     ELSE
        WRITE(*,'(a)') 'determining left null-space of preconditioner ...'
     END IF
  END IF
  
  ! rein empirisch:
  eps = 0.001*epsU
  
  IF (.NOT. read_nullspace_yes) THEN
     OPEN(10,FILE='test_log_iterations_nullspace.txt',STATUS='UNKNOWN')
     
     IF (rank == 0) WRITE(10,'(a,E27.15)') ' accuracy for null-space solution: eps =', eps
     IF (rank == 0) WRITE(10,*)
  END IF
  
  CALL solve_nullspace(eps,1,psi,4)
  
  IF (nullspace_coarse_yes) THEN
     IF (n_grids .GE. 1 ) CALL solve_nullspace(eps,1 ,psi_rel1 ,5)
     IF (n_grids .GE. 2 ) CALL solve_nullspace(eps,2 ,psi_rel2 ,5)
     IF (n_grids .GE. 3 ) CALL solve_nullspace(eps,3 ,psi_rel3 ,5)
     IF (n_grids .GE. 4 ) CALL solve_nullspace(eps,4 ,psi_rel4 ,5)
     IF (n_grids .GE. 5 ) CALL solve_nullspace(eps,5 ,psi_rel5 ,5)
     IF (n_grids .GE. 6 ) CALL solve_nullspace(eps,6 ,psi_rel6 ,5)
     IF (n_grids .GE. 7 ) CALL solve_nullspace(eps,7 ,psi_rel7 ,5)
     IF (n_grids .GE. 8 ) CALL solve_nullspace(eps,8 ,psi_rel8 ,5)
     IF (n_grids .GE. 9 ) CALL solve_nullspace(eps,9 ,psi_rel9 ,5)
     IF (n_grids .GE. 10) CALL solve_nullspace(eps,10,psi_rel10,5)
     IF (n_grids .GE. 11) CALL solve_nullspace(eps,11,psi_rel11,5)
     IF (n_grids .GE. 12) CALL solve_nullspace(eps,12,psi_rel12,5)
     IF (n_grids .GE. 13) CALL solve_nullspace(eps,13,psi_rel13,5)
     IF (n_grids .GE. 14) CALL solve_nullspace(eps,14,psi_rel14,5)
     IF (n_grids .GE. 15) CALL solve_nullspace(eps,15,psi_rel15,5)
  END IF
  
  IF (.NOT. read_nullspace_yes) CLOSE(10)
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== Null-Raum fuer vel ====================================================================================
  !===========================================================================================================
  ! Richtig wäre eine Korrektur direkt von rhs, so dass <PSI,rhs> = 0: PSI = H^-1 D^T psi mit D H^-1 G psi^T = 0
  ! Allerdings: Viel zu aufwändig!
  psi_vel = 0. ! ACHTUNG!!! Nur für timeint_mode == 1
  DO m = 1, dimens
     CALL divergence_transp  (m,psi,psi_vel(b1L,b2L,b3L,m))
     CALL bc_extrapolation_transp(m,psi_vel(b1L,b2L,b3L,m))
  END DO
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== Normieren =============================================================================================
  !===========================================================================================================
  rhs  = 0.
  IF (nullspace_ortho_yes) THEN
                      rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = psi_vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
                      rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = psi_vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
     IF (dimens == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = psi_vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
  ELSE
     ! Umspeichern auf rhs notwendig wegen Überlappung an Ecken/Kanten ...
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L .GT. 0) rhs(S11B:N11B,1        ,S31B:N31B,1) = th12(S11B:N11B,S31B:N31B,1)
     IF (BC_2U .GT. 0) rhs(S11B:N11B,N2       ,S31B:N31B,1) = th12(S11B:N11B,S31B:N31B,2)
     
     IF (BC_3L .GT. 0) rhs(S11B:N11B,S21B:N21B,1        ,1) = th13(S11B:N11B,S21B:N21B,1)
     IF (BC_3U .GT. 0) rhs(S11B:N11B,S21B:N21B,N3       ,1) = th13(S11B:N11B,S21B:N21B,2)
     
     IF (BC_1L .GT. 0) rhs(0        ,S21B:N21B,S31B:N31B,1) = th11(S21B:N21B,S31B:N31B,1)
     IF (BC_1U .GT. 0) rhs(N1       ,S21B:N21B,S31B:N31B,1) = th11(S21B:N21B,S31B:N31B,2)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) rhs(1        ,S22B:N22B,S32B:N32B,2) = th21(S22B:N22B,S32B:N32B,1)
     IF (BC_1U .GT. 0) rhs(N1       ,S22B:N22B,S32B:N32B,2) = th21(S22B:N22B,S32B:N32B,2)
     
     IF (BC_3L .GT. 0) rhs(S12B:N12B,S22B:N22B,1        ,2) = th23(S12B:N12B,S22B:N22B,1)
     IF (BC_3U .GT. 0) rhs(S12B:N12B,S22B:N22B,N3       ,2) = th23(S12B:N12B,S22B:N22B,2)
     
     IF (BC_2L .GT. 0) rhs(S12B:N12B,0        ,S32B:N32B,2) = th22(S12B:N12B,S32B:N32B,1)
     IF (BC_2U .GT. 0) rhs(S12B:N12B,N2       ,S32B:N32B,2) = th22(S12B:N12B,S32B:N32B,2)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) rhs(1        ,S23B:N23B,S33B:N33B,3) = th31(S23B:N23B,S33B:N33B,1)
     IF (BC_1U .GT. 0) rhs(N1       ,S23B:N23B,S33B:N33B,3) = th31(S23B:N23B,S33B:N33B,2)
     
     IF (BC_2L .GT. 0) rhs(S13B:N13B,1        ,S33B:N33B,3) = th32(S13B:N13B,S33B:N33B,1)
     IF (BC_2U .GT. 0) rhs(S13B:N13B,N2       ,S33B:N33B,3) = th32(S13B:N13B,S33B:N33B,2)
     
     IF (BC_3L .GT. 0) rhs(S13B:N13B,S23B:N23B,0        ,3) = th33(S13B:N13B,S23B:N23B,1)
     IF (BC_3U .GT. 0) rhs(S13B:N13B,S23B:N23B,N3       ,3) = th33(S13B:N13B,S23B:N23B,2)
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  norm = 0.
  
  m = 1
  DO k = S31B, N31B
     DO j = S21B, N21B
        DO i = S11B, N11B
           norm(1) = norm(1) + psi_vel(i,j,k,m)*rhs    (i,j,k,m)
           norm(2) = norm(2) + psi_vel(i,j,k,m)*psi_vel(i,j,k,m)
           norm(3) = norm(3) + rhs    (i,j,k,m)*rhs    (i,j,k,m)
        END DO
     END DO
  END DO
  
  m = 2
  DO k = S32B, N32B
     DO j = S22B, N22B
        DO i = S12B, N12B
           norm(1) = norm(1) + psi_vel(i,j,k,m)*rhs    (i,j,k,m)
           norm(2) = norm(2) + psi_vel(i,j,k,m)*psi_vel(i,j,k,m)
           norm(3) = norm(3) + rhs    (i,j,k,m)*rhs    (i,j,k,m)
        END DO
     END DO
  END DO
  
  IF (dimens == 3) THEN
  m = 3
  DO k = S33B, N33B
     DO j = S23B, N23B
        DO i = S13B, N13B
           norm(1) = norm(1) + psi_vel(i,j,k,m)*rhs    (i,j,k,m)
           norm(2) = norm(2) + psi_vel(i,j,k,m)*psi_vel(i,j,k,m)
           norm(3) = norm(3) + rhs    (i,j,k,m)*rhs    (i,j,k,m)
        END DO
     END DO
  END DO
  END IF
  
  CALL MPI_ALLREDUCE(norm,norm_global,3,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  
  norm(1) =      norm_global(1)
  norm(2) = SQRT(norm_global(2))
  norm(3) = SQRT(norm_global(3))
  
  IF (rank == 0) WRITE(*,*)
  IF (rank == 0) WRITE(*,*) '       |psi_vel| =', norm(2)
  IF (rank == 0) WRITE(*,*) '         |theta| =', norm(3)
  IF (rank == 0) WRITE(*,*) ' <psi_vel,theta> =', norm(1)
  IF (rank == 0) WRITE(*,*) ' (note that |<psi_vel,theta>| should not be much smaller than |psi_vel||theta|.)'
  IF (rank == 0) WRITE(*,*)
  
  IF (norm(1) == 0.) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! <psi_vel,theta> == 0.'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  norm(1) = 1./SQRT(norm_global(1))
  
  psi_vel = psi_vel * norm(1)
  
  IF (.NOT. nullspace_ortho_yes) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L .GT. 0) th12(S11B:N11B,S31B:N31B,1) = th12(S11B:N11B,S31B:N31B,1) * norm(1)
     IF (BC_2U .GT. 0) th12(S11B:N11B,S31B:N31B,2) = th12(S11B:N11B,S31B:N31B,2) * norm(1)
     
     IF (BC_3L .GT. 0) th13(S11B:N11B,S21B:N21B,1) = th13(S11B:N11B,S21B:N21B,1) * norm(1)
     IF (BC_3U .GT. 0) th13(S11B:N11B,S21B:N21B,2) = th13(S11B:N11B,S21B:N21B,2) * norm(1)
     
     IF (BC_1L .GT. 0) th11(S21B:N21B,S31B:N31B,1) = th11(S21B:N21B,S31B:N31B,1) * norm(1)
     IF (BC_1U .GT. 0) th11(S21B:N21B,S31B:N31B,2) = th11(S21B:N21B,S31B:N31B,2) * norm(1)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) th21(S22B:N22B,S32B:N32B,1) = th21(S22B:N22B,S32B:N32B,1) * norm(1)
     IF (BC_1U .GT. 0) th21(S22B:N22B,S32B:N32B,2) = th21(S22B:N22B,S32B:N32B,2) * norm(1)
     
     IF (BC_3L .GT. 0) th23(S12B:N12B,S22B:N22B,1) = th23(S12B:N12B,S22B:N22B,1) * norm(1)
     IF (BC_3U .GT. 0) th23(S12B:N12B,S22B:N22B,2) = th23(S12B:N12B,S22B:N22B,2) * norm(1)
     
     IF (BC_2L .GT. 0) th22(S12B:N12B,S32B:N32B,1) = th22(S12B:N12B,S32B:N32B,1) * norm(1)
     IF (BC_2U .GT. 0) th22(S12B:N12B,S32B:N32B,2) = th22(S12B:N12B,S32B:N32B,2) * norm(1)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) th31(S23B:N23B,S33B:N33B,1) = th31(S23B:N23B,S33B:N33B,1) * norm(1)
     IF (BC_1U .GT. 0) th31(S23B:N23B,S33B:N33B,2) = th31(S23B:N23B,S33B:N33B,2) * norm(1)
     
     IF (BC_2L .GT. 0) th32(S13B:N13B,S33B:N33B,1) = th32(S13B:N13B,S33B:N33B,1) * norm(1)
     IF (BC_2U .GT. 0) th32(S13B:N13B,S33B:N33B,2) = th32(S13B:N13B,S33B:N33B,2) * norm(1)
     
     IF (BC_3L .GT. 0) th33(S13B:N13B,S23B:N23B,1) = th33(S13B:N13B,S23B:N23B,1) * norm(1)
     IF (BC_3U .GT. 0) th33(S13B:N13B,S23B:N23B,2) = th33(S13B:N13B,S23B:N23B,2) * norm(1)
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE get_nullspace
  
  
  
  
  
  
  
  
  
  
  !> subroutine that solves the nullspace 
  SUBROUTINE solve_nullspace(eps,g,psi,problem_type)
  
  IMPLICIT NONE
  
  REAL   , INTENT(in   ) ::  eps
  INTEGER, INTENT(in   ) ::  g
  INTEGER, INTENT(in   ) ::  problem_type
  
  REAL   , INTENT(inout) ::  psi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  REAL                   ::  norm2, norm2_global
  INTEGER                ::  i, j, k
  INTEGER                ::  n_grids_orig
  
  REAL                   ::  line1(b1L:(N1+b1U),1:2)
  REAL                   ::  line2(b2L:(N2+b2U),1:2)
  REAL                   ::  line3(b3L:(N3+b3U),1:2)
  
  
  !--- KEIN Multigrid, nur Glätter auf feinem Gitter ---
  n_grids_orig = n_grids ! TEST!!!
  n_grids      = g
  
  
  IF (g == 1 .AND. problem_type == 4) THEN
     
     IF (1 == 2) THEN
        
        !=====================================================================================================
        IF (read_nullspace_yes) THEN
           !--- Lesen ---
           CALL read2_hdf('nullspace_fine'  ,'nullspace_fine'  ,S1p,S2p,S3p,N1p,N2p,N3p,0,          psi)
        ELSE
           !--- Lösen ---
           res = 0.
           psi = 1.
           CALL BiCGstab(eps,5000,.FALSE.,S1p,S2p,S3p,N1p,N2p,N3p,res,psi,problem_type,.TRUE.,.FALSE.,1)
           
           !--- Eckenbehandlung ---
           IF (corner_yes) CALL handle_corner_Lap(g,psi)
           
           !--- Normieren ---
           CALL get_norms2(g,N1,N2,N3,SNB(1,1,g),SNB(1,2,g),SNB(1,3,g),SNB(2,1,g),SNB(2,2,g),SNB(2,3,g),psi,problem_type,.FALSE.,.TRUE.,normTwo=norm2)
           psi = psi / SQRT(norm2)
           
           !--- Schreiben ---
           CALL write_hdf('nullspace_fine'  ,'nullspace_fine'  ,S1p,S2p,S3p,N1p,N2p,N3p,0,(/1,1,1/),psi)
           
           IF (rank == 0) WRITE(10,*)
           IF (rank == 0) CALL flush(10)
        END IF
        !=====================================================================================================
        
     ELSE
        
        !=====================================================================================================
        IF (read_nullspace_yes) THEN
           !--- Lesen ---
           CALL read2_1D_hdf('nullspace_fine_1','nullspace_fine_1',N1,S1p,N1p,iShift,1,line1(1,2))
        ELSE
           IF (rank == 0) WRITE(*,*) ' direction 1:'
           
           !--- Lösen ---
           line1(:,1) = 0.
           line1(:,2) = 1.
           CALL BiCGstab2_1D(eps,10*M1,1,b1L,b1U,N1,line1(b1L,1),line1(b1L,2),.TRUE.,.FALSE.,0) ! Achtung: Bei exakter Arithmetik muesste der Loeser die exakte Loesung nach spaetestens M1 Iterationen finden!
           
           !-- Normieren (soll under-/overflow vermeiden) ---
           norm2 = 0.
           DO i = SNB(1,1,g), SNB(2,1,g)
              norm2 = norm2 + line1(i,2)**2
           END DO
           CALL MPI_ALLREDUCE(norm2,norm2_global,1,MPI_REAL8,MPI_SUM,COMM_BAR1,merror)
           line1 = line1 / SQRT(norm2_global)
           
           !--- Schreiben ---
           CALL write_1D_hdf('nullspace_fine_1','nullspace_fine_1',N1,S1p,N1p,iShift,1,line1(1,2))
           
           IF (rank == 0) WRITE(10,*)
           IF (rank == 0) CALL flush(10)
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (read_nullspace_yes) THEN
           !--- Lesen ---
           CALL read2_1D_hdf('nullspace_fine_2','nullspace_fine_2',N2,S2p,N2p,jShift,2,line2(1,2))
        ELSE
           IF (rank == 0) WRITE(*,*) ' direction 2:'
           
           !--- Lösen ---
           line2(:,1) = 0.
           line2(:,2) = 1.
           CALL BiCGstab2_1D(eps,10*M2,2,b2L,b2U,N2,line2(b2L,1),line2(b2L,2),.TRUE.,.FALSE.,0)
           
           !-- Normieren (soll under-/overflow vermeiden) ---
           norm2 = 0.
           DO j = SNB(1,2,g), SNB(2,2,g)
              norm2 = norm2 + line2(j,2)**2
           END DO
           CALL MPI_ALLREDUCE(norm2,norm2_global,1,MPI_REAL8,MPI_SUM,COMM_BAR2,merror)
           line2 = line2 / SQRT(norm2_global)
           
           !--- Schreiben ---
           CALL write_1D_hdf('nullspace_fine_2','nullspace_fine_2',N2,S2p,N2p,jShift,2,line2(1,2))
           
           IF (rank == 0) WRITE(10,*)
           IF (rank == 0) CALL flush(10)
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (dimens == 3) THEN
        IF (read_nullspace_yes) THEN
           !--- Lesen ---
           CALL read2_1D_hdf('nullspace_fine_3','nullspace_fine_3',N3,S3p,N3p,kShift,3,line3(1,2))
        ELSE
           IF (rank == 0) WRITE(*,*) ' direction 3:'
           
           !--- Lösen ---
           line3(:,1) = 0.
           line3(:,2) = 1.
           CALL BiCGstab2_1D(eps,10*M3,3,b3L,b3U,N3,line3(b3L,1),line3(b3L,2),.TRUE.,.FALSE.,0)
           
           !-- Normieren (soll under-/overflow vermeiden) ---
           norm2 = 0.
           DO k = SNB(1,3,g), SNB(2,3,g)
              norm2 = norm2 + line3(k,2)**2
           END DO
           CALL MPI_ALLREDUCE(norm2,norm2_global,1,MPI_REAL8,MPI_SUM,COMM_BAR3,merror)
           line3 = line3 / SQRT(norm2_global)
           
           !--- Schreiben ---
           CALL write_1D_hdf('nullspace_fine_3','nullspace_fine_3',N3,S3p,N3p,kShift,3,line3(1,2))
           
           IF (rank == 0) WRITE(10,*)
           IF (rank == 0) CALL flush(10)
        END IF
        END IF
        !=====================================================================================================
        !--- "Dyadisches" Produkt bilden ---
        psi = 0.
        IF (dimens == 3) THEN
           DO k = SNB(1,3,g), SNB(2,3,g)
              DO j = SNB(1,2,g), SNB(2,2,g)
                 DO i = SNB(1,1,g), SNB(2,1,g)
                    psi(i,j,k) = line1(i,2)*line2(j,2)*line3(k,2)
                 END DO
              END DO
           END DO
        ELSE
           DO k = SNB(1,3,g), SNB(2,3,g)
              DO j = SNB(1,2,g), SNB(2,2,g)
                 DO i = SNB(1,1,g), SNB(2,1,g)
                    psi(i,j,k) = line1(i,2)*line2(j,2)
                 END DO
              END DO
           END DO
        END IF
        
        !--- Eckenbehandlung ---
        IF (corner_yes) CALL handle_corner_Lap(g,psi)
        !=====================================================================================================
     END IF
     
  ELSE IF (n_grids .GE. g) THEN
     
     IF (read_nullspace_yes .AND. g == 1) THEN
        !--- Lesen ---
        CALL read2_hdf('nullspace_coarse','nullspace_coarse',S1p,S2p,S3p,N1p,N2p,N3p,0,psi)
     ELSE
        !--- Lösen ---
        res = 0.
        psi = 1.
        CALL BiCGstab2(eps,200,.FALSE.,NN(1,g),NN(2,g),NN(3,g),g,SNB(1,1,g),SNB(1,2,g),SNB(1,3,g),SNB(2,1,g),SNB(2,2,g),SNB(2,3,g),res,psi,problem_type,.FALSE.,.TRUE.,1)
        
        !--- Eckenbehandlung ---
        IF (corner_yes) CALL handle_corner_Lap(g,psi)
        
        !--- Normieren ---
        CALL get_norms2(g,N1,N2,N3,SNB(1,1,g),SNB(1,2,g),SNB(1,3,g),SNB(2,1,g),SNB(2,2,g),SNB(2,3,g),psi,problem_type,.FALSE.,.TRUE.,normTwo=norm2)
        psi = psi / SQRT(norm2)
        
        !--- Schreiben ---
        IF (g == 1) CALL write_hdf('nullspace_coarse','nullspace_coarse',S1p,S2p,S3p,N1p,N2p,N3p,0,(/1,1,1/),psi)
        
        IF (rank == 0) WRITE(10,*)
        IF (rank == 0) CALL flush(10)
     END IF
     
  END IF
  
  
  n_grids = n_grids_orig ! TEST!!!
  
  
  END SUBROUTINE solve_nullspace
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! Verteilen auf andere Files, falls es sich bewaehrt ...
  !***********************************************************************
  !***********************************************************************
  SUBROUTINE BiCGstab2_1D(eps,n_it_max,m,bbL,bbU,Nmax,bb,phi,quiet_yes1,quiet_yes2,preconditioner)
  
  IMPLICIT NONE
  
  REAL   , INTENT(in)    ::  eps
  INTEGER, INTENT(in)    ::  n_it_max
  
  INTEGER, INTENT(in   ) ::  m
  INTEGER, INTENT(in   ) ::  bbL, bbU, Nmax
  
  INTEGER, INTENT(in)    ::  preconditioner
  
  REAL   , INTENT(in)    ::  bb (bbL:(Nmax+bbU))
  REAL   , INTENT(inout) ::  phi(bbL:(Nmax+bbU))
  
  REAL                   ::  pp (bbL:(Nmax+bbU))
  REAL                   ::  Ap (bbL:(Nmax+bbU))
  REAL                   ::  rr (bbL:(Nmax+bbU))
  REAL                   ::  rh (bbL:(Nmax+bbU))
  REAL                   ::  Ar (bbL:(Nmax+bbU))
  REAL                   ::  z1 (bbL:(Nmax+bbU))
  REAL                   ::  z2 (bbL:(Nmax+bbU))
  
  LOGICAL, INTENT(in)    ::  quiet_yes1
  LOGICAL, INTENT(in)    ::  quiet_yes2
  
  INTEGER                ::  i, counter
  INTEGER                ::  SS, NN
  
  INTEGER                ::  comm
  
  REAL                   ::  norm_inf, norm_inf_global
  REAL                   ::  rhr, rhr_prev, ArAr, rAr, rhAp, rhr_global, rhAp_global
  REAL                   ::  scalar_global(1:2)
  REAL                   ::  alpha, beta, omega
  
  LOGICAL                ::  exit_yes
  
  
  !===========================================================================================================
  IF      (m == 1) THEN
     SS = S1p
     NN = N1p
     comm = COMM_BAR1
  ELSE IF (m == 2) THEN
     SS = S2p
     NN = N2p
     comm = COMM_BAR2
  ELSE IF (m == 3) THEN
     SS = S3p
     NN = N3p
     comm = COMM_BAR3
  ELSE
     IF (rank == 0) WRITE(*,*) 'ERROR! m /= 1,2,3!'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  CALL product_div_grad_transp_1D(m,bbL,bbU,Nmax,phi,rr)
  
  rr(SS:NN) = bb(SS:NN) - rr(SS:NN)
  rh(SS:NN) = rr(SS:NN)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Residuum ==============================================================================================
  !===========================================================================================================
  norm_inf = 0.
  DO i = SS, NN
     norm_inf = MAX(ABS(rr(i)),norm_inf)
  END DO
  
  CALL MPI_ALLREDUCE(norm_inf,norm_inf_global,1,MPI_REAL8,MPI_MAX,comm,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
  norm_inf = norm_inf_global
  !===========================================================================================================
  
  counter = 0
  
  
  ITERATE: DO
     
     !========================================================================================================
     !=== Überprüfen des Konvergenzkriteriums ================================================================
     !========================================================================================================
     CALL status_iteration(eps,norm_inf,counter,n_it_max,exit_yes,quiet_yes1,quiet_yes2)
     IF (exit_yes) EXIT ITERATE
     !========================================================================================================
     
     
     !========================================================================================================
     !=== nächster Durchlauf =================================================================================
     !========================================================================================================
     counter = counter + 1
     !========================================================================================================
     
     
     !========================================================================================================
     rhr_prev = rhr
     rhr = 0.
     DO i = SS, NN
        rhr = rhr + rr(i)*rh(i)
     END DO
     
     CALL MPI_ALLREDUCE(rhr,rhr_global,1,MPI_REAL8,MPI_SUM,comm,merror)
     rhr = rhr_global
     !========================================================================================================
     IF (ABS(rhr) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rhr =', rhr
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'rhr =', rhr
        
        rh(SS:NN) = rr(SS:NN) ! Neuer Referenzvektor ...
        rhr = 0.
        DO i = SS, NN
           rhr = rhr + rr(i)*rh(i)
        END DO
        
        CALL MPI_ALLREDUCE(rhr,rhr_global,1,MPI_REAL8,MPI_SUM,comm,merror)
        rhr = rhr_global
     END IF
     !========================================================================================================
     IF (counter .GE. 2) THEN
        IF (omega == 0.) THEN
           IF (rank == 0) WRITE(* ,'(a,E13.5)') 'omega =', omega
           IF (rank == 0) WRITE(10,'(a,E13.5)') 'omega =', omega
           EXIT ITERATE
        END IF
        IF (rhr_prev == 0.) THEN
           IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rhr_prev =', rhr_prev
           IF (rank == 0) WRITE(10,'(a,E13.5)') 'rhr_prev =', rhr_prev
           EXIT ITERATE
        END IF
        beta = (alpha/omega)*(rhr/rhr_prev)
        omega = -beta*omega
        pp(SS:NN) = rr(SS:NN) + beta*pp(SS:NN) + omega*Ap(SS:NN)
     ELSE
        pp(SS:NN) = rr(SS:NN)
     END IF
     !========================================================================================================
     IF (preconditioner == 0) THEN
        CALL product_div_grad_transp_1D(m,bbL,bbU,Nmax,pp,Ap)
     ELSE IF (preconditioner == 1 .OR. preconditioner == 2) THEN
        !IF (preconditioner == 1) CALL multigridV(.TRUE.,g,pp,z1,problem_type)
        !IF (preconditioner == 2) CALL multigridF(.TRUE.,g,pp,z1,problem_type)
        CALL product_div_grad_transp_1D(m,bbL,bbU,Nmax,z1,Ap)
     ELSE
        IF (rank == 0) WRITE(*,'(a)') 'ERROR! Specify valid preconditioner!'
        CALL MPI_FINALIZE(merror)
        STOP
     END IF
     !========================================================================================================
     rhAp = 0.
     DO i = SS, NN
        rhAp = rhAp + rh(i)*Ap(i)
     END DO
     
     CALL MPI_ALLREDUCE(rhAp,rhAp_global,1,MPI_REAL8,MPI_SUM,comm,merror)
     rhAp = rhAp_global
     !========================================================================================================
     IF (ABS(rhAp) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rhAp =', rhAp
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'rhAp =', rhAp
        EXIT ITERATE
        alpha = 0.
     ELSE
        alpha = rhr / rhAp
     END IF
     !========================================================================================================
     rr(SS:NN) = rr(SS:NN) - alpha*Ap(SS:NN)
     !========================================================================================================
     IF (preconditioner == 0) THEN
        CALL product_div_grad_transp_1D(m,bbL,bbU,Nmax,rr,Ar)
     ELSE IF (preconditioner == 1 .OR. preconditioner == 2) THEN
        !IF (preconditioner == 1) CALL multigridV(.TRUE.,g,rr,z2,problem_type)
        !IF (preconditioner == 2) CALL multigridF(.TRUE.,g,rr,z2,problem_type)
        CALL product_div_grad_transp_1D(m,bbL,bbU,Nmax,z2,Ar)
     END IF
     !========================================================================================================
     rAr  = 0.
     ArAr = 0.
     DO i = SS, NN
        rAr  = rAr  + rr(i)*Ar(i)
        ArAr = ArAr + Ar(i)*Ar(i)
     END DO
     
     CALL MPI_ALLREDUCE((/rAr,ArAr/),scalar_global,2,MPI_REAL8,MPI_SUM,comm,merror)
     rAr  = scalar_global(1)
     ArAr = scalar_global(2)
     !========================================================================================================
     IF (ABS(rAr) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rAr =', rAr
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'rAr =', rAr
     END IF
     IF (ABS(ArAr) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'ArAr =', ArAr
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'ArAr =', ArAr
        EXIT ITERATE
        omega = 0.
     ELSE
        omega = rAr / ArAr
     END IF
     !========================================================================================================
     norm_inf = 0.
     !--------------------------------------------------------------------------------------------------------
     IF (preconditioner == 0) THEN
        DO i = SS, NN
           phi(i) = phi(i) + omega*rr(i) + alpha*pp(i)
           rr (i) = rr (i) - omega*Ar(i)
           norm_inf = MAX(ABS(rr(i)),norm_inf)
        END DO
     ELSE
        DO i = SS, NN
           phi(i) = phi(i) + omega*z2(i) + alpha*z1(i)
           rr (i) = rr (i) - omega*Ar(i)
           norm_inf = MAX(ABS(rr(i)),norm_inf)
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     CALL MPI_ALLREDUCE(norm_inf,norm_inf_global,1,MPI_REAL8,MPI_MAX,comm,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
     norm_inf = norm_inf_global
     !========================================================================================================
     
  END DO ITERATE
  
  
  END SUBROUTINE BiCGstab2_1D
  
  
  
  
  

 
  
  
  
  
  SUBROUTINE product_div_grad_transp_1D(m,bbL,bbU,Nmax,phi,Lap) ! TEST!!! woanders hin ...
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  m
  INTEGER, INTENT(in   ) ::  bbL, bbU, Nmax
  REAL   , INTENT(inout) ::  phi (bbL:(Nmax+bbU))
  REAL   , INTENT(out  ) ::  Lap (bbL:(Nmax+bbU))
  REAL                   ::  work(bbL:(Nmax+bbU))
  
  
  CALL divergence_transp_1D      (m,bbL,bbU,Nmax,phi,work)
  CALL bc_extrapolation_transp_1D(m,bbL,bbU,Nmax,work    )
  CALL gradient_transp_1D        (m,bbL,bbU,Nmax,work,Lap)
  
  
  END SUBROUTINE product_div_grad_transp_1D
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE bc_extrapolation_transp_1D(m,bbL,bbU,Nmax,phi) ! TEST!!! woanders hin ...
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  m
  INTEGER, INTENT(in   ) ::  bbL, bbU, Nmax
  REAL   , INTENT(inout) ::  phi(bbL:(Nmax+bbU))
  
  INTEGER                ::  i, ii
  
  
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 1) THEN
     IF (BC_1L .GT. 0) THEN
        i = 0
        DO ii = 0, d1U
           phi(i+ii+1) = phi(i+ii+1) - phi(i)*cIup(ii,1)/cIup(-1,1)
        END DO
        phi(i) = phi(i) / cIup(-1,1)
     END IF
     IF (BC_1U .GT. 0) THEN
        i = N1
        DO ii = d1L, -1
           phi(i+ii) = phi(i+ii) - phi(i)*cIup(ii,i)/cIup(0,i)
        END DO
        phi(i) = phi(i) / cIup(0,i)
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN
     IF (BC_2L .GT. 0) THEN
        i = 0
        DO ii = 0, d2U
           phi(i+ii+1) = phi(i+ii+1) - phi(i)*cIvp(ii,1)/cIvp(-1,1)
        END DO
        phi(i) = phi(i) / cIvp(-1,1)
     END IF
     IF (BC_2U .GT. 0) THEN
        i = N2
        DO ii = d2L, -1
           phi(i+ii) = phi(i+ii) - phi(i)*cIvp(ii,i)/cIvp(0,i)
        END DO
        phi(i) = phi(i) / cIvp(0,i)
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN
     IF (BC_3L .GT. 0) THEN
        i = 0
        DO ii = 0, d3U
           phi(i+ii+1) = phi(i+ii+1) - phi(i)*cIwp(ii,1)/cIwp(-1,1)
        END DO
        phi(i) = phi(i) / cIwp(-1,1)
     END IF
     IF (BC_3U .GT. 0) THEN
        i = N3
        DO ii = d3L, -1
           phi(i+ii) = phi(i+ii) - phi(i)*cIwp(ii,i)/cIwp(0,i)
        END DO
        phi(i) = phi(i) / cIwp(0,i)
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  
  
  END SUBROUTINE bc_extrapolation_transp_1D
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE divergence_transp_1D(m,bbL,bbU,Nmax,phi,div) ! TEST!!! woanders hin ...
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  m
  INTEGER, INTENT(in   ) ::  bbL, bbU, Nmax
  REAL   , INTENT(inout) ::  phi(bbL:(Nmax+bbU))
  REAL   , INTENT(inout) ::  div(bbL:(Nmax+bbU))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
        CALL exchange_1D(m,0,bbL,bbU,Nmax,phi)
        
        DO i = S11B, N11B
           div(i) = cDu1T(g1L,i)*phi(i+g1L)
           DO ii = g1L+1, g1U
              div(i) = div(i) + cDu1T(ii,i)*phi(i+ii)
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
        CALL exchange_1D(m,0,bbL,bbU,Nmax,phi)
        
        DO j = S22B, N22B
           div(j) = cDv2T(g2L,j)*phi(j+g2L)
           DO jj = g2L+1, g2U
              div(j) = div(j) + cDv2T(jj,j)*phi(j+jj)
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     !--------------------------------------------------------------------------------------------------------
        CALL exchange_1D(m,0,bbL,bbU,Nmax,phi)
        
        DO k = S33B, N33B
           div(k) = cDw3T(g3L,k)*phi(k+g3L)
           DO kk = g3L+1, g3U
              div(k) = div(k) + cDw3T(kk,k)*phi(k+kk)
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE divergence_transp_1D
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE gradient_transp_1D(m,bbL,bbU,Nmax,phi,grad) ! TEST!!! woanders hin ...
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  m
  INTEGER, INTENT(in   ) ::  bbL, bbU, Nmax
  REAL   , INTENT(inout) ::  phi (bbL:(Nmax+bbU))
  REAL   , INTENT(out  ) ::  grad(bbL:(Nmax+bbU))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) phi(0 ) = 0.
     IF (BC_1U .GT. 0) phi(N1) = 0.
     
     !--------------------------------------------------------------------------------------------------------
        CALL exchange_1D(m,m,bbL,bbU,Nmax,phi) ! Muss nach Randbedingungen kommen!
        
        DO i = S1p, N1p
           grad(i) = cGp1T(d1L,i)*phi(i+d1L)
           DO ii = d1L+1, d1U
              grad(i) = grad(i) + cGp1T(ii,i)*phi(i+ii)
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_2L .GT. 0) phi(0 ) = 0.
     IF (BC_2U .GT. 0) phi(N2) = 0.
     
     !--------------------------------------------------------------------------------------------------------
        CALL exchange_1D(m,m,bbL,bbU,Nmax,phi) ! Muss nach Randbedingungen kommen!
        
        DO j = S2p, N2p
           grad(j) = cGp2T(d2L,j)*phi(j+d2L)
           DO jj = d2L+1, d2U
              grad(j) = grad(j) + cGp2T(jj,j)*phi(j+jj)
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_3L .GT. 0) phi(0 ) = 0.
     IF (BC_3U .GT. 0) phi(N3) = 0.
     
     !--------------------------------------------------------------------------------------------------------
        CALL exchange_1D(m,m,bbL,bbU,Nmax,phi) ! Muss nach Randbedingungen kommen!
        
        DO k = S3p, N3p
           grad(k) = cGp3T(d3L,k)*phi(k+d3L)
           DO kk = d3L+1, d3U
              grad(k) = grad(k) + cGp3T(kk,k)*phi(k+kk)
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE gradient_transp_1D
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE exchange_1D(dir,vel_dir,bbL,bbU,Nmax,phi) ! TEST!!! woanders hin ...
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  dir
  INTEGER, INTENT(in   ) ::  vel_dir
  INTEGER, INTENT(in   ) ::  bbL, bbU, Nmax
  
  REAL   , INTENT(inout) ::  phi(bbL:(Nmax+bbU))
  
  INTEGER                ::  status(MPI_STATUS_SIZE) ! TEST!!! Warum steht das eigentlich hier lokal? gleiches gilt auch fuer die anderen Module ...
  
  REAL                   ::  ghostLR( 1:bbU)
  REAL                   ::  ghostLS( 1:bbU)
  REAL                   ::  ghostUR(bbL:-1)
  REAL                   ::  ghostUS(bbL:-1)
  
  INTEGER                ::  lengthL, lengthU
  INTEGER                ::  rankL, rankU
  INTEGER                ::  BCL, BCU
  INTEGER                ::  SS, NN
  
  
  lengthL = ABS(bbL)
  lengthU = ABS(bbU)
  
  IF      (dir == 1) THEN
     SS    = S1p
     NN    = N1p
     BCL   = BC_1L
     BCU   = BC_1U
     rankL = rank1L
     rankU = rank1U
  ELSE IF (dir == 2) THEN
     SS    = S2p
     NN    = N2p
     BCL   = BC_2L
     BCU   = BC_2U
     rankL = rank2L
     rankU = rank2U
  ELSE IF (dir == 3) THEN
     SS    = S3p
     NN    = N3p
     BCL   = BC_3L
     BCU   = BC_3U
     rankL = rank3L
     rankU = rank3U
  END IF
  
  
  !======================================================================================================================
  IF (BCL == 0) CALL MPI_IRECV(ghostUR,lengthL,MPI_REAL8,rankL,1,COMM_CART,req1L,merror) ! TEST!!! Was ist mit req1L, etc.?
  IF (BCU == 0) CALL MPI_IRECV(ghostLR,lengthU,MPI_REAL8,rankU,2,COMM_CART,req1U,merror)
  
  IF (BCU == 0) ghostUS(bbL:-1) = phi((NN+1+bbL):NN)
  IF (BCL == 0) ghostLS(1:bbU ) = phi(SS:(SS-1+bbU))
  
  IF (BCU == 0) CALL MPI_SEND(ghostUS,lengthL,MPI_REAL8,rankU,1,COMM_CART,merror)
  IF (BCL == 0) CALL MPI_SEND(ghostLS,lengthU,MPI_REAL8,rankL,2,COMM_CART,merror)
  
  IF (BCL == 0) CALL MPI_WAIT(req1L,status,merror)
  IF (BCU == 0) CALL MPI_WAIT(req1U,status,merror)
  
  IF (BCL == 0) CALL pseudocall(ghostUR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
  IF (BCU == 0) CALL pseudocall(ghostLR)
  
  IF (BCL == 0) phi((SS+bbL):(SS-1)) = ghostUR(bbL:-1)
  IF (BCU == 0) phi((NN+1):(NN+bbU)) = ghostLR(1:bbU )
  !-------------------------------------------------------------------------------------------------------------------
  IF (BCL == -1) THEN
     phi((SS+bbL):(SS-1)) = phi((NN+1+bbL):NN)
     phi((NN+1):(NN+bbU)) = phi(SS:(SS-1+bbU))
  END IF
  !-------------------------------------------------------------------------------------------------------------------
  IF (vel_dir == dir) THEN
     IF (BCL .GT. 0) phi( bbL    : -1       ) = 0.
     IF (BCU .GT. 0) phi((Nmax+1):(Nmax+bbU)) = 0.
     
     IF (BCL  == -2) phi( bbL    :  0       ) = 0.
     IF (BCU  == -2) phi( Nmax   :(Nmax+bbU)) = 0.
  ELSE
     IF (BCL .GT. 0 .OR. BCL == -2) phi( bbL    : 0        ) = 0.
     IF (BCU .GT. 0 .OR. BCU == -2) phi((Nmax+1):(Nmax+bbU)) = 0.
  END IF
  !======================================================================================================================
  
  
  END SUBROUTINE exchange_1D
  !***********************************************************************
  !***********************************************************************
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE solve_Helmholtz(m,epsU,n_it_max,init_yes,bb,phi,quiet_yes1,quiet_yes2)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  m
  
  REAL   , INTENT(in   ) ::  epsU
  INTEGER, INTENT(in   ) ::  n_it_max
  
  LOGICAL, INTENT(in   ) ::  init_yes
  
  REAL   , INTENT(in   ) ::  bb (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  LOGICAL, INTENT(in   ) ::  quiet_yes1
  LOGICAL, INTENT(in   ) ::  quiet_yes2
  
  
  IF (m == 1) THEN
     direction = 1
     IF (rank == 0 .AND. write_stout_yes .AND. .NOT. quiet_yes2) WRITE(*,'(a)') 'velocity component 1:'
     CALL BiCGstab(epsU,n_it_max,init_yes,S11B,S21B,S31B,N11B,N21B,N31B,bb,phi,1,quiet_yes1,quiet_yes2,precond_Helmh_vel)
  END IF
  
  IF (m == 2) THEN
     direction = 2
     IF (rank == 0 .AND. write_stout_yes .AND. .NOT. quiet_yes2) WRITE(*,'(a)') 'velocity component 2:'
     CALL BiCGstab(epsU,n_it_max,init_yes,S12B,S22B,S32B,N12B,N22B,N32B,bb,phi,1,quiet_yes1,quiet_yes2,precond_Helmh_vel)
  END IF
  
  IF (m == 3) THEN
     direction = 3
     IF (rank == 0 .AND. write_stout_yes .AND. .NOT. quiet_yes2) WRITE(*,'(a)') 'velocity component 3:'
     CALL BiCGstab(epsU,n_it_max,init_yes,S13B,S23B,S33B,N13B,N23B,N33B,bb,phi,1,quiet_yes1,quiet_yes2,precond_Helmh_vel)
  END IF
  
  
  END SUBROUTINE solve_Helmholtz
  
  
  
  
  
  
  
  
  
  !> performs a restriction step with a preliminary relaxation step
  SUBROUTINE relax_restrict(init_yes,nullspace_yes,g,psi,bb,phi,work,coarse,problem_type)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  init_yes
  LOGICAL, INTENT(in   ) ::  nullspace_yes
  INTEGER, INTENT(in   ) ::  g
  
  REAL   , INTENT(in   ) ::  psi   (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(inout) ::  bb    (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(inout) ::  phi   (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(inout) ::  work  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(out  ) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  
  INTEGER, INTENT(in   ) ::  problem_type
  
  
  !-----------------------------------------------------------------------------------------------------------
  IF (problem_type == 2 .OR. problem_type == 4 .OR. problem_type == 5) THEN
     IF (nullspace_yes) CALL apply_nullspace2(g,psi,bb,problem_type)
     CALL relaxation_div_grad(init_yes,n_relax_down,g,bb,phi)
     CALL product_div_grad_relax(g,phi,work)
     CALL restrict(.TRUE.,g,coarse,bb,work)
  END IF
  !-----------------------------------------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------------------------------------
  IF (problem_type == 1) THEN
     IF (g == 1) THEN
        CALL relaxation_Helmholtz(init_yes,n_relax_down,bb,phi)
        CALL product_Helmholtz_relax(phi,work)
        CALL restrict_Helmholtz(coarse,bb,work)
     ELSE
        CALL relaxation_Helmholtz_coarse(init_yes,n_relax_down,g,bb,phi)
        CALL product_Helmholtz_relax_coarse(g,phi,work)
        CALL restrict(.TRUE.,g,coarse,bb,work)
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------------------------------------
  ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
  !-----------------------------------------------------------------------------------------------------------
  
  
  END SUBROUTINE relax_restrict
  
  
  
  
  
  
  
  
  
  
  !> interpolation with a relaxation step at each interpolation 
  SUBROUTINE interpolate_relax(g,bb,phi,work,coarse,problem_type)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  g
  
  REAL   , INTENT(in   ) ::  bb    (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(inout) ::  phi   (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(out  ) ::  work  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(inout) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  
  INTEGER, INTENT(in   ) ::  problem_type
  
  
  !-----------------------------------------------------------------------------------------------------------
  IF (problem_type == 2 .OR. problem_type == 4 .OR. problem_type == 5) THEN
     CALL interpolate(.TRUE.,g,coarse,phi,work)
     CALL relaxation_div_grad_inv(.FALSE.,n_relax_up,g,bb,phi)
  END IF
  !-----------------------------------------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------------------------------------
  IF (problem_type == 1) THEN
     IF (g == 1) THEN
        CALL interpolate_Helmholtz(coarse,phi,work)
        CALL relaxation_Helmholtz       (.FALSE.,n_relax_up,  bb,phi)
     ELSE
        CALL interpolate(.TRUE.,g,coarse,phi,work)
        CALL relaxation_Helmholtz_coarse(.FALSE.,n_relax_up,g,bb,phi)
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------------------------------------
  ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
  !-----------------------------------------------------------------------------------------------------------
  
  END SUBROUTINE interpolate_relax
  
  
  
  
  
  
  
  
  
  
  !> performs a relaxation step at the coarsest level
  SUBROUTINE relax_bottom(init_yes,nullspace_yes,g,psi,bb,phi,problem_type)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  init_yes
  LOGICAL, INTENT(in   ) ::  nullspace_yes
  INTEGER, INTENT(in   ) ::  g
  
  REAL   , INTENT(in   ) ::  psi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL   , INTENT(inout) ::  bb (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL   , INTENT(inout) ::  phi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  INTEGER, INTENT(in   ) ::  problem_type
  
  
  !-----------------------------------------------------------------------------------------------------------
  IF (problem_type == 2 .OR. problem_type == 4 .OR. problem_type == 5) THEN
     IF (nullspace_yes) CALL apply_nullspace2(g,psi,bb,problem_type)
     CALL relaxation_div_grad    (init_yes,n_relax_bottom,g,bb,phi)
     CALL relaxation_div_grad_inv(.FALSE. ,n_relax_bottom,g,bb,phi)
  END IF
  !-----------------------------------------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------------------------------------
  IF (problem_type == 1) THEN
     IF (g  ==  1) CALL relaxation_Helmholtz       (init_yes,n_relax_bottom,  bb,phi)
     IF (g .GT. 1) CALL relaxation_Helmholtz_coarse(init_yes,n_relax_bottom,g,bb,phi)
  END IF
  !-----------------------------------------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------------------------------------
  ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
  !-----------------------------------------------------------------------------------------------------------
  
  
  END SUBROUTINE relax_bottom
  
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE plain_restrict(nullspace_yes,g,psi,bb,work,coarse,problem_type)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  nullspace_yes
  INTEGER, INTENT(in   ) ::  g
  
  REAL   , INTENT(in   ) ::  psi   (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(inout) ::  bb    (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(inout) ::  work  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(out  ) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  
  INTEGER, INTENT(in   ) ::  problem_type
  
  
  !-----------------------------------------------------------------------------------------------------------
  IF (problem_type == 2 .OR. problem_type == 4 .OR. problem_type == 5) THEN
     IF (nullspace_yes) CALL apply_nullspace2(g,psi,bb,problem_type)
     CALL restrict(.FALSE.,g,coarse,bb,work)
  END IF
  !-----------------------------------------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------------------------------------
  IF (problem_type == 1) THEN
     IF (g == 1) THEN
        CALL restrict_Helmholtz(coarse,bb,work)
     ELSE
        CALL restrict(.FALSE.,g,coarse,bb,work)
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------------------------------------
  ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
  !-----------------------------------------------------------------------------------------------------------
  
  
  END SUBROUTINE plain_restrict
  
  
  
  
  
  
  
  
  
  
  
  !> interpolation with a v-cycle at each interpolation step 
  SUBROUTINE interpolate_mg(g,bb,phi,work,coarse,problem_type)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  g
  
  REAL   , INTENT(inout) ::  bb    (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(inout) ::  phi   (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(inout) ::  work  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(inout) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  
  INTEGER, INTENT(in   ) ::  problem_type
  
  
  !-----------------------------------------------------------------------------------------------------------
  IF (problem_type == 2 .OR. problem_type == 4 .OR. problem_type == 5) THEN
     CALL interpolate(.FALSE.,g,coarse,work,phi)
     CALL multigridV (.FALSE.,g,bb,phi,problem_type)
  END IF
  !-----------------------------------------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------------------------------------
  IF (problem_type == 1) THEN
     IF (g == 1) THEN
        CALL interpolate_Helmholtz(coarse,phi,work)
        !!CALL relaxation_Helmholtz(.FALSE.,n_relax_up,  bb,phi)
        !CALL multigridV(.FALSE.,g,bb,phi,problem_type) ! TEST!!! Test schreiben, damit das nicht aufgerufen werden kann ...
     ELSE
        CALL interpolate(.FALSE.,g,coarse,work,phi)
        CALL multigridV (.FALSE.,g,bb,phi,problem_type)
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------------------------------------
  ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
  !-----------------------------------------------------------------------------------------------------------
  
  
  END SUBROUTINE interpolate_mg
  
  
  
  
  
  
  
  
  
  
  !> computes a V-cycle MG smoothing (preconditioning) 
  SUBROUTINE multigridV(init_yes,gstart,bb,phi,problem_type)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  init_yes
  INTEGER, INTENT(in   ) ::  gstart
  INTEGER, INTENT(in   ) ::  problem_type
  
  REAL   , INTENT(inout) ::  bb (b1L:(NN(1,gstart)+b1U),b2L:(NN(2,gstart)+b2U),b3L:(NN(3,gstart)+b3U))
  REAL   , INTENT(inout) ::  phi(b1L:(NN(1,gstart)+b1U),b2L:(NN(2,gstart)+b2U),b3L:(NN(3,gstart)+b3U))
  
  INTEGER                ::  g
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Übergebene Felder sind generell zu gross (b1L:(N1+b1U) vs. 0:(N1+1)), andererseits stört  !
  !                der Umstand i.A. auch nicht!                                                              !
  !              - vec1C bzw. vecC können eingespart werden, wenn Produkt und Restriktion zusammengelegt     !
  !                werden. Darauf wird hier der Übersicht halber verzichtet, zumal vecC nicht auf dem        !
  !                feinsten Gitter existiert und somit der Speicher nicht überstrapaziert wird. Diese Felder !
  !                werden zudem an die Interpolationsroutinen als Arbeitsfelder übergeben, wären aber auch   !
  !                dort mit entsprechendem Programmieraufwand eliminierbar.                                  !
  !              - INTENT(inout) für bb ist wegen restrict_Helmholtz notwendig.                              !
  !              - Es wird vorausgesetzt, dass das feinstes Gitterlevel keine Nullraum-Korrektur benötigt    !
  !                (sollte normalerweise auch erfüllt sein).                                                 !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  DO g = gstart, n_grids-1
     
     IF (participate_yes(g)) THEN
        IF (g == 1 )    CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel1 ,bb    ,phi   ,vec1C ,vec2A ,problem_type)
        IF (g == gstart) THEN
           IF (g == 2 ) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel2 ,bb    ,phi   ,vec2C ,vec3A ,problem_type)
           IF (g == 3 ) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel3 ,bb    ,phi   ,vec3C ,vec4A ,problem_type)
           IF (g == 4 ) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel4 ,bb    ,phi   ,vec4C ,vec5A ,problem_type)
           IF (g == 5 ) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel5 ,bb    ,phi   ,vec5C ,vec6A ,problem_type)
           IF (g == 6 ) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel6 ,bb    ,phi   ,vec6C ,vec7A ,problem_type)
           IF (g == 7 ) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel7 ,bb    ,phi   ,vec7C ,vec8A ,problem_type)
           IF (g == 8 ) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel8 ,bb    ,phi   ,vec8C ,vec9A ,problem_type)
           IF (g == 9 ) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel9 ,bb    ,phi   ,vec9C ,vec10A,problem_type)
           IF (g == 10) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel10,bb    ,phi   ,vec10C,vec11A,problem_type)
           IF (g == 11) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel11,bb    ,phi   ,vec11C,vec12A,problem_type)
           IF (g == 12) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel12,bb    ,phi   ,vec12C,vec13A,problem_type)
           IF (g == 13) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel13,bb    ,phi   ,vec13C,vec14A,problem_type)
           IF (g == 14) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel14,bb    ,phi   ,vec14C,vec15A,problem_type)
        ELSE
           IF (g == 2 ) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel2 ,vec2A ,vec2B ,vec2C ,vec3A ,problem_type)
           IF (g == 3 ) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel3 ,vec3A ,vec3B ,vec3C ,vec4A ,problem_type)
           IF (g == 4 ) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel4 ,vec4A ,vec4B ,vec4C ,vec5A ,problem_type)
           IF (g == 5 ) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel5 ,vec5A ,vec5B ,vec5C ,vec6A ,problem_type)
           IF (g == 6 ) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel6 ,vec6A ,vec6B ,vec6C ,vec7A ,problem_type)
           IF (g == 7 ) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel7 ,vec7A ,vec7B ,vec7C ,vec8A ,problem_type)
           IF (g == 8 ) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel8 ,vec8A ,vec8B ,vec8C ,vec9A ,problem_type)
           IF (g == 9 ) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel9 ,vec9A ,vec9B ,vec9C ,vec10A,problem_type)
           IF (g == 10) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel10,vec10A,vec10B,vec10C,vec11A,problem_type)
           IF (g == 11) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel11,vec11A,vec11B,vec11C,vec12A,problem_type)
           IF (g == 12) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel12,vec12A,vec12B,vec12C,vec13A,problem_type)
           IF (g == 13) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel13,vec13A,vec13B,vec13C,vec14A,problem_type)
           IF (g == 14) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel14,vec14A,vec14B,vec14C,vec15A,problem_type)
        END IF
     END IF
     
  END DO
  !===========================================================================================================
  
  !--- Grob-Gitter Lösung -------------------------------------
  IF (participate_yes(n_grids)) THEN
     IF (gstart == n_grids) THEN
                           CALL relax_bottom(init_yes,.FALSE.             ,n_grids,psi_rel1 ,bb    ,phi   ,problem_type) ! Achtung: psi_rel1 ist i.A. zu gross! (nullspace_coarse == .FALSE. <-- unschön!)
     ELSE
        IF (n_grids == 2 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel2 ,vec2A ,vec2B ,problem_type)
        IF (n_grids == 3 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel3 ,vec3A ,vec3B ,problem_type)
        IF (n_grids == 4 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel4 ,vec4A ,vec4B ,problem_type)
        IF (n_grids == 5 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel5 ,vec5A ,vec5B ,problem_type)
        IF (n_grids == 6 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel6 ,vec6A ,vec6B ,problem_type)
        IF (n_grids == 7 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel7 ,vec7A ,vec7B ,problem_type)
        IF (n_grids == 8 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel8 ,vec8A ,vec8B ,problem_type)
        IF (n_grids == 9 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel9 ,vec9A ,vec9B ,problem_type)
        IF (n_grids == 10) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel10,vec10A,vec10B,problem_type)
        IF (n_grids == 11) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel11,vec11A,vec11B,problem_type)
        IF (n_grids == 12) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel12,vec12A,vec12B,problem_type)
        IF (n_grids == 13) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel13,vec13A,vec13B,problem_type)
        IF (n_grids == 14) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel14,vec14A,vec14B,problem_type)
        IF (n_grids == 15) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel15,vec15A,vec15B,problem_type)
     END IF
  END IF
  
  !===========================================================================================================
  DO g = n_grids-1, gstart, -1
     
     IF (participate_yes(g)) THEN
        IF (g == gstart) THEN
           IF (g == 14) CALL interpolate_relax(g,bb    ,phi   ,vec14C,vec15B,problem_type)
           IF (g == 13) CALL interpolate_relax(g,bb    ,phi   ,vec13C,vec14B,problem_type)
           IF (g == 12) CALL interpolate_relax(g,bb    ,phi   ,vec12C,vec13B,problem_type)
           IF (g == 11) CALL interpolate_relax(g,bb    ,phi   ,vec11C,vec12B,problem_type)
           IF (g == 10) CALL interpolate_relax(g,bb    ,phi   ,vec10C,vec11B,problem_type)
           IF (g == 9 ) CALL interpolate_relax(g,bb    ,phi   ,vec9C ,vec10B,problem_type)
           IF (g == 8 ) CALL interpolate_relax(g,bb    ,phi   ,vec8C ,vec9B ,problem_type)
           IF (g == 7 ) CALL interpolate_relax(g,bb    ,phi   ,vec7C ,vec8B ,problem_type)
           IF (g == 6 ) CALL interpolate_relax(g,bb    ,phi   ,vec6C ,vec7B ,problem_type)
           IF (g == 5 ) CALL interpolate_relax(g,bb    ,phi   ,vec5C ,vec6B ,problem_type)
           IF (g == 4 ) CALL interpolate_relax(g,bb    ,phi   ,vec4C ,vec5B ,problem_type)
           IF (g == 3 ) CALL interpolate_relax(g,bb    ,phi   ,vec3C ,vec4B ,problem_type)
           IF (g == 2 ) CALL interpolate_relax(g,bb    ,phi   ,vec2C ,vec3B ,problem_type)
        ELSE
           IF (g == 14) CALL interpolate_relax(g,vec14A,vec14B,vec14C,vec15B,problem_type)
           IF (g == 13) CALL interpolate_relax(g,vec13A,vec13B,vec13C,vec14B,problem_type)
           IF (g == 12) CALL interpolate_relax(g,vec12A,vec12B,vec12C,vec13B,problem_type)
           IF (g == 11) CALL interpolate_relax(g,vec11A,vec11B,vec11C,vec12B,problem_type)
           IF (g == 10) CALL interpolate_relax(g,vec10A,vec10B,vec10C,vec11B,problem_type)
           IF (g == 9 ) CALL interpolate_relax(g,vec9A ,vec9B ,vec9C ,vec10B,problem_type)
           IF (g == 8 ) CALL interpolate_relax(g,vec8A ,vec8B ,vec8C ,vec9B ,problem_type)
           IF (g == 7 ) CALL interpolate_relax(g,vec7A ,vec7B ,vec7C ,vec8B ,problem_type)
           IF (g == 6 ) CALL interpolate_relax(g,vec6A ,vec6B ,vec6C ,vec7B ,problem_type)
           IF (g == 5 ) CALL interpolate_relax(g,vec5A ,vec5B ,vec5C ,vec6B ,problem_type)
           IF (g == 4 ) CALL interpolate_relax(g,vec4A ,vec4B ,vec4C ,vec5B ,problem_type)
           IF (g == 3 ) CALL interpolate_relax(g,vec3A ,vec3B ,vec3C ,vec4B ,problem_type)
           IF (g == 2 ) CALL interpolate_relax(g,vec2A ,vec2B ,vec2C ,vec3B ,problem_type)
        END IF
        IF (g == 1 )    CALL interpolate_relax(g,bb    ,phi   ,vec1C ,vec2B ,problem_type)
     END IF
     
  END DO
  !===========================================================================================================
  
  
  END SUBROUTINE multigridV
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE multigridF(init_yes,gstart,bb,phi,problem_type)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  init_yes
  INTEGER, INTENT(in   ) ::  gstart
  INTEGER, INTENT(in   ) ::  problem_type
  
  REAL   , INTENT(inout) ::  bb (b1L:(NN(1,gstart)+b1U),b2L:(NN(2,gstart)+b2U),b3L:(NN(3,gstart)+b3U))
  REAL   , INTENT(inout) ::  phi(b1L:(NN(1,gstart)+b1U),b2L:(NN(2,gstart)+b2U),b3L:(NN(3,gstart)+b3U))
  
  INTEGER                ::  g
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Es wird vorausgesetzt, dass das feinstes Gitterlevel keine Nullraum-Korrektur benötigt    !
  !                (sollte normalerweise auch erfüllt sein)                                                  !
  !              - Es wird immer initialisiert, mit Ausnahme von gstart == n_grids.                          !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  DO g = gstart, n_grids-1
     
     IF (participate_yes(g)) THEN
        IF (g == 1 )    CALL plain_restrict(.FALSE.             ,g,psi_rel1 ,bb    ,vec1C ,vec2A ,problem_type)
        IF (g == gstart) THEN
           IF (g == 2 ) CALL plain_restrict(.FALSE.             ,g,psi_rel2 ,bb    ,vec2C ,vec3A ,problem_type)
           IF (g == 3 ) CALL plain_restrict(.FALSE.             ,g,psi_rel3 ,bb    ,vec3C ,vec4A ,problem_type)
           IF (g == 4 ) CALL plain_restrict(.FALSE.             ,g,psi_rel4 ,bb    ,vec4C ,vec5A ,problem_type)
           IF (g == 5 ) CALL plain_restrict(.FALSE.             ,g,psi_rel5 ,bb    ,vec5C ,vec6A ,problem_type)
           IF (g == 6 ) CALL plain_restrict(.FALSE.             ,g,psi_rel6 ,bb    ,vec6C ,vec7A ,problem_type)
           IF (g == 7 ) CALL plain_restrict(.FALSE.             ,g,psi_rel7 ,bb    ,vec7C ,vec8A ,problem_type)
           IF (g == 8 ) CALL plain_restrict(.FALSE.             ,g,psi_rel8 ,bb    ,vec8C ,vec9A ,problem_type)
           IF (g == 9 ) CALL plain_restrict(.FALSE.             ,g,psi_rel9 ,bb    ,vec9C ,vec10A,problem_type)
           IF (g == 10) CALL plain_restrict(.FALSE.             ,g,psi_rel10,bb    ,vec10C,vec11A,problem_type)
           IF (g == 11) CALL plain_restrict(.FALSE.             ,g,psi_rel11,bb    ,vec11C,vec12A,problem_type)
           IF (g == 12) CALL plain_restrict(.FALSE.             ,g,psi_rel12,bb    ,vec12C,vec13A,problem_type)
           IF (g == 13) CALL plain_restrict(.FALSE.             ,g,psi_rel13,bb    ,vec13C,vec14A,problem_type)
           IF (g == 14) CALL plain_restrict(.FALSE.             ,g,psi_rel14,bb    ,vec14C,vec15A,problem_type)
        ELSE
           IF (g == 2 ) CALL plain_restrict(nullspace_coarse_yes,g,psi_rel2 ,vec2A ,vec2C ,vec3A ,problem_type)
           IF (g == 3 ) CALL plain_restrict(nullspace_coarse_yes,g,psi_rel3 ,vec3A ,vec3C ,vec4A ,problem_type)
           IF (g == 4 ) CALL plain_restrict(nullspace_coarse_yes,g,psi_rel4 ,vec4A ,vec4C ,vec5A ,problem_type)
           IF (g == 5 ) CALL plain_restrict(nullspace_coarse_yes,g,psi_rel5 ,vec5A ,vec5C ,vec6A ,problem_type)
           IF (g == 6 ) CALL plain_restrict(nullspace_coarse_yes,g,psi_rel6 ,vec6A ,vec6C ,vec7A ,problem_type)
           IF (g == 7 ) CALL plain_restrict(nullspace_coarse_yes,g,psi_rel7 ,vec7A ,vec7C ,vec8A ,problem_type)
           IF (g == 8 ) CALL plain_restrict(nullspace_coarse_yes,g,psi_rel8 ,vec8A ,vec8C ,vec9A ,problem_type)
           IF (g == 9 ) CALL plain_restrict(nullspace_coarse_yes,g,psi_rel9 ,vec9A ,vec9C ,vec10A,problem_type)
           IF (g == 10) CALL plain_restrict(nullspace_coarse_yes,g,psi_rel10,vec10A,vec10C,vec11A,problem_type)
           IF (g == 11) CALL plain_restrict(nullspace_coarse_yes,g,psi_rel11,vec11A,vec11C,vec12A,problem_type)
           IF (g == 12) CALL plain_restrict(nullspace_coarse_yes,g,psi_rel12,vec12A,vec12C,vec13A,problem_type)
           IF (g == 13) CALL plain_restrict(nullspace_coarse_yes,g,psi_rel13,vec13A,vec13C,vec14A,problem_type)
           IF (g == 14) CALL plain_restrict(nullspace_coarse_yes,g,psi_rel14,vec14A,vec14C,vec15A,problem_type)
        END IF
     END IF
     
  END DO
  !===========================================================================================================
  
  !--- Grob-Gitter Lösung -------------------------------------
  IF (participate_yes(n_grids)) THEN
     IF (gstart == n_grids) THEN
                           CALL relax_bottom(init_yes,.FALSE.             ,n_grids,psi_rel1 ,bb    ,phi   ,problem_type) ! Achtung: psi_rel1 ist i.A. zu gross! (nullspace_coarse == .FALSE. <-- unschön!)
     ELSE
        IF (n_grids == 2 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel2 ,vec2A ,vec2B ,problem_type)
        IF (n_grids == 3 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel3 ,vec3A ,vec3B ,problem_type)
        IF (n_grids == 4 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel4 ,vec4A ,vec4B ,problem_type)
        IF (n_grids == 5 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel5 ,vec5A ,vec5B ,problem_type)
        IF (n_grids == 6 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel6 ,vec6A ,vec6B ,problem_type)
        IF (n_grids == 7 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel7 ,vec7A ,vec7B ,problem_type)
        IF (n_grids == 8 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel8 ,vec8A ,vec8B ,problem_type)
        IF (n_grids == 9 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel9 ,vec9A ,vec9B ,problem_type)
        IF (n_grids == 10) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel10,vec10A,vec10B,problem_type)
        IF (n_grids == 11) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel11,vec11A,vec11B,problem_type)
        IF (n_grids == 12) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel12,vec12A,vec12B,problem_type)
        IF (n_grids == 13) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel13,vec13A,vec13B,problem_type)
        IF (n_grids == 14) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel14,vec14A,vec14B,problem_type)
        IF (n_grids == 15) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel15,vec15A,vec15B,problem_type)
     END IF
  END IF
  
  !===========================================================================================================
  DO g = n_grids-1, gstart, -1
     
     IF (participate_yes(g)) THEN
        IF (g == gstart) THEN
           IF (g == 14) CALL interpolate_mg(g,bb    ,phi   ,vec14C,vec15B,problem_type)
           IF (g == 13) CALL interpolate_mg(g,bb    ,phi   ,vec13C,vec14B,problem_type)
           IF (g == 12) CALL interpolate_mg(g,bb    ,phi   ,vec12C,vec13B,problem_type)
           IF (g == 11) CALL interpolate_mg(g,bb    ,phi   ,vec11C,vec12B,problem_type)
           IF (g == 10) CALL interpolate_mg(g,bb    ,phi   ,vec10C,vec11B,problem_type)
           IF (g == 9 ) CALL interpolate_mg(g,bb    ,phi   ,vec9C ,vec10B,problem_type)
           IF (g == 8 ) CALL interpolate_mg(g,bb    ,phi   ,vec8C ,vec9B ,problem_type)
           IF (g == 7 ) CALL interpolate_mg(g,bb    ,phi   ,vec7C ,vec8B ,problem_type)
           IF (g == 6 ) CALL interpolate_mg(g,bb    ,phi   ,vec6C ,vec7B ,problem_type)
           IF (g == 5 ) CALL interpolate_mg(g,bb    ,phi   ,vec5C ,vec6B ,problem_type)
           IF (g == 4 ) CALL interpolate_mg(g,bb    ,phi   ,vec4C ,vec5B ,problem_type)
           IF (g == 3 ) CALL interpolate_mg(g,bb    ,phi   ,vec3C ,vec4B ,problem_type)
           IF (g == 2 ) CALL interpolate_mg(g,bb    ,phi   ,vec2C ,vec3B ,problem_type)
        ELSE
           IF (g == 14) CALL interpolate_mg(g,vec14A,vec14B,vec14C,vec15B,problem_type)
           IF (g == 13) CALL interpolate_mg(g,vec13A,vec13B,vec13C,vec14B,problem_type)
           IF (g == 12) CALL interpolate_mg(g,vec12A,vec12B,vec12C,vec13B,problem_type)
           IF (g == 11) CALL interpolate_mg(g,vec11A,vec11B,vec11C,vec12B,problem_type)
           IF (g == 10) CALL interpolate_mg(g,vec10A,vec10B,vec10C,vec11B,problem_type)
           IF (g == 9 ) CALL interpolate_mg(g,vec9A ,vec9B ,vec9C ,vec10B,problem_type)
           IF (g == 8 ) CALL interpolate_mg(g,vec8A ,vec8B ,vec8C ,vec9B ,problem_type)
           IF (g == 7 ) CALL interpolate_mg(g,vec7A ,vec7B ,vec7C ,vec8B ,problem_type)
           IF (g == 6 ) CALL interpolate_mg(g,vec6A ,vec6B ,vec6C ,vec7B ,problem_type)
           IF (g == 5 ) CALL interpolate_mg(g,vec5A ,vec5B ,vec5C ,vec6B ,problem_type)
           IF (g == 4 ) CALL interpolate_mg(g,vec4A ,vec4B ,vec4C ,vec5B ,problem_type)
           IF (g == 3 ) CALL interpolate_mg(g,vec3A ,vec3B ,vec3C ,vec4B ,problem_type)
           IF (g == 2 ) CALL interpolate_mg(g,vec2A ,vec2B ,vec2C ,vec3B ,problem_type)
        END IF
        IF (g == 1 )    CALL interpolate_mg(g,bb    ,phi   ,vec1C ,vec2B ,problem_type)
     END IF
     
  END DO
  !===========================================================================================================
  
  
  END SUBROUTINE multigridF
  
  ! TEST!!! Hier lassen sich evtl. Operationen einsparen! (fine1-fine2) nur einmal berechnen, ist aber fraglich, ob das auch schneller ist ...
      SUBROUTINE restrict(add_yes,g,coarse,fine1,fine2) ! TEST!!! aufräumen ...
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  add_yes
  INTEGER, INTENT(in   ) ::  g
  
  REAL   , INTENT(out  ) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  REAL   , INTENT(inout) ::  fine1 (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(inout) ::  fine2 (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  
  INTEGER                ::  i, ii, di, imax, iimax
  INTEGER                ::  j, jj, dj, jmax, jjmax
  INTEGER                ::  k, kk, dk, kmax, kkmax
  
  INTEGER                ::  sizsg(1:3), offsg(1:3), dispg    
  REAL   , ALLOCATABLE   ::  sendbuf(:,:,:)
  REAL                   ::  recvbuf(1:NN(1,g+1)*NN(2,g+1)*NN(3,g+1))
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - für allgemeine di, dj, dk geeignet!                                                       !
  !              - überlappende Schicht in Blöcken wird (der Einfachheit halber) ebenfalls ausgetauscht, ist !
  !                aber im Prinzip redundant (genauer: coarse(S1R:N1R,S2R:N2R,S3R:N3R) = ...).               !
  !              - Motivation für diese kurze Routine ist die Möglichkeit, auch Varianten wie Full-Weighting !
  !                etc. ggf. einzubauen, ansonsten könnte sie auch eingespaart werden.                       !
  !              - Die Block-überlappenden Stirnflächen werden ebenfalls mitverarbeitet, aber eigentlich     !
  !                nicht gebraucht (erleichtert die Programmierung), so dass eine Initialisierung notwendig  !
  !                ist. Dies wiederum bedingt die INTENT(inout)-Deklaration.                                 !
  !----------------------------------------------------------------------------------------------------------!
  
  
  imax    = NN(1,g)
  jmax    = NN(2,g)
  kmax    = NN(3,g)
  
  iimax   = (NN(1,g+1)-1)/n_gather(1,g+1)+1
  jjmax   = (NN(2,g+1)-1)/n_gather(2,g+1)+1
  kkmax   = (NN(3,g+1)-1)/n_gather(3,g+1)+1
  
  di      = (NN(1,g)-1)/(iimax-1)
  dj      = (NN(2,g)-1)/(jjmax-1)
  dk      = (NN(3,g)-1)/(kkmax-1)
  
  
  IF (BC(1,1,g) .GT. 0 .AND. BC(1,2,g) .GT. 0) fine1(  1      ,  1      ,1:NN(3,g)) = (fine1(2        ,1      ,1:NN(3,g)) + fine1(1      ,2        ,1:NN(3,g)))/2.
  IF (BC(1,1,g) .GT. 0 .AND. BC(2,2,g) .GT. 0) fine1(  1      ,  NN(2,g),1:NN(3,g)) = (fine1(2        ,NN(2,g),1:NN(3,g)) + fine1(1      ,NN(2,g)-1,1:NN(3,g)))/2.
  IF (BC(2,1,g) .GT. 0 .AND. BC(1,2,g) .GT. 0) fine1(  NN(1,g),  1      ,1:NN(3,g)) = (fine1(NN(1,g)-1,1      ,1:NN(3,g)) + fine1(NN(1,g),2        ,1:NN(3,g)))/2.
  IF (BC(2,1,g) .GT. 0 .AND. BC(2,2,g) .GT. 0) fine1(  NN(1,g),  NN(2,g),1:NN(3,g)) = (fine1(NN(1,g)-1,NN(2,g),1:NN(3,g)) + fine1(NN(1,g),NN(2,g)-1,1:NN(3,g)))/2.
  
  IF (BC(1,1,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) fine1(  1      ,1:NN(2,g),  1      ) = (fine1(2        ,1:NN(2,g),1      ) + fine1(1      ,1:NN(2,g),2        ))/2.
  IF (BC(1,1,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) fine1(  1      ,1:NN(2,g),  NN(3,g)) = (fine1(2        ,1:NN(2,g),NN(3,g)) + fine1(1      ,1:NN(2,g),NN(3,g)-1))/2.
  IF (BC(2,1,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) fine1(  NN(1,g),1:NN(2,g),  1      ) = (fine1(NN(1,g)-1,1:NN(2,g),1      ) + fine1(NN(1,g),1:NN(2,g),2        ))/2.
  IF (BC(2,1,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) fine1(  NN(1,g),1:NN(2,g),  NN(3,g)) = (fine1(NN(1,g)-1,1:NN(2,g),NN(3,g)) + fine1(NN(1,g),1:NN(2,g),NN(3,g)-1))/2.

  IF (BC(1,2,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) fine1(1:NN(1,g),  1      ,  1      ) = (fine1(1:NN(1,g),2        ,1      ) + fine1(1:NN(1,g),1      ,2        ))/2.
  IF (BC(1,2,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) fine1(1:NN(1,g),  1      ,  NN(3,g)) = (fine1(1:NN(1,g),2        ,NN(3,g)) + fine1(1:NN(1,g),1      ,NN(3,g)-1))/2.
  IF (BC(2,2,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) fine1(1:NN(1,g),  NN(2,g),  1      ) = (fine1(1:NN(1,g),NN(2,g)-1,1      ) + fine1(1:NN(1,g),NN(2,g),2        ))/2.
  IF (BC(2,2,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) fine1(1:NN(1,g),  NN(2,g),  NN(3,g)) = (fine1(1:NN(1,g),NN(2,g)-1,NN(3,g)) + fine1(1:NN(1,g),NN(2,g),NN(3,g)-1))/2.
  
  
  IF (BC(1,1,g) .GT. 0 .AND. BC(1,2,g) .GT. 0) fine2(  1      ,  1      ,1:NN(3,g)) = (fine2(2        ,1      ,1:NN(3,g)) + fine2(1      ,2        ,1:NN(3,g)))/2.
  IF (BC(1,1,g) .GT. 0 .AND. BC(2,2,g) .GT. 0) fine2(  1      ,  NN(2,g),1:NN(3,g)) = (fine2(2        ,NN(2,g),1:NN(3,g)) + fine2(1      ,NN(2,g)-1,1:NN(3,g)))/2.
  IF (BC(2,1,g) .GT. 0 .AND. BC(1,2,g) .GT. 0) fine2(  NN(1,g),  1      ,1:NN(3,g)) = (fine2(NN(1,g)-1,1      ,1:NN(3,g)) + fine2(NN(1,g),2        ,1:NN(3,g)))/2.
  IF (BC(2,1,g) .GT. 0 .AND. BC(2,2,g) .GT. 0) fine2(  NN(1,g),  NN(2,g),1:NN(3,g)) = (fine2(NN(1,g)-1,NN(2,g),1:NN(3,g)) + fine2(NN(1,g),NN(2,g)-1,1:NN(3,g)))/2.
  
  IF (BC(1,1,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) fine2(  1      ,1:NN(2,g),  1      ) = (fine2(2        ,1:NN(2,g),1      ) + fine2(1      ,1:NN(2,g),2        ))/2.
  IF (BC(1,1,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) fine2(  1      ,1:NN(2,g),  NN(3,g)) = (fine2(2        ,1:NN(2,g),NN(3,g)) + fine2(1      ,1:NN(2,g),NN(3,g)-1))/2.
  IF (BC(2,1,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) fine2(  NN(1,g),1:NN(2,g),  1      ) = (fine2(NN(1,g)-1,1:NN(2,g),1      ) + fine2(NN(1,g),1:NN(2,g),2        ))/2.
  IF (BC(2,1,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) fine2(  NN(1,g),1:NN(2,g),  NN(3,g)) = (fine2(NN(1,g)-1,1:NN(2,g),NN(3,g)) + fine2(NN(1,g),1:NN(2,g),NN(3,g)-1))/2.

  IF (BC(1,2,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) fine2(1:NN(1,g),  1      ,  1      ) = (fine2(1:NN(1,g),2        ,1      ) + fine2(1:NN(1,g),1      ,2        ))/2.
  IF (BC(1,2,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) fine2(1:NN(1,g),  1      ,  NN(3,g)) = (fine2(1:NN(1,g),2        ,NN(3,g)) + fine2(1:NN(1,g),1      ,NN(3,g)-1))/2.
  IF (BC(2,2,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) fine2(1:NN(1,g),  NN(2,g),  1      ) = (fine2(1:NN(1,g),NN(2,g)-1,1      ) + fine2(1:NN(1,g),NN(2,g),2        ))/2.
  IF (BC(2,2,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) fine2(1:NN(1,g),  NN(2,g),  NN(3,g)) = (fine2(1:NN(1,g),NN(2,g)-1,NN(3,g)) + fine2(1:NN(1,g),NN(2,g),NN(3,g)-1))/2.
  
  
  IF (ls1 ==  0 .AND. (BC(1,1,g) == 0 .OR. BC(1,1,g) == -1)) iimax = iimax-1
  IF (ls1 == -1 .AND. (BC(2,1,g) == 0 .OR. BC(2,1,g) == -1)) iimax = iimax-1
  
  IF (ls2 ==  0 .AND. (BC(1,2,g) == 0 .OR. BC(1,2,g) == -1)) jjmax = jjmax-1
  IF (ls2 == -1 .AND. (BC(2,2,g) == 0 .OR. BC(2,2,g) == -1)) jjmax = jjmax-1
  
  IF (ls3 ==  0 .AND. (BC(1,3,g) == 0 .OR. BC(1,3,g) == -1)) kkmax = kkmax-1
  IF (ls3 == -1 .AND. (BC(2,3,g) == 0 .OR. BC(2,3,g) == -1)) kkmax = kkmax-1
  
  
  IF (1 == 2) THEN ! TEST!!!
     IF (add_yes) THEN
!pgi$ unroll = n:8
        coarse(1:iimax,1:jjmax,1:kkmax) = fine1(1:imax:di,1:jmax:dj,1:kmax:dk) - fine2(1:imax:di,1:jmax:dj,1:kmax:dk)
     ELSE
        coarse(1:iimax,1:jjmax,1:kkmax) = fine1(1:imax:di,1:jmax:dj,1:kmax:dk)
     END IF
     
  ELSE
     
     IF (add_yes) THEN ! TEST!!! etwas seriöser einbauen ...
     
     CALL exchange_relax(g,0,0,0,0,.TRUE.,fine1)
     CALL exchange_relax(g,0,0,0,0,.TRUE.,fine2)
     
     IF (dimens == 3) THEN
        DO kk = 1, kkmax
           k = dk*(kk-1)+1
           DO jj = 1, jjmax
              j = dj*(jj-1)+1
              DO ii = 1, iimax
                 i = di*(ii-1)+1
                 coarse(ii,jj,kk) = ((cR1(0,ii,g+1)+cR2(0,jj,g+1)+cR3(0,kk,g+1))*(fine1(i,j,k)-fine2(i,j,k)) +   &
                         &              cR1(-1,ii,g+1)*(fine1(i-1,j,k)-fine2(i-1,j,k)) +  &
                         &              cR1( 1,ii,g+1)*(fine1(i+1,j,k)-fine2(i+1,j,k)) +  &
                         &              cR2(-1,jj,g+1)*(fine1(i,j-1,k)-fine2(i,j-1,k)) +  &
                         &              cR2( 1,jj,g+1)*(fine1(i,j+1,k)-fine2(i,j+1,k)) +  &
                         &              cR3(-1,kk,g+1)*(fine1(i,j,k-1)-fine2(i,j,k-1)) +  &
                         &              cR3( 1,kk,g+1)*(fine1(i,j,k+1)-fine2(i,j,k+1))) / 3.
              END DO
           END DO
        END DO
     ELSE
        k  = 1
        kk = 1
        DO jj = 1, jjmax
           j = dj*(jj-1)+1
           DO ii = 1, iimax
              i = di*(ii-1)+1
              coarse(ii,jj,kk) = ((cR1(0,ii,g+1)+cR2(0,jj,g+1))*(fine1(i,j,k)-fine2(i,j,k)) +   &
                      &              cR1(-1,ii,g+1)*(fine1(i-1,j,k)-fine2(i-1,j,k)) +  &
                      &              cR1( 1,ii,g+1)*(fine1(i+1,j,k)-fine2(i+1,j,k)) +  &
                      &              cR2(-1,jj,g+1)*(fine1(i,j-1,k)-fine2(i,j-1,k)) +  &
                      &              cR2( 1,jj,g+1)*(fine1(i,j+1,k)-fine2(i,j+1,k))) / 2.
           END DO
        END DO
     END IF
     
     ELSE
     
     CALL exchange_relax(g,0,0,0,0,.TRUE.,fine1)
     
     IF (dimens == 3) THEN
        DO kk = 1, kkmax
           k = dk*(kk-1)+1
           DO jj = 1, jjmax
              j = dj*(jj-1)+1
              DO ii = 1, iimax
                 i = di*(ii-1)+1
                 coarse(ii,jj,kk) = ((cR1(0,ii,g+1)+cR2(0,jj,g+1)+cR3(0,kk,g+1))*fine1(i,j,k) +   &
                         &              cR1(-1,ii,g+1)*fine1(i-1,j,k) +  &
                         &              cR1( 1,ii,g+1)*fine1(i+1,j,k) +  &
                         &              cR2(-1,jj,g+1)*fine1(i,j-1,k) +  &
                         &              cR2( 1,jj,g+1)*fine1(i,j+1,k) +  &
                         &              cR3(-1,kk,g+1)*fine1(i,j,k-1) +  &
                         &              cR3( 1,kk,g+1)*fine1(i,j,k+1)) / 3.
              END DO
           END DO
        END DO
     ELSE
        k  = 1
        kk = 1
        DO jj = 1, jjmax
           j = dj*(jj-1)+1
           DO ii = 1, iimax
              i = di*(ii-1)+1
              coarse(ii,jj,kk) = ((cR1(0,ii,g+1)+cR2(0,jj,g+1))*fine1(i,j,k) +   &
                      &              cR1(-1,ii,g+1)*fine1(i-1,j,k) +  &
                      &              cR1( 1,ii,g+1)*fine1(i+1,j,k) +  &
                      &              cR2(-1,jj,g+1)*fine1(i,j-1,k) +  &
                      &              cR2( 1,jj,g+1)*fine1(i,j+1,k)) / 2.
           END DO
        END DO
     END IF
     
     END IF
     
  END IF
  
  
  
  IF (n_gather(1,g+1)*n_gather(2,g+1)*n_gather(3,g+1) .GT. 1) THEN
     
     ALLOCATE(sendbuf(1:iimax,1:jjmax,1:kkmax)) ! Anmerkung: Besser nicht fest allocieren um Speicherplatz zu sparen, ODER gleich "coarse" verwenden!
     
     DO kk = 1, kkmax
        DO jj = 1, jjmax
           DO ii = 1, iimax
              sendbuf(ii,jj,kk) = coarse(ii,jj,kk)
           END DO
        END DO
     END DO
     
     CALL MPI_GATHERv(sendbuf,iimax*jjmax*kkmax,MPI_REAL8,recvbuf,recvR(1,g+1),dispR(1,g+1),MPI_REAL8,rankc2(g+1),comm2(g+1),merror)
     
     DEALLOCATE(sendbuf)
     
     
     IF (participate_yes(g+1)) THEN
        DO k = 1, n_gather(3,g+1)
           DO j = 1, n_gather(2,g+1)
              DO i = 1, n_gather(1,g+1)
                 
                 sizsg(1:3) = sizsR(1:3,i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                 offsg(1:3) = offsR(1:3,i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                 dispg      = dispR(    i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                 
                 DO kk = 1, sizsg(3)
                    DO jj = 1, sizsg(2)
                       DO ii = 1, sizsg(1)
                          coarse(ii+offsg(1),jj+offsg(2),kk+offsg(3)) = recvbuf(dispg+ii+(jj-1)*sizsg(1)+(kk-1)*sizsg(1)*sizsg(2))
                       END DO
                    END DO
                 END DO
                 
              END DO
           END DO
        END DO
     END IF
     
  END IF
  
  
  END SUBROUTINE restrict
  
  
  SUBROUTINE restrict_Helmholtz(coarse,fine1,fine2)
  
  IMPLICIT NONE
  
  INTEGER, PARAMETER     ::  g = 1
  
  REAL   , INTENT(out  ) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  REAL   , INTENT(inout) ::  fine1 (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(inout) ::  fine2 (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  
  INTEGER                ::  i, ii, di
  INTEGER                ::  j, jj, dj
  INTEGER                ::  k, kk, dk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Null-Setzen am Rand nicht notwendig!                                                      !
  !              - Da nur in Richtung der jeweiligen Geschwindigkeitskomponente gemittelt wird, muss nicht   !
  !                die spezialisierte Helmholtz-Variante aufgerufen werden.                                  !
  !              - Austauschrichtung ist invers zu ex1, ex2, ex3. Bei mehreren Blöcken wird auch der jeweils !
  !                redundante "überlappende" Punkt aufgrund der zu grossen Intervallgrenzen (1:iimax) zwar   !
  !                berechnet, aber aufgrund des Einweg-Austauschs falsch berechnet! Dieses Vorgehen wurde    !
  !                bislang aus Übersichtsgründen vorgezogen, macht aber eine Initialisierung notwendig.      !
  !                Dabei werden Intervalle der Form 0:imax anstelle von 1:imax bearbeitet, da hier nur die   !
  !                das feinste Geschwindigkeitsgitter behandelt wird!                                        !
  !              - INTENT(inout) ist bei den feinen Gittern notwendig, da Ghost-Werte ausgetauscht werden    !
  !                müssen.                                                                                   !
  !              - Zuviele Daten werden ausgetauscht; eigentlich müsste in der Grenzfläche nur jeder 4.      !
  !                Punkt behandelt werden (4x zuviel!). Leider etwas unschön, könnte aber durch eine         !
  !                spezialisierte Austauschroutine behandelt werden, da das Übergeben von Feldern mit        !
  !                Intervallen von b1L:(iimax+b1U) nur sehr schlecht funktionieren würde (d.h. mit Um-       !
  !                kopieren).                                                                                !
  !----------------------------------------------------------------------------------------------------------!
  
  ! TEST!!! Test schreiben, um n_gather(:,2) .GT. 1 hier zu vermeiden! Gleiches gilt natürlich für die Interpolation.
  
  
  di = (NN(1,g)-1)/(NN(1,g+1)-1)
  dj = (NN(2,g)-1)/(NN(2,g+1)-1)
  dk = (NN(3,g)-1)/(NN(3,g+1)-1)
  
  CALL exchange_relax(g,-ex1,-ex2,-ex3,direction,.FALSE.,fine1)
  CALL exchange_relax(g,-ex1,-ex2,-ex3,direction,.FALSE.,fine2)
  
  
  IF (ls1 ==  0 .AND. (BC(1,1,g) == 0 .OR. BC(1,1,g) == -1)) fine1(1      ,0:NN(2,g),0:NN(3,g)) = 0.
  IF (ls1 ==  0 .AND. (BC(1,1,g) == 0 .OR. BC(1,1,g) == -1)) fine2(1      ,0:NN(2,g),0:NN(3,g)) = 0.
  IF (ls1 == -1 .AND. (BC(2,1,g) == 0 .OR. BC(2,1,g) == -1)) fine1(NN(1,g),0:NN(2,g),0:NN(3,g)) = 0.
  IF (ls1 == -1 .AND. (BC(2,1,g) == 0 .OR. BC(2,1,g) == -1)) fine2(NN(1,g),0:NN(2,g),0:NN(3,g)) = 0.
  
  IF (ls2 ==  0 .AND. (BC(1,2,g) == 0 .OR. BC(1,2,g) == -1)) fine1(0:NN(1,g),1      ,0:NN(3,g)) = 0.
  IF (ls2 ==  0 .AND. (BC(1,2,g) == 0 .OR. BC(1,2,g) == -1)) fine2(0:NN(1,g),1      ,0:NN(3,g)) = 0.
  IF (ls2 == -1 .AND. (BC(2,2,g) == 0 .OR. BC(2,2,g) == -1)) fine1(0:NN(1,g),NN(2,g),0:NN(3,g)) = 0.
  IF (ls2 == -1 .AND. (BC(2,2,g) == 0 .OR. BC(2,2,g) == -1)) fine2(0:NN(1,g),NN(2,g),0:NN(3,g)) = 0.
  
  IF (ls3 ==  0 .AND. (BC(1,3,g) == 0 .OR. BC(1,3,g) == -1)) fine1(0:NN(1,g),0:NN(2,g),1      ) = 0.
  IF (ls3 ==  0 .AND. (BC(1,3,g) == 0 .OR. BC(1,3,g) == -1)) fine2(0:NN(1,g),0:NN(2,g),1      ) = 0.
  IF (ls3 == -1 .AND. (BC(2,3,g) == 0 .OR. BC(2,3,g) == -1)) fine1(0:NN(1,g),0:NN(2,g),NN(3,g)) = 0.
  IF (ls3 == -1 .AND. (BC(2,3,g) == 0 .OR. BC(2,3,g) == -1)) fine2(0:NN(1,g),0:NN(2,g),NN(3,g)) = 0.
  
  
  !===========================================================================================================
  IF (direction == 1) THEN
     
     DO kk = 1, NN(3,g+1)
        k = dk*kk-1
        DO jj = 1, NN(2,g+1)
           j = dj*jj-1
           DO ii = 1, NN(1,g+1)
              i = di*ii-1
              coarse(ii,jj,kk) = cRH1(1,i)*(fine1(i-1,j,k)-fine2(i-1,j,k)) + cRH1(2,i)*(fine1(i,j,k)-fine2(i,j,k))
           END DO
        END DO
     END DO
     
  END IF
  !===========================================================================================================
  IF (direction == 2) THEN
     
     DO kk = 1, NN(3,g+1)
        k = dk*kk-1
        DO jj = 1, NN(2,g+1)
           j = dj*jj-1
           DO ii = 1, NN(1,g+1)
              i = di*ii-1
              coarse(ii,jj,kk) = cRH2(1,j)*(fine1(i,j-1,k)-fine2(i,j-1,k)) + cRH2(2,j)*(fine1(i,j,k)-fine2(i,j,k))
           END DO
        END DO
     END DO
     
  END IF
  !===========================================================================================================
  IF (direction == 3) THEN
     
     DO kk = 1, NN(3,g+1)
        k = dk*kk-1
        DO jj = 1, NN(2,g+1)
           j = dj*jj-1
           DO ii = 1, NN(1,g+1)
              i = di*ii-1
              coarse(ii,jj,kk) = cRH3(1,k)*(fine1(i,j,k-1)-fine2(i,j,k-1)) + cRH3(2,k)*(fine1(i,j,k)-fine2(i,j,k))
           END DO
        END DO
     END DO
     
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE restrict_Helmholtz
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE interpolate(add_yes,g,coarse,fine,work)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  add_yes
  INTEGER, INTENT(in   ) ::  g
  REAL   , INTENT(inout) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  REAL   , INTENT(inout) ::  fine  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(out  ) ::  work  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  
  INTEGER                ::  i, di, imax, iimax, iiShift
  INTEGER                ::  j, dj, jmax, jjmax, jjShift
  INTEGER                ::  k, dk, kmax, kkmax, kkShift
  
  !******************************************************************
  INTEGER                ::  ii, jj, kk
  INTEGER                ::  dispg, offsg(1:3)
  REAL   , ALLOCATABLE   ::  sendbuf(:,:,:)
  REAL                   ::  recvbuf(1:(NN(1,g+1)+NB(1,g)-1)*(NN(2,g+1)+NB(2,g)-1)*(NN(3,g+1)+NB(3,g)-1)) ! TEST!!! Ist das richtig so?
  !******************************************************************
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - für allgemeine di, dj, dk geeignet                                                        !
  !              - di /= 1 <===> N1 /= 1                                                                     !
  !              - es wird nur in eine Richung ausgetauscht                                                  !
  !              - Null-Setzen am Rand nicht notwendig                                                       !
  !              - Es wird sequentiell über alle Raumrichtungen interpoliert, um keinen individuellen        !
  !                Interpolationsstencil für jeden Punkt im Raum speichern zu müssen.                        !
  !              - Durch das sequentielle Interpolieren kann der Interpolationsstencil klein und damit der   !
  !                Gesamtaufwand minimiert werden (Alternative: 8- bzw. 26-Punkt Stencil (!)). Nachteilig    !
  !                ist dabei das zusätzliche Arbeitsfeld auf dem feineren Gitterniveau (wird der Multigrid-  !
  !                Routine entliehen).                                               .                       !
  !              - Interpolationskoeffizienten werden auf dem jeweils feineren Gitter gespeichert, um nicht  !
  !                auf die entsprechenden Indizes des gröberen Gitters umrechnen zu müssen.                  !
  !              - Die Block-überlappenden Stirnflächen werden ebenfalls mitverarbeitet, aber eigentlich     !
  !                nicht gebraucht (erleichtert die Programmierung), so dass eigentlich eine Initialisierung !
  !                notwendig wäre. Dies wird jedoch zuvor schon in der korrespondierenden Restriktions-      !
  !                Routine erledigt, so dass dies hier nicht mehr notwendig ist.                             !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! ACHTUNG!!! Verwendung von parametrischen Strides di,dj,dk in Schleifen verhindert die Vektorisierung bzw.
  !            das Prefetching! Andererseits ist der Geschwindigkeitsgewinn nur sehr gering (schon getestet).
  !
  ! - Geschwindigkeit ist trotz Prefetching / Vektorisierung stark durch die Speicherzugriffszeit limitiert
  ! - Das wird z.B. deutlich bei Single- vs. Dualcorebetrieb
  
  
  imax    = NN(1,g)
  jmax    = NN(2,g)
  kmax    = NN(3,g)
  
  iimax   = (NN(1,g+1)-1)/n_gather(1,g+1)+1
  jjmax   = (NN(2,g+1)-1)/n_gather(2,g+1)+1
  kkmax   = (NN(3,g+1)-1)/n_gather(3,g+1)+1
  
  di      = (NN(1,g)-1)/(iimax-1)
  dj      = (NN(2,g)-1)/(jjmax-1)
  dk      = (NN(3,g)-1)/(kkmax-1)
  
  iiShift = (iimax-1)*MOD(iB(1,g)-1,n_gather(1,g+1))
  jjShift = (jjmax-1)*MOD(iB(2,g)-1,n_gather(2,g+1))
  kkShift = (kkmax-1)*MOD(iB(3,g)-1,n_gather(3,g+1))
  
  
  IF (BC(1,1,g+1) .GT. 0 .AND. BC(1,2,g+1) .GT. 0) coarse(1        ,1        ,1:NN(3,g+1)) = (coarse(2          ,1        ,1:NN(3,g+1)) + coarse(1        ,2          ,1:NN(3,g+1)) + coarse(2          ,2          ,1:NN(3,g+1)))/3.
  IF (BC(1,1,g+1) .GT. 0 .AND. BC(2,2,g+1) .GT. 0) coarse(1        ,NN(2,g+1),1:NN(3,g+1)) = (coarse(2          ,NN(2,g+1),1:NN(3,g+1)) + coarse(1        ,NN(2,g+1)-1,1:NN(3,g+1)) + coarse(2          ,NN(2,g+1)-1,1:NN(3,g+1)))/3.
  IF (BC(2,1,g+1) .GT. 0 .AND. BC(1,2,g+1) .GT. 0) coarse(NN(1,g+1),1        ,1:NN(3,g+1)) = (coarse(NN(1,g+1)-1,1        ,1:NN(3,g+1)) + coarse(NN(1,g+1),2          ,1:NN(3,g+1)) + coarse(NN(1,g+1)-1,2          ,1:NN(3,g+1)))/3.
  IF (BC(2,1,g+1) .GT. 0 .AND. BC(2,2,g+1) .GT. 0) coarse(NN(1,g+1),NN(2,g+1),1:NN(3,g+1)) = (coarse(NN(1,g+1)-1,NN(2,g+1),1:NN(3,g+1)) + coarse(NN(1,g+1),NN(2,g+1)-1,1:NN(3,g+1)) + coarse(NN(1,g+1)-1,NN(2,g+1)-1,1:NN(3,g+1)))/3.
  
  IF (BC(1,1,g+1) .GT. 0 .AND. BC(1,3,g+1) .GT. 0) coarse(1        ,1:NN(2,g+1),1        ) = (coarse(2          ,1:NN(2,g+1),1        ) + coarse(1        ,1:NN(2,g+1),2          ) + coarse(2          ,1:NN(2,g+1),2          ))/3.
  IF (BC(1,1,g+1) .GT. 0 .AND. BC(2,3,g+1) .GT. 0) coarse(1        ,1:NN(2,g+1),NN(3,g+1)) = (coarse(2          ,1:NN(2,g+1),NN(3,g+1)) + coarse(1        ,1:NN(2,g+1),NN(3,g+1)-1) + coarse(2          ,1:NN(2,g+1),NN(3,g+1)-1))/3.
  IF (BC(2,1,g+1) .GT. 0 .AND. BC(1,3,g+1) .GT. 0) coarse(NN(1,g+1),1:NN(2,g+1),1        ) = (coarse(NN(1,g+1)-1,1:NN(2,g+1),1        ) + coarse(NN(1,g+1),1:NN(2,g+1),2          ) + coarse(NN(1,g+1)-1,1:NN(2,g+1),2          ))/3.
  IF (BC(2,1,g+1) .GT. 0 .AND. BC(2,3,g+1) .GT. 0) coarse(NN(1,g+1),1:NN(2,g+1),NN(3,g+1)) = (coarse(NN(1,g+1)-1,1:NN(2,g+1),NN(3,g+1)) + coarse(NN(1,g+1),1:NN(2,g+1),NN(3,g+1)-1) + coarse(NN(1,g+1)-1,1:NN(2,g+1),NN(3,g+1)-1))/3.

  IF (BC(1,2,g+1) .GT. 0 .AND. BC(1,3,g+1) .GT. 0) coarse(1:NN(1,g+1),1        ,1        ) = (coarse(1:NN(1,g+1),2          ,1        ) + coarse(1:NN(1,g+1),1        ,2          ) + coarse(1:NN(1,g+1),2          ,2          ))/3.
  IF (BC(1,2,g+1) .GT. 0 .AND. BC(2,3,g+1) .GT. 0) coarse(1:NN(1,g+1),1        ,NN(3,g+1)) = (coarse(1:NN(1,g+1),2          ,NN(3,g+1)) + coarse(1:NN(1,g+1),1        ,NN(3,g+1)-1) + coarse(1:NN(1,g+1),2          ,NN(3,g+1)-1))/3.
  IF (BC(2,2,g+1) .GT. 0 .AND. BC(1,3,g+1) .GT. 0) coarse(1:NN(1,g+1),NN(2,g+1),1        ) = (coarse(1:NN(1,g+1),NN(2,g+1)-1,1        ) + coarse(1:NN(1,g+1),NN(2,g+1),2          ) + coarse(1:NN(1,g+1),NN(2,g+1)-1,2          ))/3.
  IF (BC(2,2,g+1) .GT. 0 .AND. BC(2,3,g+1) .GT. 0) coarse(1:NN(1,g+1),NN(2,g+1),NN(3,g+1)) = (coarse(1:NN(1,g+1),NN(2,g+1)-1,NN(3,g+1)) + coarse(1:NN(1,g+1),NN(2,g+1),NN(3,g+1)-1) + coarse(1:NN(1,g+1),NN(2,g+1)-1,NN(3,g+1)-1))/3.
  
  
  CALL exchange_relax(g+1,ex1,ex2,ex3,0,.FALSE.,coarse) ! Anmerkung: .FALSE. ist ok ...
  
  
  !***********************************************************************************************************
  IF (n_gather(1,g+1)*n_gather(2,g+1)*n_gather(3,g+1) .GT. 1) THEN
     
     IF (participate_yes(g+1)) THEN
        DO k = 1, n_gather(3,g+1)
           DO j = 1, n_gather(2,g+1)
              DO i = 1, n_gather(1,g+1)
                 
                 dispg      = dispI(    i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                 offsg(1:3) = offsI(1:3,i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                 
                 DO kk = 1, kkmax
                    DO jj = 1, jjmax
                       DO ii = 1, iimax
                          recvbuf(dispg+ii+(jj-1)*iimax+(kk-1)*iimax*jjmax) = coarse(ii+offsg(1),jj+offsg(2),kk+offsg(3))
                       END DO
                    END DO
                 END DO
                 
              END DO
           END DO
        END DO
     END IF
     
     
     ALLOCATE(sendbuf(1:iimax,1:jjmax,1:kkmax)) ! Anmerkung: Besser nicht fest allocieren um Speicherplatz zu sparen, ODER gleich "coarse" verwenden!
     
     CALL MPI_SCATTER(recvbuf,iimax*jjmax*kkmax,MPI_REAL8,sendbuf,iimax*jjmax*kkmax,MPI_REAL8,rankc2(g+1),comm2(g+1),merror)
     
     DO kk = 1, kkmax
        k = dk*(kk-1)+1
        DO jj = 1, jjmax
           j = dj*(jj-1)+1
!pgi$ unroll = n:8
           DO ii = 1, iimax
              i = di*(ii-1)+1
              work(i,j,k) = sendbuf(ii,jj,kk)
           END DO
        END DO
     END DO
     
     DEALLOCATE(sendbuf)
     
  ELSE
     
!pgi$ unroll = n:8
     work(1:imax:di,1:jmax:dj,1:kmax:dk) = coarse((1+iiShift):(iimax+iiShift),(1+jjShift):(jjmax+jjShift),(1+kkShift):(kkmax+kkShift))
     
  END IF
  !***********************************************************************************************************
  
  
  
  !===========================================================================================================
  IF (dk /= 1) THEN ! (dimens == 2) <==> (dk == 1) automatisch erfüllt!
     
     DO k = 2, kmax-1, dk
        DO j = 1, jmax, dj
!pgi$ unroll = n:8
           DO i = 1, imax, di
              work(i,j,k) = cI3(1,k,g)*work(i,j,k-1) + cI3(2,k,g)*work(i,j,k+1)
           END DO
        END DO
     END DO
     
  END IF
  !===========================================================================================================
  IF (dj /= 1) THEN ! TEST!!! in 2D wird hier doppelte Arbeit geleistet! (kmax == 2??)
     
     DO k = 1, kmax
        DO j = 2, jmax-1, dj
!pgi$ unroll = n:8
           DO i = 1, imax, di
              work(i,j,k) = cI2(1,j,g)*work(i,j-1,k) + cI2(2,j,g)*work(i,j+1,k)
           END DO
        END DO
     END DO
     
  END IF
  !===========================================================================================================
  IF (di /= 1) THEN
     
     DO k = 1, kmax
        DO j = 1, jmax
!pgi$ unroll = n:8
           DO i = 2, imax-1, di
              work(i,j,k) = cI1(1,i,g)*work(i-1,j,k) + cI1(2,i,g)*work(i+1,j,k)
           END DO
        END DO
     END DO
     
  END IF
  !===========================================================================================================
  
  
  IF (add_yes) THEN
!pgi$ unroll = n:8
     fine(1:imax,1:jmax,1:kmax) = fine(1:imax,1:jmax,1:kmax) + work(1:imax,1:jmax,1:kmax)
  END IF
  
  
  END SUBROUTINE interpolate
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE interpolate_Helmholtz(coarse,fine,work)
  
  IMPLICIT NONE
  
  INTEGER, PARAMETER     ::  g = 1
  
  REAL   , INTENT(inout) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  REAL   , INTENT(inout) ::  fine  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL   , INTENT(out  ) ::  work  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  
  INTEGER                ::  i, di, imax, iimax
  INTEGER                ::  j, dj, jmax, jjmax
  INTEGER                ::  k, dk, kmax, kkmax
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Null-Setzen am Rand nicht notwendig!                                                      !
  !              - di /= 1 <===> N1 /= 1                                                                     !
  !              - Da imax=N1 usw., könnten auf dem feinen Gitter auch problemlos die bekannten engeren      !
  !                Intervallgrenzen S11:N11 usw. benutzt werden. Wurde bislang aus Übersichtsgründen nicht   !
  !                vollzogen.                                                                                !
  !              - Die Block-überlappenden Stirnflächen werden ebenfalls mitverarbeitet, aber eigentlich     !
  !                nicht gebraucht (erleichtert die Programmierung), so dass eigentlich eine Initialisierung !
  !                notwendig wäre. Dies wird jedoch zuvor schon in der korrespondierenden Restriktions-      !
  !                Routine erledigt, so dass dies hier nicht mehr notwendig ist.                             !
  !              - Es wird sequentiell über alle Raumrichtungen interpoliert, um keinen individuellen        !
  !                Interpolationsstencil für jeden Punkt im Raum speichern zu müssen.                        !
  !              - Durch das sequentielle Interpolieren kann der Interpolationsstencil klein und damit der   !
  !                Gesamtaufwand minimiert werden (Alternative: 8- bzw. 26-Punkt Stencil (!)). Nachteilig    !
  !                ist dabei das zusätzliche Arbeitsfeld auf dem feineren Gitterniveau (wird der Multigrid-  !
  !                Routine entliehen).                                               .                       !
  !              - Interpolationskoeffizienten werden auf dem jeweils feineren Gitter gespeichert, um nicht  !
  !                auf die entsprechenden Indizes des gröberen Gitters umrechnen zu müssen.                  !
  !----------------------------------------------------------------------------------------------------------!
  
  
  iimax = NN(1,g+1) ! TEST!!! iimax, etc. substituieren ...
  jjmax = NN(2,g+1)
  kkmax = NN(3,g+1)
  
  imax  = NN(1,g)
  jmax  = NN(2,g)
  kmax  = NN(3,g)
  
  di = (imax-1)/(iimax-1)
  dj = (jmax-1)/(jjmax-1)
  dk = (kmax-1)/(kkmax-1)
  
  
  ! ACHTUNG!!! Verwendung von parametrischen Strides di,dj,dk in Schleifen verhindert die Vektorisierung bzw.
  !            das Prefetching!
  
  ! Anmerkung: .FALSE. ist ok!
  CALL exchange_relax(g+1,ex1,ex2,ex3,0,.FALSE.,coarse)
  
  
!pgi$ unroll = n:8
  work(1:imax:di,1:jmax:dj,1:kmax:dk) = coarse(1:iimax,1:jjmax,1:kkmax)
  
  
  !===========================================================================================================
  IF (direction == 1) THEN
     
     IF (dk /= 1) THEN ! (dimens == 2) <==> (dk == 1) automatisch erfüllt!
        DO k = 2, kmax-1, dk
           DO j = 1, jmax, dj
!pgi$ unroll = n:8
              DO i = 1, imax, di
                 work(i,j,k) = cI3(1,k,g)*work(i,j,k-1) + cI3(2,k,g)*work(i,j,k+1)
              END DO
           END DO
        END DO
     END IF
     
     !--------------------------------------------------------------------------------------------------------
     IF (dj /= 1) THEN
        DO k = 1, kmax
           DO j = 2, jmax-1, dj
!pgi$ unroll = n:8
              DO i = 1, imax, di
                 work(i,j,k) = cI2(1,j,g)*work(i,j-1,k) + cI2(2,j,g)*work(i,j+1,k)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
     IF (di /= 1) THEN
        DO k = 1, kmax
           DO j = 1, jmax
!pgi$ unroll = n:8
              DO i = 1, imax-2, di
                 fine(i  ,j,k) = fine(i  ,j,k) + cIH1(1,i  )*work(i,j,k) + cIH1(2,i  )*work(i+2,j,k)
                 fine(i+1,j,k) = fine(i+1,j,k) + cIH1(1,i+1)*work(i,j,k) + cIH1(2,i+1)*work(i+2,j,k)
              END DO
           END DO
        END DO
        
        !--- Extrapolation am Rand ---
        IF (BC_1L .GT. 0) THEN
           DO k = 1, kmax
!pgi$ unroll = n:8
              DO j = 1, jmax
                 fine(0   ,j,k) = fine(0   ,j,k) + cIH1(1,0   )*work(1      ,j,k) + cIH1(2,0   )*work(1+di,j,k)
              END DO
           END DO
        END IF
        IF (BC_1U .GT. 0) THEN
           DO k = 1, kmax
!pgi$ unroll = n:8
              DO j = 1, jmax
                 fine(imax,j,k) = fine(imax,j,k) + cIH1(1,imax)*work(imax-di,j,k) + cIH1(2,imax)*work(imax,j,k)
              END DO
           END DO
        END IF
     ELSE
        IF (rank == 0) WRITE(*,*) 'ERROR! This choice of sub-grids for the Helmholtz-problem is not intended!'
        CALL MPI_FINALIZE(merror)
        STOP
     END IF
     
  END IF
  !===========================================================================================================
  IF (direction == 2) THEN
     
     IF (dk /= 1) THEN ! (dimens == 2) <==> (dk == 1) automatisch erfüllt!
        DO k = 2, kmax-1, dk
           DO j = 1, jmax, dj
!pgi$ unroll = n:8
              DO i = 1, imax, di
                 work(i,j,k) = cI3(1,k,g)*work(i,j,k-1) + cI3(2,k,g)*work(i,j,k+1)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (di /= 1) THEN
        DO k = 1, kmax
           DO j = 1, jmax, dj
!pgi$ unroll = n:8
              DO i = 2, imax-1, di
                 work(i,j,k) = cI1(1,i,g)*work(i-1,j,k) + cI1(2,i,g)*work(i+1,j,k)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (dj /= 1) THEN
        DO k = 1, kmax
           DO i = 1, imax
!pgi$ unroll = n:8
              DO j = 1, jmax-2, dj
                 fine(i,j  ,k) = fine(i,j  ,k) + cIH2(1,j  )*work(i,j,k) + cIH2(2,j  )*work(i,j+2,k)
                 fine(i,j+1,k) = fine(i,j+1,k) + cIH2(1,j+1)*work(i,j,k) + cIH2(2,j+1)*work(i,j+2,k)
              END DO
           END DO
        END DO
        
        !--- Extrapolation am Rand ---
        IF (BC_2L .GT. 0) THEN
           DO k = 1, kmax
!pgi$ unroll = n:8
              DO i = 1, imax
                 fine(i,0   ,k) = fine(i,0   ,k) + cIH2(1,0   )*work(i,1      ,k) + cIH2(2,0   )*work(i,1+dj,k)
              END DO
           END DO
        END IF
        IF (BC_2U .GT. 0) THEN
           DO k = 1, kmax
!pgi$ unroll = n:8
              DO i = 1, imax
                 fine(i,jmax,k) = fine(i,jmax,k) + cIH2(1,jmax)*work(i,jmax-dj,k) + cIH2(2,jmax)*work(i,jmax,k)
              END DO
           END DO
        END IF
     ELSE
        IF (rank == 0) WRITE(*,*) 'ERROR! This choice of sub-grids for the Helmholtz-problem is not intended!'
        CALL MPI_FINALIZE(merror)
        STOP
     END IF
     
  END IF
  !===========================================================================================================
  IF (direction == 3) THEN ! (dimens == 2) <==> (direction /= 3) automatisch erfüllt!
     
     IF (dj /= 1) THEN
        DO k = 1, kmax, dk
           DO j = 2, jmax-1, dj
!pgi$ unroll = n:8
              DO i = 1, imax, di
                 work(i,j,k) = cI2(1,j,g)*work(i,j-1,k) + cI2(2,j,g)*work(i,j+1,k)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (di /= 1) THEN
        DO k = 1, kmax, dk
           DO j = 1, jmax
!pgi$ unroll = n:8
              DO i = 2, imax-1, di
                 work(i,j,k) = cI1(1,i,g)*work(i-1,j,k) + cI1(2,i,g)*work(i+1,j,k)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (dk /= 1) THEN
        DO j = 1, jmax
           DO i = 1, imax
!pgi$ unroll = n:8
              DO k = 1, kmax-2, dk
                 fine(i,j,k  ) = fine(i,j,k  ) + cIH3(1,k  )*work(i,j,k) + cIH3(2,k  )*work(i,j,k+2)
                 fine(i,j,k+1) = fine(i,j,k+1) + cIH3(1,k+1)*work(i,j,k) + cIH3(2,k+1)*work(i,j,k+2)
              END DO
           END DO
        END DO
        
        !--- Extrapolation am Rand ---
        IF (BC_3L .GT. 0) THEN
           DO j = 1, jmax
!pgi$ unroll = n:8
              DO i = 1, imax
                 fine(i,j,0   ) = fine(i,j,0   ) + cIH3(1,0   )*work(i,j,1      ) + cIH3(2,0   )*work(i,j,1+dk)
              END DO
           END DO
        END IF
        IF (BC_3U .GT. 0) THEN
           DO j = 1, jmax
!pgi$ unroll = n:8
              DO i = 1, imax
                 fine(i,j,kmax) = fine(i,j,kmax) + cIH3(1,kmax)*work(i,j,kmax-dk) + cIH3(2,kmax)*work(i,j,kmax)
              END DO
           END DO
        END IF
     ELSE
        IF (rank == 0) WRITE(*,*) 'ERROR! This choice of sub-grids for the Helmholtz-problem is not intended!'
        CALL MPI_FINALIZE(merror)
        STOP
     END IF
     
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE interpolate_Helmholtz
  
  
  
  
  
  
  
  
  
  
  !> subroutine for solving LDEs with the iterative bi-conjugate gradient stabilized method.
  !! The subroutine is used both for the preconditioned pressure iteration (Poisson problem) as well as for the 
  !! velocity iteration (Helmholtz problems).
  !! The routine is only called with problem_type={1,2,3,4} unlike BiCGstab2 which is called only with problem_type=5.
  !! @param[in] preconditioner referring to the preconditioner for the Krylov-subspace methods (biCGstab), not the preconditioner
  !!        of the LDE matrices. 0=:no preconditioning, 1:=V-cycle MG smoothing, 2:=F-cycle MG smoothing.
  !! @note: why aren't BiCGstab and BiCGstab2 combined into a single routine?
   SUBROUTINE BiCGstab(eps,n_it_max,init_yes,SS1,SS2,SS3,NN1,NN2,NN3,bb,phi,problem_type,quiet_yes1,quiet_yes2,preconditioner)
  
  IMPLICIT NONE
  
  REAL   , INTENT(in)    ::  eps      !< termination criterion 
  INTEGER, INTENT(in)    ::  n_it_max !< max iteration steps
  LOGICAL, INTENT(in)    ::  init_yes !< whether to initialize (what?)
  
  INTEGER, INTENT(in)    ::  SS1      !< lower index bound dimension 1
  INTEGER, INTENT(in)    ::  SS2      !< lower index bound dimension 2
  INTEGER, INTENT(in)    ::  SS3      !< lower index bound dimension 3
  
  INTEGER, INTENT(in)    ::  NN1      !< upper index bound dimension 1
  INTEGER, INTENT(in)    ::  NN2      !< upper index bound dimension 2
  INTEGER, INTENT(in)    ::  NN3      !< upper index bound dimension 3
  
  INTEGER, INTENT(in)    ::  preconditioner ! see above 
  INTEGER, INTENT(in)    ::  problem_type   !< problem type to be solved.
  
  REAL   , INTENT(in)    ::  bb (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< residuum 
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< solution vector
  
  LOGICAL, INTENT(in)    ::  quiet_yes1 !< output suppression? 
  LOGICAL, INTENT(in)    ::  quiet_yes2 !< output suppression?
  
  INTEGER                ::  counter
  INTEGER                ::  i, j, k
  
  REAL                   ::  norm2   !, norm2_global
  REAL                   ::  norm_inf, norm_inf_global, norm_inf_prev
  REAL                   ::  rhr, rhr_prev, ArAr, rAr, rhAp
  REAL                   ::  scalar_global(1:2)
  REAL                   ::  alpha, beta, omega
  
  LOGICAL                ::  exit_yes
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Gitterstreckung ist ein Problem bei der Gewichtung der Residuen, was beim Helmholtz-      !
  !                Problem aufgrund der Diagonaldominanz weit weniger ausgeprägt sein sollte als beim reinen !
  !                Poisson-Problem. Da der Multigrid-Glätter (Jacobi, Gauss-Seidel) davon generell unbeein-  !
  !                flusst sind, kann dieses Problem hier mehr oder weniger vernachlässigt werden. Alternativ !
  !                könnten die Stencils analog zu "weight" auch für das gesamte Feld gespeichert werden, was !
  !                allerdings einen sehr grossen Speicheraufwand und lange Speicherzugriffszeiten erwarten   !
  !                lassen, so dass hier eine evtl. geringere Konvergenzrate mit gutem Gewissen in Kauf       !
  !                genommen wird.                                                                            !
  !              - Es wird hier generell rh = rr angenommen, da somit immer <rh,rr> /= 0 gilt!               !
  !              - dAXPY ist nicht wirklich schneller als die direkt programmierten Operationen.             !
  !              - rhr muss nicht initialisiert werden!                                                      !
  !              - Um pp nur bei Bedarf berechnen zu müssen, wird dazu eine Fallunterscheidung nach erst     !
  !                hinter "CALL status_iteration" eingeführt.                                                !
  !              - z2 könnte zu Ungunsten der Geschwindigkeit eingespaart werden.                            !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! - Geschwindigkeit ist trotz Prefetching / Vektorisierung stark durch die Speicherzugriffszeit limitiert
  ! - Das wird z.B. deutlich bei Single- vs. Dualcorebetrieb
  
  
  !===========================================================================================================
  IF (.NOT. init_yes) THEN
     IF (problem_type == 1) CALL product_Helmholtz(phi,rr)
     IF (problem_type == 2) CALL product_div_grad (phi,rr)
     !-----------------------------------------------------------------------------------------------------------
     ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
     !-----------------------------------------------------------------------------------------------------------
     IF (problem_type == 4) CALL product_div_grad_transp(phi,rr)
     
     rr(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3) - rr(SS1:NN1,SS2:NN2,SS3:NN3)
     rh(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3)
  END IF
  !===========================================================================================================
  
  !===========================================================================================================
  !=== Residuum ==============================================================================================
  !===========================================================================================================
  IF (init_yes) THEN
     CALL get_norms(SS1,SS2,SS3,NN1,NN2,NN3,bb,problem_type,.TRUE.,.TRUE. ,normInf=norm_inf,normTwo=norm2)
  ELSE
     CALL get_norms(SS1,SS2,SS3,NN1,NN2,NN3,rr,problem_type,.TRUE.,.FALSE.,normInf=norm_inf)
  END IF
  norm_inf_prev = norm_inf
  !===========================================================================================================
  
  counter = 0
  
  
  ITERATE: DO
     
     !========================================================================================================
     !=== Überprüfen des Konvergenzkriteriums ================================================================
     !========================================================================================================
     CALL status_iteration(eps,norm_inf,counter,n_it_max,exit_yes,quiet_yes1,quiet_yes2)
     
     ! soll Konvergenz sicherstellen (oder offsetprec runtersetzen ...).
     !    norm_inf == 0. << eps ist ein Sonderfall, bei dem es keinen Sinn macht, eine weitere Iteration
     !    zu rechnen. 
     IF (problem_type == 1 .AND.                           counter .LT. 1 .AND. norm_inf /= 0.) exit_yes = .FALSE.
     IF (problem_type == 2 .AND. number_poisson == 1 .AND. counter .LT. 1 .AND. norm_inf /= 0.) exit_yes = .FALSE.
     IF (problem_type == 2 .AND. number_poisson == 2 .AND. counter .LT. 0 .AND. norm_inf /= 0.) exit_yes = .FALSE.
     
     IF (exit_yes .AND. counter == 0 .AND. init_yes) phi(SS1:NN1,SS2:NN2,SS3:NN3) = 0.
     IF (exit_yes) EXIT ITERATE
     !========================================================================================================
     
     
     !========================================================================================================
     !=== nächster Durchlauf =================================================================================
     !========================================================================================================
     counter = counter + 1
     !========================================================================================================
     
     
     !========================================================================================================
     rhr_prev = rhr
     IF (init_yes) THEN
        IF (counter == 1) THEN
           rhr = norm2
        ELSE
           CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,rr,bb,rhr)
        END IF
     ELSE
        CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,rr,rh,rhr)
     END IF
     !========================================================================================================
     IF (ABS(rhr) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rhr =', rhr
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'rhr =', rhr
        
        rh(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3) ! Neuer Referenzvektor ...
        CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,rr,rh,rhr)
     END IF
     !========================================================================================================
     IF (counter .GE. 2) THEN
        IF (omega == 0.) THEN
           IF (rank == 0) WRITE(* ,'(a,E13.5)') 'omega =', omega
           IF (rank == 0) WRITE(10,'(a,E13.5)') 'omega =', omega
           IF (problem_type == 2 .OR. problem_type == 4) THEN
              EXIT ITERATE
           ELSE
              CALL MPI_FINALIZE(merror)
              STOP
           END IF
        END IF
        IF (rhr_prev == 0.) THEN
           IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rhr_prev =', rhr_prev
           IF (rank == 0) WRITE(10,'(a,E13.5)') 'rhr_prev =', rhr_prev
           IF (problem_type == 2 .OR. problem_type == 4) THEN
              EXIT ITERATE
           ELSE
              CALL MPI_FINALIZE(merror)
              STOP
           END IF
        END IF
        beta = (alpha/omega)*(rhr/rhr_prev)
        omega = -beta*omega
        pp(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3) + beta*pp(SS1:NN1,SS2:NN2,SS3:NN3) + omega*Ap(SS1:NN1,SS2:NN2,SS3:NN3)
     ELSE
        IF (init_yes .AND. counter == 1) THEN
           pp(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3)
        ELSE
           pp(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3)
        END IF
     END IF
     !========================================================================================================
     IF (preconditioner == 0) THEN
        IF (problem_type == 1) CALL product_Helmholtz(pp,Ap)
        IF (problem_type == 2) CALL product_div_grad (pp,Ap)
        !-----------------------------------------------------------------------------------------------------------
        ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
        !-----------------------------------------------------------------------------------------------------------
        IF (problem_type == 4) CALL product_div_grad_transp(pp,Ap)
     ELSE IF (preconditioner == 1 .OR. preconditioner == 2) THEN
        IF (preconditioner == 1) CALL multigridV(.TRUE.,1,pp,z1,problem_type)
        IF (preconditioner == 2) CALL multigridF(.TRUE.,1,pp,z1,problem_type)
        IF (problem_type == 1) CALL product_Helmholtz(z1,Ap)
        IF (problem_type == 2) CALL product_div_grad (z1,Ap)
        !-----------------------------------------------------------------------------------------------------------
        ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
        !-----------------------------------------------------------------------------------------------------------
        IF (problem_type == 4) CALL product_div_grad_transp(z1,Ap)
     ELSE
        IF (rank == 0) WRITE(*,'(a)') 'ERROR! Specify valid preconditioner!'
        CALL MPI_FINALIZE(merror)
        STOP
     END IF
     !========================================================================================================
     IF (init_yes) THEN
        CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ap,bb,rhAp)
     ELSE
        CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ap,rh,rhAp)
     END IF
     !========================================================================================================
     IF (ABS(rhAp) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rhAp =', rhAp
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'rhAp =', rhAp
        IF (problem_type == 2 .OR. problem_type == 4) THEN
           EXIT ITERATE
        ELSE
           IF (rhr /= 0.) THEN
              CALL MPI_FINALIZE(merror)
              STOP
           END IF
        END IF
        alpha = 0.
     ELSE
        alpha = rhr / rhAp
     END IF
     !========================================================================================================
     IF (init_yes .AND. counter == 1) THEN
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 rr(i,j,k) = bb(i,j,k) - alpha*Ap(i,j,k)
              END DO
           END DO
        END DO
     ELSE
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 rr(i,j,k) = rr(i,j,k) - alpha*Ap(i,j,k)
              END DO
           END DO
        END DO
     END IF
     !========================================================================================================
     IF (preconditioner == 0) THEN
        IF (problem_type == 1) CALL product_Helmholtz(rr,Ar)
        IF (problem_type == 2) CALL product_div_grad (rr,Ar)
        !-----------------------------------------------------------------------------------------------------------
        ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
        !-----------------------------------------------------------------------------------------------------------
        IF (problem_type == 4) CALL product_div_grad_transp(rr,Ar)
     ELSE IF (preconditioner == 1 .OR. preconditioner == 2) THEN
        IF (preconditioner == 1) CALL multigridV(.TRUE.,1,rr,z2,problem_type)
        IF (preconditioner == 2) CALL multigridF(.TRUE.,1,rr,z2,problem_type)
        IF (problem_type == 1) CALL product_Helmholtz(z2,Ar)
        IF (problem_type == 2) CALL product_div_grad (z2,Ar)
        !-----------------------------------------------------------------------------------------------------------
        ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
        !-----------------------------------------------------------------------------------------------------------
        IF (problem_type == 4) CALL product_div_grad_transp(z2,Ar)
     END IF
     !========================================================================================================
     rAr  = 0.
     ArAr = 0.
     DO k = SS3, NN3
        DO j = SS2, NN2
!pgi$ unroll = n:8
           DO i = SS1, NN1
              rAr  = rAr  + rr(i,j,k)*Ar(i,j,k)
              ArAr = ArAr + Ar(i,j,k)*Ar(i,j,k)
           END DO
        END DO
     END DO
     
     CALL MPI_ALLREDUCE((/rAr,ArAr/),scalar_global,2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     rAr  = scalar_global(1)
     ArAr = scalar_global(2)
     
     ! ACHTUNG!!! Zu teuer im Vergleich zur obigen Variante:
     !CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ar,rr,rAr)
     !CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ar,Ar,ArAr)
     !========================================================================================================
     IF (ABS(rAr) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rAr =', rAr
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'rAr =', rAr
     END IF
     IF (ABS(ArAr) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'ArAr =', ArAr
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'ArAr =', ArAr
        IF (problem_type == 2 .OR. problem_type == 4) THEN
           EXIT ITERATE
        ELSE
           IF (rAr /= 0.) THEN
              CALL MPI_FINALIZE(merror)
              STOP
           END IF
        END IF
        omega = 0.
     ELSE
        omega = rAr / ArAr
     END IF
     !========================================================================================================
     norm_inf = 0.
     !--------------------------------------------------------------------------------------------------------
     IF (counter == 1 .AND. init_yes) THEN
        IF (preconditioner == 0) THEN
           IF (problem_type == 2 .AND. weighting_yes) THEN
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = alpha*pp(i,j,k) + omega*rr(i,j,k)
                       rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           ELSE
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = alpha*pp(i,j,k) + omega*rr(i,j,k)
                       rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           END IF
        !-----------------------------------------------------------------------------------------------------
        ELSE
           IF (problem_type == 2 .AND. weighting_yes) THEN
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = alpha*z1(i,j,k) + omega*z2(i,j,k)
                       rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           ELSE
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = alpha*z1(i,j,k) + omega*z2(i,j,k)
                       rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           END IF
        END IF
     !--------------------------------------------------------------------------------------------------------
     ELSE
        IF (preconditioner == 0) THEN
           IF (problem_type == 2 .AND. weighting_yes) THEN
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = phi(i,j,k) + omega*rr(i,j,k) + alpha*pp(i,j,k)
                       rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           ELSE
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = phi(i,j,k) + omega*rr(i,j,k) + alpha*pp(i,j,k)
                       rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           END IF
        !-----------------------------------------------------------------------------------------------------
        ELSE
           IF (problem_type == 2 .AND. weighting_yes) THEN
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = phi(i,j,k) + omega*z2(i,j,k) + alpha*z1(i,j,k)
                       rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           ELSE
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = phi(i,j,k) + omega*z2(i,j,k) + alpha*z1(i,j,k)
                       rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           END IF
        END IF
     END IF
     !--------------------------------------------------------------------------------------------------------
     CALL MPI_ALLREDUCE(norm_inf,norm_inf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
     norm_inf = norm_inf_global
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Konvergenzstatistik ================================================================================
     !========================================================================================================
     IF (problem_type == 1) ratioH(substep,direction     ) = ratioH(substep,direction     ) + LOG10(norm_inf/norm_inf_prev)
     IF (problem_type == 2) ratioP(substep,number_poisson) = ratioP(substep,number_poisson) + LOG10(norm_inf/norm_inf_prev)
     
     IF (problem_type == 1) countH(substep,direction     ) = countH(substep,direction     ) + 1
     IF (problem_type == 2) countP(substep,number_poisson) = countP(substep,number_poisson) + 1
     
     norm_inf_prev = norm_inf
     !========================================================================================================
     
  END DO ITERATE
  
  
  END SUBROUTINE BiCGstab
  
  
  
  
  
  
  
  
  
  ! TEST!!! gleiches Spiel; g und SNB uebergeben ....
  ! ACHTUNG!!! Routine sollte nur bei g .GE. 2 aufgerufen werden, weil sonst viele grosse Felder allociert werden ...
  ! TEST!!! weight in Residuum?
  ! N1, N2, N3, g werden übergeben, um z.B. psi auch auf gröberen Gittern berechnen zu können ...
  SUBROUTINE BiCGstab2(eps,n_it_max,init_yes,N1,N2,N3,g,SS1,SS2,SS3,NN1,NN2,NN3,bb,phi,problem_type,quiet_yes1,quiet_yes2,preconditioner)
  
  IMPLICIT NONE
  
  REAL   , INTENT(in)    ::  eps
  INTEGER, INTENT(in)    ::  n_it_max
  LOGICAL, INTENT(in)    ::  init_yes
  
  INTEGER, INTENT(in)    ::  N1
  INTEGER, INTENT(in)    ::  N2
  INTEGER, INTENT(in)    ::  N3
  
  INTEGER, INTENT(in)    ::  g
  
  INTEGER, INTENT(in)    ::  SS1
  INTEGER, INTENT(in)    ::  SS2
  INTEGER, INTENT(in)    ::  SS3
  
  INTEGER, INTENT(in)    ::  NN1
  INTEGER, INTENT(in)    ::  NN2
  INTEGER, INTENT(in)    ::  NN3
  
  INTEGER, INTENT(in)    ::  preconditioner
  INTEGER, INTENT(in)    ::  problem_type
  
  REAL   , INTENT(in)    ::  bb (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  REAL                   ::  pp(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) ! TEST!!! ! FELD!!!
  REAL                   ::  Ap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  rr(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  rh(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  Ar(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  z1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  z2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  LOGICAL, INTENT(in)    ::  quiet_yes1
  LOGICAL, INTENT(in)    ::  quiet_yes2
  
  INTEGER                ::  counter
  INTEGER                ::  i, j, k
  
  REAL                   ::  norm2   !, norm2_global
  REAL                   ::  norm_inf, norm_inf_global, norm_inf_prev
  REAL                   ::  rhr, rhr_prev, ArAr, rAr, rhAp
  REAL                   ::  scalar_global(1:2)
  REAL                   ::  alpha, beta, omega
  
  LOGICAL                ::  exit_yes
  
  
  !*****************************************************************************************
  !SNB(:,:,g) = (/SS1,NN1/), etc. ! TEST!!!
  !*****************************************************************************************
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Gitterstreckung ist ein Problem bei der Gewichtung der Residuen, was beim Helmholtz-      !
  !                Problem aufgrund der Diagonaldominanz weit weniger ausgeprägt sein sollte als beim reinen !
  !                Poisson-Problem. Da der Multigrid-Glätter (Jacobi, Gauss-Seidel) davon generell unbeein-  !
  !                flusst sind, kann dieses Problem hier mehr oder weniger vernachlässigt werden. Alternativ !
  !                könnten die Stencils analog zu "weight" auch für das gesamte Feld gespeichert werden, was !
  !                allerdings einen sehr grossen Speicheraufwand und lange Speicherzugriffszeiten erwarten   !
  !                lassen, so dass hier eine evtl. geringere Konvergenzrate mit gutem Gewissen in Kauf       !
  !                genommen wird.                                                                            !
  !              - Es wird hier generell rh = rr angenommen, da somit immer <rh,rr> /= 0 gilt!               !
  !              - dAXPY ist nicht wirklich schneller als die direkt programmierten Operationen.             !
  !              - rhr muss nicht initialisiert werden!                                                      !
  !              - Um pp nur bei Bedarf berechnen zu müssen, wird dazu eine Fallunterscheidung nach erst     !
  !                hinter "CALL status_iteration" eingeführt.                                                !
  !              - z2 könnte zu Ungunsten der Geschwindigkeit eingespaart werden.                            !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! - Geschwindigkeit ist trotz Prefetching / Vektorisierung stark durch die Speicherzugriffszeit limitiert
  ! - Das wird z.B. deutlich bei Single- vs. Dualcorebetrieb
  
  
  !===========================================================================================================
  IF (.NOT. init_yes) THEN
     !IF (problem_type == 1) CALL product_Helmholtz(phi,rr)
     !IF (problem_type == 2) CALL product_div_grad (phi,rr)
     !IF (problem_type == 3) CALL product_Helmholtz_conc(phi,rr)
     !IF (problem_type == 4) CALL product_div_grad_transp(phi,rr)
     IF (problem_type == 5) CALL product_div_grad_relax(g,phi,rr)
     
     rr(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3) - rr(SS1:NN1,SS2:NN2,SS3:NN3)
     rh(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3)
  END IF
  !===========================================================================================================
  
  !===========================================================================================================
  !=== Residuum ==============================================================================================
  !===========================================================================================================
  IF (init_yes) THEN
     CALL get_norms2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,bb,problem_type,.TRUE.,.TRUE. ,normInf=norm_inf,normTwo=norm2)
  ELSE
     CALL get_norms2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,rr,problem_type,.TRUE.,.FALSE.,normInf=norm_inf)
  END IF
  norm_inf_prev = norm_inf
  !===========================================================================================================
  
  counter = 0
  
  
  ITERATE: DO
     
     !========================================================================================================
     !=== Überprüfen des Konvergenzkriteriums ================================================================
     !========================================================================================================
     CALL status_iteration(eps,norm_inf,counter,n_it_max,exit_yes,quiet_yes1,quiet_yes2)
     
     ! soll Konvergenz sicherstellen (oder offsetprec runtersetzen ...).
     !    norm_inf == 0. << eps ist ein Sonderfall, bei dem es keinen Sinn macht, eine weitere Iteration
     !    zu rechnen. 
     IF (problem_type == 1 .AND.                           counter .LT. 1 .AND. norm_inf /= 0.) exit_yes = .FALSE.
     IF (problem_type == 2 .AND. number_poisson == 1 .AND. counter .LT. 1 .AND. norm_inf /= 0.) exit_yes = .FALSE.
     IF (problem_type == 2 .AND. number_poisson == 2 .AND. counter .LT. 0 .AND. norm_inf /= 0.) exit_yes = .FALSE.
     
     IF (exit_yes .AND. counter == 0 .AND. init_yes) phi(SS1:NN1,SS2:NN2,SS3:NN3) = 0.
     IF (exit_yes) EXIT ITERATE
     !========================================================================================================
     
     
     !========================================================================================================
     !=== nächster Durchlauf =================================================================================
     !========================================================================================================
     counter = counter + 1
     !========================================================================================================
     
     
     !========================================================================================================
     rhr_prev = rhr
     IF (init_yes) THEN
        IF (counter == 1) THEN
           rhr = norm2
        ELSE
           CALL product_scalar2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,rr,bb,rhr)
        END IF
     ELSE
        CALL product_scalar2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,rr,rh,rhr)
     END IF
     !========================================================================================================
     IF (ABS(rhr) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rhr =', rhr
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'rhr =', rhr
        
        rh(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3) ! Neuer Referenzvektor ...
        CALL product_scalar2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,rr,rh,rhr)
     END IF
     !========================================================================================================
     IF (counter .GE. 2) THEN
        IF (omega == 0.) THEN
           IF (rank == 0) WRITE(* ,'(a,E13.5)') 'omega =', omega
           IF (rank == 0) WRITE(10,'(a,E13.5)') 'omega =', omega
           IF (problem_type == 2 .OR. problem_type == 4 .OR. problem_type == 5) THEN
              EXIT ITERATE
           ELSE
              CALL MPI_FINALIZE(merror)
              STOP
           END IF
        END IF
        IF (rhr_prev == 0.) THEN
           IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rhr_prev =', rhr_prev
           IF (rank == 0) WRITE(10,'(a,E13.5)') 'rhr_prev =', rhr_prev
           IF (problem_type == 2 .OR. problem_type == 4 .OR. problem_type == 5) THEN
              EXIT ITERATE
           ELSE
              CALL MPI_FINALIZE(merror)
              STOP
           END IF
        END IF
        beta = (alpha/omega)*(rhr/rhr_prev)
        omega = -beta*omega
        pp(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3) + beta*pp(SS1:NN1,SS2:NN2,SS3:NN3) + omega*Ap(SS1:NN1,SS2:NN2,SS3:NN3)
     ELSE
        IF (init_yes .AND. counter == 1) THEN
           pp(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3)
        ELSE
           pp(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3)
        END IF
     END IF
     !========================================================================================================
     IF (preconditioner == 0) THEN
        !IF (problem_type == 1) CALL product_Helmholtz(pp,Ap)
        !IF (problem_type == 2) CALL product_div_grad (pp,Ap)
        !IF (problem_type == 3) CALL product_Helmholtz_conc(pp,Ap)
        !IF (problem_type == 4) CALL product_div_grad_transp(pp,Ap)
        IF (problem_type == 5) CALL product_div_grad_relax(g,pp,Ap)
     ELSE IF (preconditioner == 1 .OR. preconditioner == 2) THEN
        IF (preconditioner == 1) CALL multigridV(.TRUE.,g,pp,z1,problem_type)
        IF (preconditioner == 2) CALL multigridF(.TRUE.,g,pp,z1,problem_type)
        !IF (problem_type == 1) CALL product_Helmholtz(z1,Ap)
        !IF (problem_type == 2) CALL product_div_grad (z1,Ap)
        !IF (problem_type == 3) CALL product_Helmholtz_conc(z1,Ap)
        !IF (problem_type == 4) CALL product_div_grad_transp(z1,Ap)
        IF (problem_type == 5) CALL product_div_grad_relax(g,z1,Ap)
     ELSE
        IF (rank == 0) WRITE(*,'(a)') 'ERROR! Specify valid preconditioner!'
        CALL MPI_FINALIZE(merror)
        STOP
     END IF
     !========================================================================================================
     IF (init_yes) THEN
        CALL product_scalar2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,Ap,bb,rhAp)
     ELSE
        CALL product_scalar2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,Ap,rh,rhAp)
     END IF
     !========================================================================================================
     IF (ABS(rhAp) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rhAp =', rhAp
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'rhAp =', rhAp
        IF (problem_type == 2 .OR. problem_type == 4 .OR. problem_type == 5) THEN
           EXIT ITERATE
        ELSE
           IF (rhr /= 0.) THEN
              CALL MPI_FINALIZE(merror)
              STOP
           END IF
        END IF
        alpha = 0.
     ELSE
        alpha = rhr / rhAp
     END IF
     !========================================================================================================
     IF (init_yes .AND. counter == 1) THEN
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 rr(i,j,k) = bb(i,j,k) - alpha*Ap(i,j,k)
              END DO
           END DO
        END DO
     ELSE
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 rr(i,j,k) = rr(i,j,k) - alpha*Ap(i,j,k)
              END DO
           END DO
        END DO
     END IF
     !========================================================================================================
     IF (preconditioner == 0) THEN
        !IF (problem_type == 1) CALL product_Helmholtz(rr,Ar)
        !IF (problem_type == 2) CALL product_div_grad (rr,Ar)
        !IF (problem_type == 3) CALL product_Helmholtz_conc(rr,Ar)
        !IF (problem_type == 4) CALL product_div_grad_transp(rr,Ar)
        IF (problem_type == 5) CALL product_div_grad_relax(g,rr,Ar)
     ELSE IF (preconditioner == 1 .OR. preconditioner == 2) THEN
        IF (preconditioner == 1) CALL multigridV(.TRUE.,g,rr,z2,problem_type)
        IF (preconditioner == 2) CALL multigridF(.TRUE.,g,rr,z2,problem_type)
        !IF (problem_type == 1) CALL product_Helmholtz(z2,Ar)
        !IF (problem_type == 2) CALL product_div_grad (z2,Ar)
        !IF (problem_type == 3) CALL product_Helmholtz_conc(z2,Ar)
        !IF (problem_type == 4) CALL product_div_grad_transp(z2,Ar)
        IF (problem_type == 5) CALL product_div_grad_relax(g,z2,Ar)
     END IF
     !========================================================================================================
     rAr  = 0.
     ArAr = 0.
     DO k = SS3, NN3
        DO j = SS2, NN2
!pgi$ unroll = n:8
           DO i = SS1, NN1
              rAr  = rAr  + rr(i,j,k)*Ar(i,j,k)
              ArAr = ArAr + Ar(i,j,k)*Ar(i,j,k)
           END DO
        END DO
     END DO
     
     CALL MPI_ALLREDUCE((/rAr,ArAr/),scalar_global,2,MPI_REAL8,MPI_SUM,comm1(g),merror)
     rAr  = scalar_global(1)
     ArAr = scalar_global(2)
     
     ! ACHTUNG!!! Zu teuer im Vergleich zur obigen Variante:
     !CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ar,rr,rAr)
     !CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ar,Ar,ArAr)
     !========================================================================================================
     IF (ABS(rAr) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rAr =', rAr
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'rAr =', rAr
     END IF
     IF (ABS(ArAr) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'ArAr =', ArAr
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'ArAr =', ArAr
        IF (problem_type == 2 .OR. problem_type == 4 .OR. problem_type == 5) THEN
           EXIT ITERATE
        ELSE
           IF (rAr /= 0.) THEN
              CALL MPI_FINALIZE(merror)
              STOP
           END IF
        END IF
        omega = 0.
     ELSE
        omega = rAr / ArAr
     END IF
     !========================================================================================================
     norm_inf = 0.
     !--------------------------------------------------------------------------------------------------------
     IF (counter == 1 .AND. init_yes) THEN
        IF (preconditioner == 0) THEN
           IF (problem_type == 2 .AND. weighting_yes) THEN
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = alpha*pp(i,j,k) + omega*rr(i,j,k)
                       rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           ELSE
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = alpha*pp(i,j,k) + omega*rr(i,j,k)
                       rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           END IF
        !-----------------------------------------------------------------------------------------------------
        ELSE
           IF (problem_type == 2 .AND. weighting_yes) THEN
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = alpha*z1(i,j,k) + omega*z2(i,j,k)
                       rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           ELSE
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = alpha*z1(i,j,k) + omega*z2(i,j,k)
                       rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           END IF
        END IF
     !--------------------------------------------------------------------------------------------------------
     ELSE
        IF (preconditioner == 0) THEN
           IF (problem_type == 2 .AND. weighting_yes) THEN
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = phi(i,j,k) + omega*rr(i,j,k) + alpha*pp(i,j,k)
                       rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           ELSE
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = phi(i,j,k) + omega*rr(i,j,k) + alpha*pp(i,j,k)
                       rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           END IF
        !-----------------------------------------------------------------------------------------------------
        ELSE
           IF (problem_type == 2 .AND. weighting_yes) THEN
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = phi(i,j,k) + omega*z2(i,j,k) + alpha*z1(i,j,k)
                       rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           ELSE
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = phi(i,j,k) + omega*z2(i,j,k) + alpha*z1(i,j,k)
                       rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           END IF
        END IF
     END IF
     !--------------------------------------------------------------------------------------------------------
     CALL MPI_ALLREDUCE(norm_inf,norm_inf_global,1,MPI_REAL8,MPI_MAX,comm1(g),merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
     norm_inf = norm_inf_global
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Konvergenzstatistik ================================================================================
     !========================================================================================================
     IF (problem_type == 1) ratioH(substep,direction     ) = ratioH(substep,direction     ) + LOG10(norm_inf/norm_inf_prev)
     IF (problem_type == 2) ratioP(substep,number_poisson) = ratioP(substep,number_poisson) + LOG10(norm_inf/norm_inf_prev)
     
     IF (problem_type == 1) countH(substep,direction     ) = countH(substep,direction     ) + 1
     IF (problem_type == 2) countP(substep,number_poisson) = countP(substep,number_poisson) + 1
     
     norm_inf_prev = norm_inf
     !========================================================================================================
     
  END DO ITERATE
  
  
  END SUBROUTINE BiCGstab2
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE Richardson(eps,n_it_max,init_yes,SS1,SS2,SS3,NN1,NN2,NN3,bb,phi,problem_type,quiet_yes1,quiet_yes2,preconditioner)
  
  IMPLICIT NONE
  
  REAL   , INTENT(in)    ::  eps
  INTEGER, INTENT(in)    ::  n_it_max
  LOGICAL, INTENT(in)    ::  init_yes
  
  INTEGER, INTENT(in)    ::  SS1
  INTEGER, INTENT(in)    ::  SS2
  INTEGER, INTENT(in)    ::  SS3
  
  INTEGER, INTENT(in)    ::  NN1
  INTEGER, INTENT(in)    ::  NN2
  INTEGER, INTENT(in)    ::  NN3
  
  INTEGER, INTENT(in)    ::  preconditioner
  INTEGER, INTENT(in)    ::  problem_type
  
  REAL   , INTENT(in)    ::  bb (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  LOGICAL, INTENT(in)    ::  quiet_yes1
  LOGICAL, INTENT(in)    ::  quiet_yes2
  
  INTEGER                ::  counter
  REAL                   ::  norm_inf, norm_inf0
  LOGICAL                ::  exit_yes
  
  
  
  !------------------------------------------------------------
  IF (init_yes) THEN
!pgi$ unroll = n:8
     rr(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3)
     phi = 0.
  ELSE
     IF (problem_type == 1) CALL product_Helmholtz(phi,rr)
     IF (problem_type == 2) CALL product_div_grad (phi,rr)
     !-----------------------------------------------------------------------------------------------------------
     ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
     !-----------------------------------------------------------------------------------------------------------
     IF (problem_type == 4) CALL product_div_grad_transp(phi,rr)
!pgi$ unroll = n:8
     rr(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3) - rr(SS1:NN1,SS2:NN2,SS3:NN3)
  END IF
  
  !--- grösstes Residuum --------------------------------------
  CALL get_norms(SS1,SS2,SS3,NN1,NN2,NN3,rr,problem_type,.TRUE.,.FALSE.,normInf=norm_inf)
  norm_inf0 = norm_inf
  
  !------------------------------------------------------------
  counter = 0
  !------------------------------------------------------------
  Ap(SS1:NN1,SS2:NN2,SS3:NN3) = 0.
  pp(SS1:NN1,SS2:NN2,SS3:NN3) = 0.
  !------------------------------------------------------------
  
  ITERATE: DO
     
     !--- Überprüfen des Konvergenzkriteriums -----------------
     CALL status_iteration(eps,norm_inf,counter,n_it_max,exit_yes,quiet_yes1,quiet_yes2)
     
     IF (exit_yes) EXIT ITERATE
     
     
     !--- nächster Restart ------------------------------------
     counter = counter + 1
     
     
     IF (preconditioner == 0) THEN
        pp = 0.25*rr ! TEST!!! nicht sehr schön ...
        !pp = 1.0*rr
     ELSE IF (preconditioner == 1 .OR. preconditioner == 2) THEN
        IF (preconditioner == 1) CALL multigridV(.TRUE.,1,rr,pp,problem_type)
        IF (preconditioner == 2) CALL multigridF(.TRUE.,1,rr,pp,problem_type)
     ELSE
        IF (rank == 0) WRITE(*,'(a)') 'ERROR! Specify valid preconditioner!'
        CALL MPI_FINALIZE(merror)
        STOP
     END IF
     
     
     IF (1 == 1) THEN
        
        IF (problem_type == 1) CALL product_Helmholtz      (pp,Ap)
        IF (problem_type == 2) CALL product_div_grad       (pp,Ap)
        !-----------------------------------------------------------------------------------------------------------
        ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
        !-----------------------------------------------------------------------------------------------------------
        IF (problem_type == 4) CALL product_div_grad_transp(pp,Ap)
        
        rr (SS1:NN1,SS2:NN2,SS3:NN3) = rr (SS1:NN1,SS2:NN2,SS3:NN3) - Ap(SS1:NN1,SS2:NN2,SS3:NN3)
        phi(SS1:NN1,SS2:NN2,SS3:NN3) = phi(SS1:NN1,SS2:NN2,SS3:NN3) + pp(SS1:NN1,SS2:NN2,SS3:NN3)
        
     ELSE
        
        phi(SS1:NN1,SS2:NN2,SS3:NN3) = phi(SS1:NN1,SS2:NN2,SS3:NN3) + pp(SS1:NN1,SS2:NN2,SS3:NN3)
        
        IF (problem_type == 1) CALL product_Helmholtz      (phi,Ap)
        IF (problem_type == 2) CALL product_div_grad       (phi,Ap)
	!-----------------------------------------------------------------------------------------------------------
	! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
	!-----------------------------------------------------------------------------------------------------------
        IF (problem_type == 4) CALL product_div_grad_transp(phi,Ap)
        
        rr(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3) - Ap(SS1:NN1,SS2:NN2,SS3:NN3)
     END IF
     
     
     !--- grösstes Residuum --------------------------------------
     CALL get_norms(SS1,SS2,SS3,NN1,NN2,NN3,rr,problem_type,.TRUE.,.FALSE.,normInf=norm_inf)
     
     
     !--- Konvergenzstatistik ---------------------------------
     IF (problem_type == 1) ratioH(substep,direction     ) = ratioH(substep,direction     ) + LOG10(norm_inf/norm_inf0)
     IF (problem_type == 2) ratioP(substep,number_poisson) = ratioP(substep,number_poisson) + LOG10(norm_inf/norm_inf0)
     
     IF (problem_type == 1) countH(substep,direction     ) = countH(substep,direction     ) + 1
     IF (problem_type == 2) countP(substep,number_poisson) = countP(substep,number_poisson) + 1
     
     norm_inf0 = norm_inf
     
  END DO ITERATE
  
  
  END SUBROUTINE Richardson
  
  
  
  
  
  
  
  
  
  
  !> computed the element-wise product of two fields (arrays) 
  SUBROUTINE product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,phi1,phi2,scalar)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)    ::  SS1 !< lower bound dimension one
  INTEGER, INTENT(in)    ::  SS2 !< lower bound dimension two
  INTEGER, INTENT(in)    ::  SS3 !< lower bound dimension three
  
  INTEGER, INTENT(in)    ::  NN1 !< upper bound dimension one
  INTEGER, INTENT(in)    ::  NN2 !< upper bound dimension two
  INTEGER, INTENT(in)    ::  NN3 !< upper bound dimension three
  
  REAL   , INTENT(in)    ::  phi1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< array 1
  REAL   , INTENT(in)    ::  phi2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< array 2
  
  REAL   , INTENT(out)   ::  scalar !< element-wise product
  REAL                   ::  scalar_global
  INTEGER                ::  i, j, k
  
  
  scalar = 0.
  
  DO k = SS3, NN3
     DO j = SS2, NN2
!pgi$ unroll = n:8
        DO i = SS1, NN1
           scalar = scalar + phi1(i,j,k)*phi2(i,j,k)
        END DO
     END DO
  END DO
  
  CALL MPI_ALLREDUCE(scalar,scalar_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  scalar = scalar_global
  
  
  END SUBROUTINE product_scalar
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! N1, N2, N3 werden ebenfalls uebergeben ...
  SUBROUTINE product_scalar2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,phi1,phi2,scalar)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)    ::  g
  
  INTEGER, INTENT(in)    ::  N1
  INTEGER, INTENT(in)    ::  N2
  INTEGER, INTENT(in)    ::  N3
  
  INTEGER, INTENT(in)    ::  SS1
  INTEGER, INTENT(in)    ::  SS2
  INTEGER, INTENT(in)    ::  SS3
  
  INTEGER, INTENT(in)    ::  NN1
  INTEGER, INTENT(in)    ::  NN2
  INTEGER, INTENT(in)    ::  NN3
  
  REAL   , INTENT(in)    ::  phi1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(in)    ::  phi2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  REAL   , INTENT(out)   ::  scalar
  REAL                   ::  scalar_global
  INTEGER                ::  i, j, k
  
  
  scalar = 0.
  
  DO k = SS3, NN3
     DO j = SS2, NN2
!pgi$ unroll = n:8
        DO i = SS1, NN1
           scalar = scalar + phi1(i,j,k)*phi2(i,j,k)
        END DO
     END DO
  END DO
  
  CALL MPI_ALLREDUCE(scalar,scalar_global,1,MPI_REAL8,MPI_SUM,comm1(g),merror)
  scalar = scalar_global
  
  
  END SUBROUTINE product_scalar2
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE multadd1(SS1,SS2,SS3,NN1,NN2,NN3,mult,vec,phi)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)    ::  SS1
  INTEGER, INTENT(in)    ::  SS2
  INTEGER, INTENT(in)    ::  SS3
  
  INTEGER, INTENT(in)    ::  NN1
  INTEGER, INTENT(in)    ::  NN2
  INTEGER, INTENT(in)    ::  NN3
  
  REAL   , INTENT(in)    ::  mult
  
  REAL   , INTENT(inout) ::  vec(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(in)    ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, j, k
  
  
  DO k = SS3, NN3
     DO j = SS2, NN2
!pgi$ unroll = n:8
        DO i = SS1, NN1
           vec(i,j,k) = vec(i,j,k) + mult*phi(i,j,k)
        END DO
     END DO
  END DO
  
  
  END SUBROUTINE multadd1
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE multadd2(SS1,SS2,SS3,NN1,NN2,NN3,mult1,mult2,vec,phi1,phi2,init_yes)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)    ::  SS1
  INTEGER, INTENT(in)    ::  SS2
  INTEGER, INTENT(in)    ::  SS3
  
  INTEGER, INTENT(in)    ::  NN1
  INTEGER, INTENT(in)    ::  NN2
  INTEGER, INTENT(in)    ::  NN3
  
  REAL   , INTENT(in)    ::  mult1
  REAL   , INTENT(in)    ::  mult2
  
  REAL   , INTENT(inout) ::  vec (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(in)    ::  phi1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(in)    ::  phi2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  LOGICAL, INTENT(in)    ::  init_yes
  
  INTEGER                ::  i, j, k
  
  
  IF (init_yes) THEN
     
     DO k = SS3, NN3
        DO j = SS2, NN2
!pgi$ unroll = n:8
           DO i = SS1, NN1
              vec(i,j,k) = mult1*phi1(i,j,k) + mult2*phi2(i,j,k)
           END DO
        END DO
     END DO
     
  ELSE
     
     DO k = SS3, NN3
        DO j = SS2, NN2
!pgi$ unroll = n:8
           DO i = SS1, NN1
              vec(i,j,k) = vec(i,j,k) + mult1*phi1(i,j,k) + mult2*phi2(i,j,k)
           END DO
        END DO
     END DO
     
  END IF
  
  
  END SUBROUTINE multadd2
  
  
  
  
  
  
  
  
  
  
  !> subroutine that compute infinity and/or two norms.  
  SUBROUTINE get_norms(SS1,SS2,SS3,NN1,NN2,NN3,phi,problem_type,inf_yes,two_yes,normInf,normTwo)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)    ::  SS1 !< lower bound first dimension
  INTEGER, INTENT(in)    ::  SS2 !< lower bound second dimension
  INTEGER, INTENT(in)    ::  SS3 !< lower bound third dimension
  
  INTEGER, INTENT(in)    ::  NN1 !< upper bound first dimension
  INTEGER, INTENT(in)    ::  NN2 !< upper bound second dimension
  INTEGER, INTENT(in)    ::  NN3 !< upper bound third dimension
  
  REAL   , INTENT(in)    ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< field to compute the norm of
  
  INTEGER, INTENT(in)    ::  problem_type !< type of problem
  LOGICAL, INTENT(in)    ::  inf_yes      !< whether to compute infinity norm
  LOGICAL, INTENT(in)    ::  two_yes      !< whether to compute two-norm
  
  REAL   , OPTIONAL, INTENT(out) ::  normInf
  REAL   , OPTIONAL, INTENT(out) ::  normTwo
  
  REAL                   ::  normInf_global, normTwo_global
  INTEGER                ::  i, j, k
  
  
  IF (inf_yes .AND. two_yes) THEN
     
     normInf = 0.
     normTwo = 0.
     
     IF (problem_type == 2 .AND. weighting_yes) THEN
        
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)*weight(i,j,k)),normInf)
                 normTwo = normTwo + phi(i,j,k)**2
              END DO
           END DO
        END DO
        
     ELSE
        
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)),normInf)
                 normTwo = normTwo + phi(i,j,k)**2
              END DO
           END DO
        END DO
        
     END IF
     
     ! Lassen sich wegen MPI_SUM / MPI_MAX nicht zusammenlegen:
     CALL MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror)
     CALL MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     normInf = normInf_global
     normTwo = normTwo_global
     
  ELSE IF (inf_yes) THEN
     
     normInf = 0.
     
     IF (problem_type == 2 .AND. weighting_yes) THEN
        
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)*weight(i,j,k)),normInf)
              END DO
           END DO
        END DO
        
     ELSE
        
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)),normInf)
              END DO
           END DO
        END DO
        
     END IF
     
     CALL MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
     normInf = normInf_global
     
  ELSE IF (two_yes) THEN
     
     normTwo = 0.
     
     DO k = SS3, NN3
        DO j = SS2, NN2
!pgi$ unroll = n:8
           DO i = SS1, NN1
              normTwo = normTwo + phi(i,j,k)**2
           END DO
        END DO
     END DO
     
     CALL MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     normTwo = normTwo_global
     
  END IF
  
  
  END SUBROUTINE get_norms
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! N1, N2, N3 werden ebenfalls uebergeben ...
  SUBROUTINE get_norms2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,phi,problem_type,inf_yes,two_yes,normInf,normTwo)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)    ::  g
  
  INTEGER, INTENT(in)    ::  N1
  INTEGER, INTENT(in)    ::  N2
  INTEGER, INTENT(in)    ::  N3
  
  INTEGER, INTENT(in)    ::  SS1
  INTEGER, INTENT(in)    ::  SS2
  INTEGER, INTENT(in)    ::  SS3
  
  INTEGER, INTENT(in)    ::  NN1
  INTEGER, INTENT(in)    ::  NN2
  INTEGER, INTENT(in)    ::  NN3
  
  REAL   , INTENT(in)    ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER, INTENT(in)    ::  problem_type
  LOGICAL, INTENT(in)    ::  inf_yes
  LOGICAL, INTENT(in)    ::  two_yes
  
  REAL   , OPTIONAL, INTENT(out) ::  normInf
  REAL   , OPTIONAL, INTENT(out) ::  normTwo
  
  REAL                   ::  normInf_global, normTwo_global
  INTEGER                ::  i, j, k
  
  
  IF (inf_yes .AND. two_yes) THEN
     
     normInf = 0.
     normTwo = 0.
     
     IF (problem_type == 2 .AND. weighting_yes) THEN
        
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)*weight(i,j,k)),normInf)
                 normTwo = normTwo + phi(i,j,k)**2
              END DO
           END DO
        END DO
        
     ELSE
        
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)),normInf)
                 normTwo = normTwo + phi(i,j,k)**2
              END DO
           END DO
        END DO
        
     END IF
     
     ! Lassen sich wegen MPI_SUM / MPI_MAX nicht zusammenlegen:
     CALL MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,comm1(g),merror)
     CALL MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,comm1(g),merror)
     normInf = normInf_global
     normTwo = normTwo_global
     
  ELSE IF (inf_yes) THEN
     
     normInf = 0.
     
     IF (problem_type == 2 .AND. weighting_yes) THEN
        
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)*weight(i,j,k)),normInf)
              END DO
           END DO
        END DO
        
     ELSE
        
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)),normInf)
              END DO
           END DO
        END DO
        
     END IF
     
     CALL MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,comm1(g),merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
     normInf = normInf_global
     
  ELSE IF (two_yes) THEN
     
     normTwo = 0.
     
     DO k = SS3, NN3
        DO j = SS2, NN2
!pgi$ unroll = n:8
           DO i = SS1, NN1
              normTwo = normTwo + phi(i,j,k)**2
           END DO
        END DO
     END DO
     
     CALL MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,comm1(g),merror)
     normTwo = normTwo_global
     
  END IF
  
  
  END SUBROUTINE get_norms2
  
  
  
  
  
  
  
  
  
  
  !> checks whether iteration criterion was met or whether maximum numbers of iterations have been reached.
  SUBROUTINE status_iteration(eps,norm,counter,N_restarts,exit_yes,quiet_yes1,quiet_yes2)
  
  IMPLICIT NONE
  
  REAL   , INTENT(in)    ::  eps       !< convergence criterion
  REAL   , INTENT(in)    ::  norm      !< (two, infinity) norm
  INTEGER, INTENT(in)    ::  counter   !< iteration counter (current)
  INTEGER, INTENT(in)    ::  N_restarts!< max iteration count
  LOGICAL, INTENT(out)   ::  exit_yes  !< whether to stop the iteration
  LOGICAL, INTENT(in)    ::  quiet_yes1!< whether to suppress output during iteration
  LOGICAL, INTENT(in)    ::  quiet_yes2!< whether to suppress output for termination
  
  
  exit_yes = .FALSE.
  
  IF (norm .LT. eps) THEN
     IF (rank == 0 .AND. write_stout_yes .AND. .NOT. quiet_yes2) WRITE(* ,'(a,i5,a,E24.17,a)') '  Iteration',counter,'; ||res|| =',norm,'  (Termination criterion satisfied)'
     IF (rank == 0 .AND. log_iteration_yes                     ) WRITE(10,'(a,i5,a,E24.17,a)') '  Iteration',counter,'; ||res|| =',norm,'  (Termination criterion satisfied)'
     exit_yes = .TRUE.
  END IF
  
  IF ((.NOT. exit_yes) .AND. counter == N_restarts) THEN
     IF (rank == 0 .AND. write_stout_yes .AND. .NOT. quiet_yes2) WRITE(* ,'(a,i5,a,E24.17,a)') '  Iteration',counter,'; ||res|| =',norm,'  WARNING! Too many iterations!'
     IF (rank == 0 .AND. log_iteration_yes                     ) WRITE(10,'(a,i5,a,E24.17,a)') '  Iteration',counter,'; ||res|| =',norm,'  WARNING! Too many iterations!'
     exit_yes = .TRUE.
  END IF
  
  IF (.NOT. exit_yes) THEN
     IF (rank == 0 .AND. write_stout_yes .AND. .NOT. quiet_yes1) WRITE(* ,'(a,i5,a,E24.17  )') '  Iteration',counter,'; ||res|| =',norm
     IF (rank == 0 .AND. log_iteration_yes                     ) WRITE(10,'(a,i5,a,E24.17  )') '  Iteration',counter,'; ||res|| =',norm
  END IF
  
  
  END SUBROUTINE status_iteration
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE handle_corner_rhs(g,phi) ! TEST!!! Ersetzen mit einheitlicher Routine?? ! TEST!!! validieren ...
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  g
  REAL   , INTENT(inout) ::  phi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Siehe Subroutine "handle_corner_Lap"!                                                     !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (BC(1,1,g) .GT. 0 .AND. BC(1,2,g) .GT. 0) phi(1      ,1      ,1:NN(3,g)) = 0.
  IF (BC(1,1,g) .GT. 0 .AND. BC(2,2,g) .GT. 0) phi(1      ,NN(2,g),1:NN(3,g)) = 0.
  IF (BC(2,1,g) .GT. 0 .AND. BC(1,2,g) .GT. 0) phi(NN(1,g),1      ,1:NN(3,g)) = 0.
  IF (BC(2,1,g) .GT. 0 .AND. BC(2,2,g) .GT. 0) phi(NN(1,g),NN(2,g),1:NN(3,g)) = 0.
  
  IF (BC(1,1,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) phi(1      ,1:NN(2,g),1      ) = 0.
  IF (BC(1,1,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) phi(1      ,1:NN(2,g),NN(3,g)) = 0.
  IF (BC(2,1,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) phi(NN(1,g),1:NN(2,g),1      ) = 0.
  IF (BC(2,1,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) phi(NN(1,g),1:NN(2,g),NN(3,g)) = 0.
  
  IF (BC(1,2,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) phi(1:NN(1,g),1      ,1      ) = 0.
  IF (BC(1,2,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) phi(1:NN(1,g),1      ,NN(3,g)) = 0.
  IF (BC(2,2,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) phi(1:NN(1,g),NN(2,g),1      ) = 0.
  IF (BC(2,2,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) phi(1:NN(1,g),NN(2,g),NN(3,g)) = 0.
  
  
  END SUBROUTINE handle_corner_rhs
  
  
  
  
END MODULE mod_solvers
