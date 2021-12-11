!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becske, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!*************************************************************************************************************

!> module containing subroutines for computing the rhs. It uses modules mod_dims, mod_vars, mod_diff, 
!! mod_particles and mod_les
MODULE mod_rhs
  
  
  USE mod_dims
  USE mod_vars
  USE mod_diff
  USE mpi  
  
  PRIVATE
  
  PUBLIC rhs_vel
  
  !INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> computes the rhs of the velocity
  SUBROUTINE rhs_vel
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k, m
  REAL                   ::  mult
  REAL                   ::  flux, flux_global
  
  
  ! TEST!!! generell muesste eigentlich nur der "innere" Bereich (d.h. ohne Raender) von rhs" belegt werden!
  
  !===========================================================================================================
  !=== rhs = rhs + Helmholtz(vel) ============================================================================
  !===========================================================================================================
  IF (timeint_mode == 1 .OR. thetaL == 1.) THEN
     DO m = 1, dimens
        IF (m == 1) rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
        IF (m == 2) rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
        IF (m == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
        
     END DO
  ELSE
     multL = -(1.-thetaL)*(aRK(substep)+bRK(substep))*dtime/Re
     DO m = 1, dimens
        CALL Helmholtz(m,.FALSE.,vel(b1L,b2L,b3L,m),rhs(b1L,b2L,b3L,m))
     END DO
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== rhs = rhs - nl_old ====================================================================================
  !===========================================================================================================
  IF (substep /= 1) THEN
     mult = dtime*bRK(substep)
     
     DO m = 1, dimens
        IF (m == 1) rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) - mult*nl(S11B:N11B,S21B:N21B,S31B:N31B,1)
        IF (m == 2) rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) - mult*nl(S12B:N12B,S22B:N22B,S32B:N32B,2)
        IF (m == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) - mult*nl(S13B:N13B,S23B:N23B,S33B:N33B,3)
     END DO
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== rhs = rhs - nl ========================================================================================
  !===========================================================================================================
  !--- advective terms ---------------------------------------------------------------------------------------
  IF (Stokes_yes) THEN
                      nl(S11B:N11B,S21B:N21B,S31B:N31B,1) = 0.
                      nl(S12B:N12B,S22B:N22B,S32B:N32B,2) = 0.
     IF (dimens == 3) nl(S13B:N13B,S23B:N23B,S33B:N33B,3) = 0.
  ELSE
     CALL nonlinear(.FALSE.)
  END IF
    
  !--- viscose terms -----------------------------------------------------------------------------------------
  IF (timeint_mode == 1 .AND. (.NOT. Euler_yes)) THEN
     multL = 1./Re
     CALL Helmholtz_explicit(.FALSE.)
  END IF
  
  !--- forcing -----------------------------------------------------------------------------------------------
  CALL forcing_vel

  !--- rhs = rhs - nl ----------------------------------------------------------------------------------------
  mult = dtime*aRK(substep)
  
  DO m = 1, dimens
     IF (m == 1) rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) - mult*nl(S11B:N11B,S21B:N21B,S31B:N31B,1)
     IF (m == 2) rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) - mult*nl(S12B:N12B,S22B:N22B,S32B:N32B,2)
     IF (m == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) - mult*nl(S13B:N13B,S23B:N23B,S33B:N33B,3)
  END DO
  !===========================================================================================================
  
  
  
  ! TEST!!! BC-Abfrage unklar bzw. inkonsistent ...
  !===========================================================================================================
  !=== Randbedingungen =======================================================================================
  !===========================================================================================================
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Momentan sind nur Dirichlet-Randbedingungen für Geschwindigkeiten wählbar.                !
  !              - "vel" darf NICHT überschrieben werden, da sonst die Konzentrationen nicht korrekt berech- !
  !                net werden können.                                                                        !
  !----------------------------------------------------------------------------------------------------------!
  
  !----------------------------------------------------------------------------------------------------------!
  ! ACHTUNG: Es ist sehr genau auf die "Stetigkeit" der RB über Ecken/Kanten zu achten, bzw. die Reihenfolge ! 
  !          der Operationen und die Intervallgrenzen müssen UNBEDINGT sorgfaltig gewählt werden (Grund:     !
  !          Divergenz wird auch auf dem Rand gebildet). Siehe dazu auch die Reihenfolge der RB in           !
  !          product_Helmholtz!                                                                              !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== Randbedingungen (Zeitintegration) =====================================================================
  !===========================================================================================================
  IF (substep /= 1) THEN
     !--------------------------------------------------------------------------------------------------------
     !--- rhs = rhs - nl_old ---------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     mult = dtime*bRK(substep)
     !--------------------------------------------------------------------------------------------------------
     IF (outlet(1,1,1) .AND. (BC_1L == 1 .OR. BC_1L == 2)) bc11(S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1) - mult*nlbc11(S21B:N21B,S31B:N31B,1)
     IF (outlet(1,2,1) .AND. (BC_1U == 1 .OR. BC_1U == 2)) bc11(S21B:N21B,S31B:N31B,2) = bc11(S21B:N21B,S31B:N31B,2) - mult*nlbc11(S21B:N21B,S31B:N31B,2)
     
     IF (outlet(1,1,2) .AND. (BC_1L == 1 .OR. BC_1L == 3)) bc21(S22B:N22B,S32B:N32B,1) = bc21(S22B:N22B,S32B:N32B,1) - mult*nlbc21(S22B:N22B,S32B:N32B,1)
     IF (outlet(1,2,2) .AND. (BC_1U == 1 .OR. BC_1U == 3)) bc21(S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2) - mult*nlbc21(S22B:N22B,S32B:N32B,2)
     
     IF (outlet(1,1,3) .AND. (BC_1L == 1 .OR. BC_1L == 3)) bc31(S23B:N23B,S33B:N33B,1) = bc31(S23B:N23B,S33B:N33B,1) - mult*nlbc31(S23B:N23B,S33B:N33B,1)
     IF (outlet(1,2,3) .AND. (BC_1U == 1 .OR. BC_1U == 3)) bc31(S23B:N23B,S33B:N33B,2) = bc31(S23B:N23B,S33B:N33B,2) - mult*nlbc31(S23B:N23B,S33B:N33B,2)
     !--------------------------------------------------------------------------------------------------------
     IF (outlet(2,1,1) .AND. (BC_2L == 1 .OR. BC_2L == 3)) bc12(S11B:N11B,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1) - mult*nlbc12(S11B:N11B,S31B:N31B,1)
     IF (outlet(2,2,1) .AND. (BC_2U == 1 .OR. BC_2U == 3)) bc12(S11B:N11B,S31B:N31B,2) = bc12(S11B:N11B,S31B:N31B,2) - mult*nlbc12(S11B:N11B,S31B:N31B,2)
     
     IF (outlet(2,1,2) .AND. (BC_2L == 1 .OR. BC_2L == 2)) bc22(S12B:N12B,S32B:N32B,1) = bc22(S12B:N12B,S32B:N32B,1) - mult*nlbc22(S12B:N12B,S32B:N32B,1)
     IF (outlet(2,2,2) .AND. (BC_2U == 1 .OR. BC_2U == 2)) bc22(S12B:N12B,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2) - mult*nlbc22(S12B:N12B,S32B:N32B,2)
     
     IF (outlet(2,1,3) .AND. (BC_2L == 1 .OR. BC_2L == 3)) bc32(S13B:N13B,S33B:N33B,1) = bc32(S13B:N13B,S33B:N33B,1) - mult*nlbc32(S13B:N13B,S33B:N33B,1)
     IF (outlet(2,2,3) .AND. (BC_2U == 1 .OR. BC_2U == 3)) bc32(S13B:N13B,S33B:N33B,2) = bc32(S13B:N13B,S33B:N33B,2) - mult*nlbc32(S13B:N13B,S33B:N33B,2)
     !--------------------------------------------------------------------------------------------------------
     IF (outlet(3,1,1) .AND. (BC_3L == 1 .OR. BC_3L == 3)) bc13(S11B:N11B,S21B:N21B,1) = bc13(S11B:N11B,S21B:N21B,1) - mult*nlbc13(S11B:N11B,S21B:N21B,1)
     IF (outlet(3,2,1) .AND. (BC_3U == 1 .OR. BC_3U == 3)) bc13(S11B:N11B,S21B:N21B,2) = bc13(S11B:N11B,S21B:N21B,2) - mult*nlbc13(S11B:N11B,S21B:N21B,2)
     
     IF (outlet(3,1,2) .AND. (BC_3L == 1 .OR. BC_3L == 3)) bc23(S12B:N12B,S22B:N22B,1) = bc23(S12B:N12B,S22B:N22B,1) - mult*nlbc23(S12B:N12B,S22B:N22B,1)
     IF (outlet(3,2,2) .AND. (BC_3U == 1 .OR. BC_3U == 3)) bc23(S12B:N12B,S22B:N22B,2) = bc23(S12B:N12B,S22B:N22B,2) - mult*nlbc23(S12B:N12B,S22B:N22B,2)
     
     IF (outlet(3,1,3) .AND. (BC_3L == 1 .OR. BC_3L == 2)) bc33(S13B:N13B,S23B:N23B,1) = bc33(S13B:N13B,S23B:N23B,1) - mult*nlbc33(S13B:N13B,S23B:N23B,1)
     IF (outlet(3,2,3) .AND. (BC_3U == 1 .OR. BC_3U == 2)) bc33(S13B:N13B,S23B:N23B,2) = bc33(S13B:N13B,S23B:N23B,2) - mult*nlbc33(S13B:N13B,S23B:N23B,2)
     !--------------------------------------------------------------------------------------------------------
  END IF
  
  !--- Ausfluss-RB -------------------------------------------------------------------------------------------
  CALL outflow_bc
  
  !--- andere RB (instationaer + Zeitintegration) ------------------------------------------------------------
  CALL boundary_vel_tint
  
  !-----------------------------------------------------------------------------------------------------------
  !--- rhs = rhs - nl ----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  mult = dtime*aRK(substep)
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(1,1,1) .AND. (BC_1L == 1 .OR. BC_1L == 2)) bc11(S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1) - mult*nlbc11(S21B:N21B,S31B:N31B,1)
  IF (outlet(1,2,1) .AND. (BC_1U == 1 .OR. BC_1U == 2)) bc11(S21B:N21B,S31B:N31B,2) = bc11(S21B:N21B,S31B:N31B,2) - mult*nlbc11(S21B:N21B,S31B:N31B,2)
  
  IF (outlet(1,1,2) .AND. (BC_1L == 1 .OR. BC_1L == 3)) bc21(S22B:N22B,S32B:N32B,1) = bc21(S22B:N22B,S32B:N32B,1) - mult*nlbc21(S22B:N22B,S32B:N32B,1)
  IF (outlet(1,2,2) .AND. (BC_1U == 1 .OR. BC_1U == 3)) bc21(S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2) - mult*nlbc21(S22B:N22B,S32B:N32B,2)
  
  IF (outlet(1,1,3) .AND. (BC_1L == 1 .OR. BC_1L == 3)) bc31(S23B:N23B,S33B:N33B,1) = bc31(S23B:N23B,S33B:N33B,1) - mult*nlbc31(S23B:N23B,S33B:N33B,1)
  IF (outlet(1,2,3) .AND. (BC_1U == 1 .OR. BC_1U == 3)) bc31(S23B:N23B,S33B:N33B,2) = bc31(S23B:N23B,S33B:N33B,2) - mult*nlbc31(S23B:N23B,S33B:N33B,2)
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(2,1,1) .AND. (BC_2L == 1 .OR. BC_2L == 3)) bc12(S11B:N11B,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1) - mult*nlbc12(S11B:N11B,S31B:N31B,1)
  IF (outlet(2,2,1) .AND. (BC_2U == 1 .OR. BC_2U == 3)) bc12(S11B:N11B,S31B:N31B,2) = bc12(S11B:N11B,S31B:N31B,2) - mult*nlbc12(S11B:N11B,S31B:N31B,2)
  
  IF (outlet(2,1,2) .AND. (BC_2L == 1 .OR. BC_2L == 2)) bc22(S12B:N12B,S32B:N32B,1) = bc22(S12B:N12B,S32B:N32B,1) - mult*nlbc22(S12B:N12B,S32B:N32B,1)
  IF (outlet(2,2,2) .AND. (BC_2U == 1 .OR. BC_2U == 2)) bc22(S12B:N12B,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2) - mult*nlbc22(S12B:N12B,S32B:N32B,2)
  
  IF (outlet(2,1,3) .AND. (BC_2L == 1 .OR. BC_2L == 3)) bc32(S13B:N13B,S33B:N33B,1) = bc32(S13B:N13B,S33B:N33B,1) - mult*nlbc32(S13B:N13B,S33B:N33B,1)
  IF (outlet(2,2,3) .AND. (BC_2U == 1 .OR. BC_2U == 3)) bc32(S13B:N13B,S33B:N33B,2) = bc32(S13B:N13B,S33B:N33B,2) - mult*nlbc32(S13B:N13B,S33B:N33B,2)
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(3,1,1) .AND. (BC_3L == 1 .OR. BC_3L == 3)) bc13(S11B:N11B,S21B:N21B,1) = bc13(S11B:N11B,S21B:N21B,1) - mult*nlbc13(S11B:N11B,S21B:N21B,1)
  IF (outlet(3,2,1) .AND. (BC_3U == 1 .OR. BC_3U == 3)) bc13(S11B:N11B,S21B:N21B,2) = bc13(S11B:N11B,S21B:N21B,2) - mult*nlbc13(S11B:N11B,S21B:N21B,2)
  
  IF (outlet(3,1,2) .AND. (BC_3L == 1 .OR. BC_3L == 3)) bc23(S12B:N12B,S22B:N22B,1) = bc23(S12B:N12B,S22B:N22B,1) - mult*nlbc23(S12B:N12B,S22B:N22B,1)
  IF (outlet(3,2,2) .AND. (BC_3U == 1 .OR. BC_3U == 3)) bc23(S12B:N12B,S22B:N22B,2) = bc23(S12B:N12B,S22B:N22B,2) - mult*nlbc23(S12B:N12B,S22B:N22B,2)
  
  IF (outlet(3,1,3) .AND. (BC_3L == 1 .OR. BC_3L == 2)) bc33(S13B:N13B,S23B:N23B,1) = bc33(S13B:N13B,S23B:N23B,1) - mult*nlbc33(S13B:N13B,S23B:N23B,1)
  IF (outlet(3,2,3) .AND. (BC_3U == 1 .OR. BC_3U == 2)) bc33(S13B:N13B,S23B:N23B,2) = bc33(S13B:N13B,S23B:N23B,2) - mult*nlbc33(S13B:N13B,S23B:N23B,2)
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Randbedingungen (instantan) ===========================================================================
  !===========================================================================================================
  CALL boundary_vel_stat
  !===========================================================================================================
  
  
  
  
  !===========================================================================================================
  !=== RB einsetzen ==========================================================================================
  !===========================================================================================================
  ! Anmerkung: Reihenfolge bestimmt Konvention (Linien, Kanten)!
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_2L .GT. 0) rhs(S11B:N11B,1        ,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1)
  IF (BC_2U .GT. 0) rhs(S11B:N11B,N2       ,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,2)
  
  IF (BC_3L .GT. 0) rhs(S11B:N11B,S21B:N21B,1        ,1) = bc13(S11B:N11B,S21B:N21B,1)
  IF (BC_3U .GT. 0) rhs(S11B:N11B,S21B:N21B,N3       ,1) = bc13(S11B:N11B,S21B:N21B,2)
  
  IF (BC_1L .GT. 0) rhs(0        ,S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1)
  IF (BC_1U .GT. 0) rhs(N1       ,S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,2)
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_1L .GT. 0) rhs(1        ,S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,1)
  IF (BC_1U .GT. 0) rhs(N1       ,S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2)
  
  IF (BC_3L .GT. 0) rhs(S12B:N12B,S22B:N22B,1        ,2) = bc23(S12B:N12B,S22B:N22B,1)
  IF (BC_3U .GT. 0) rhs(S12B:N12B,S22B:N22B,N3       ,2) = bc23(S12B:N12B,S22B:N22B,2)
  
  IF (BC_2L .GT. 0) rhs(S12B:N12B,0        ,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,1)
  IF (BC_2U .GT. 0) rhs(S12B:N12B,N2       ,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2)
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_1L .GT. 0) rhs(1        ,S23B:N23B,S33B:N33B,3) = bc31(S23B:N23B,S33B:N33B,1)
  IF (BC_1U .GT. 0) rhs(N1       ,S23B:N23B,S33B:N33B,3) = bc31(S23B:N23B,S33B:N33B,2)
  
  IF (BC_2L .GT. 0) rhs(S13B:N13B,1        ,S33B:N33B,3) = bc32(S13B:N13B,S33B:N33B,1)
  IF (BC_2U .GT. 0) rhs(S13B:N13B,N2       ,S33B:N33B,3) = bc32(S13B:N13B,S33B:N33B,2)
  
  IF (BC_3L .GT. 0) rhs(S13B:N13B,S23B:N23B,0        ,3) = bc33(S13B:N13B,S23B:N23B,1)
  IF (BC_3U .GT. 0) rhs(S13B:N13B,S23B:N23B,N3       ,3) = bc33(S13B:N13B,S23B:N23B,2)
  !===========================================================================================================
  
  
  
  
  !===========================================================================================================
  !=== Fluss-Korrektur =======================================================================================
  !===========================================================================================================
  IF (nullspace_yes) THEN
     flux = 0.
     
     m = 1
     DO k = S31B, N31B
        DO j = S21B, N21B
!pgi$ unroll = n:8
           DO i = S11B, N11B
              flux = flux + psi_vel(i,j,k,m)*rhs(i,j,k,m)
           END DO
        END DO
     END DO
     
     m = 2
     DO k = S32B, N32B
        DO j = S22B, N22B
!pgi$ unroll = n:8
           DO i = S12B, N12B
              flux = flux + psi_vel(i,j,k,m)*rhs(i,j,k,m)
           END DO
        END DO
     END DO
     
     IF (dimens == 3) THEN
     m = 3
     DO k = S33B, N33B
        DO j = S23B, N23B
!pgi$ unroll = n:8
           DO i = S13B, N13B
              flux = flux + psi_vel(i,j,k,m)*rhs(i,j,k,m)
           END DO
        END DO
     END DO
     END IF
     
     
     CALL MPI_ALLREDUCE(flux,flux_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     flux = flux_global
     
     IF (rank == 0) WRITE(10,'(a,E25.17)') 'flux =', flux
     
     
     ! d/dt ist unstetig in der Zeit ==> kein Versuch unternommen, in Zeitintegration einzubauen ...
     IF (nullspace_ortho_yes) THEN
        
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L .GT. 0) bc12(S11B:N11B,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1) - flux*psi_vel(S11B:N11B,1        ,S31B:N31B,1)
        IF (BC_2U .GT. 0) bc12(S11B:N11B,S31B:N31B,2) = bc12(S11B:N11B,S31B:N31B,2) - flux*psi_vel(S11B:N11B,N2       ,S31B:N31B,1)
        
        IF (BC_3L .GT. 0) bc13(S11B:N11B,S21B:N21B,1) = bc13(S11B:N11B,S21B:N21B,1) - flux*psi_vel(S11B:N11B,S21B:N21B,1        ,1)
        IF (BC_3U .GT. 0) bc13(S11B:N11B,S21B:N21B,2) = bc13(S11B:N11B,S21B:N21B,2) - flux*psi_vel(S11B:N11B,S21B:N21B,N3       ,1)
        
        IF (BC_1L .GT. 0) bc11(S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1) - flux*psi_vel(0        ,S21B:N21B,S31B:N31B,1)
        IF (BC_1U .GT. 0) bc11(S21B:N21B,S31B:N31B,2) = bc11(S21B:N21B,S31B:N31B,2) - flux*psi_vel(N1       ,S21B:N21B,S31B:N31B,1)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L .GT. 0) bc21(S22B:N22B,S32B:N32B,1) = bc21(S22B:N22B,S32B:N32B,1) - flux*psi_vel(1        ,S22B:N22B,S32B:N32B,2)
        IF (BC_1U .GT. 0) bc21(S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2) - flux*psi_vel(N1       ,S22B:N22B,S32B:N32B,2)
        
        IF (BC_3L .GT. 0) bc23(S12B:N12B,S22B:N22B,1) = bc23(S12B:N12B,S22B:N22B,1) - flux*psi_vel(S12B:N12B,S22B:N22B,1        ,2)
        IF (BC_3U .GT. 0) bc23(S12B:N12B,S22B:N22B,2) = bc23(S12B:N12B,S22B:N22B,2) - flux*psi_vel(S12B:N12B,S22B:N22B,N3       ,2)
        
        IF (BC_2L .GT. 0) bc22(S12B:N12B,S32B:N32B,1) = bc22(S12B:N12B,S32B:N32B,1) - flux*psi_vel(S12B:N12B,0        ,S32B:N32B,2)
        IF (BC_2U .GT. 0) bc22(S12B:N12B,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2) - flux*psi_vel(S12B:N12B,N2       ,S32B:N32B,2)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L .GT. 0) bc31(S23B:N23B,S33B:N33B,1) = bc31(S23B:N23B,S33B:N33B,1) - flux*psi_vel(1        ,S23B:N23B,S33B:N33B,3)
        IF (BC_1U .GT. 0) bc31(S23B:N23B,S33B:N33B,2) = bc31(S23B:N23B,S33B:N33B,2) - flux*psi_vel(N1       ,S23B:N23B,S33B:N33B,3)
        
        IF (BC_2L .GT. 0) bc32(S13B:N13B,S33B:N33B,1) = bc32(S13B:N13B,S33B:N33B,1) - flux*psi_vel(S13B:N13B,1        ,S33B:N33B,3)
        IF (BC_2U .GT. 0) bc32(S13B:N13B,S33B:N33B,2) = bc32(S13B:N13B,S33B:N33B,2) - flux*psi_vel(S13B:N13B,N2       ,S33B:N33B,3)
        
        IF (BC_3L .GT. 0) bc33(S13B:N13B,S23B:N23B,1) = bc33(S13B:N13B,S23B:N23B,1) - flux*psi_vel(S13B:N13B,S23B:N23B,0        ,3)
        IF (BC_3U .GT. 0) bc33(S13B:N13B,S23B:N23B,2) = bc33(S13B:N13B,S23B:N23B,2) - flux*psi_vel(S13B:N13B,S23B:N23B,N3       ,3)
        !-----------------------------------------------------------------------------------------------------
        
                         rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) - flux*psi_vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
                         rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) - flux*psi_vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
        IF (dimens == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) - flux*psi_vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
        
     ELSE
        
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L .GT. 0) bc12(S11B:N11B,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1) - flux*th12(S11B:N11B,S31B:N31B,1)
        IF (BC_2U .GT. 0) bc12(S11B:N11B,S31B:N31B,2) = bc12(S11B:N11B,S31B:N31B,2) - flux*th12(S11B:N11B,S31B:N31B,2)
        
        IF (BC_3L .GT. 0) bc13(S11B:N11B,S21B:N21B,1) = bc13(S11B:N11B,S21B:N21B,1) - flux*th13(S11B:N11B,S21B:N21B,1)
        IF (BC_3U .GT. 0) bc13(S11B:N11B,S21B:N21B,2) = bc13(S11B:N11B,S21B:N21B,2) - flux*th13(S11B:N11B,S21B:N21B,2)
        
        IF (BC_1L .GT. 0) bc11(S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1) - flux*th11(S21B:N21B,S31B:N31B,1)
        IF (BC_1U .GT. 0) bc11(S21B:N21B,S31B:N31B,2) = bc11(S21B:N21B,S31B:N31B,2) - flux*th11(S21B:N21B,S31B:N31B,2)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L .GT. 0) bc21(S22B:N22B,S32B:N32B,1) = bc21(S22B:N22B,S32B:N32B,1) - flux*th21(S22B:N22B,S32B:N32B,1)
        IF (BC_1U .GT. 0) bc21(S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2) - flux*th21(S22B:N22B,S32B:N32B,2)
        
        IF (BC_3L .GT. 0) bc23(S12B:N12B,S22B:N22B,1) = bc23(S12B:N12B,S22B:N22B,1) - flux*th23(S12B:N12B,S22B:N22B,1)
        IF (BC_3U .GT. 0) bc23(S12B:N12B,S22B:N22B,2) = bc23(S12B:N12B,S22B:N22B,2) - flux*th23(S12B:N12B,S22B:N22B,2)
        
        IF (BC_2L .GT. 0) bc22(S12B:N12B,S32B:N32B,1) = bc22(S12B:N12B,S32B:N32B,1) - flux*th22(S12B:N12B,S32B:N32B,1)
        IF (BC_2U .GT. 0) bc22(S12B:N12B,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2) - flux*th22(S12B:N12B,S32B:N32B,2)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L .GT. 0) bc31(S23B:N23B,S33B:N33B,1) = bc31(S23B:N23B,S33B:N33B,1) - flux*th31(S23B:N23B,S33B:N33B,1)
        IF (BC_1U .GT. 0) bc31(S23B:N23B,S33B:N33B,2) = bc31(S23B:N23B,S33B:N33B,2) - flux*th31(S23B:N23B,S33B:N33B,2)
        
        IF (BC_2L .GT. 0) bc32(S13B:N13B,S33B:N33B,1) = bc32(S13B:N13B,S33B:N33B,1) - flux*th32(S13B:N13B,S33B:N33B,1)
        IF (BC_2U .GT. 0) bc32(S13B:N13B,S33B:N33B,2) = bc32(S13B:N13B,S33B:N33B,2) - flux*th32(S13B:N13B,S33B:N33B,2)
        
        IF (BC_3L .GT. 0) bc33(S13B:N13B,S23B:N23B,1) = bc33(S13B:N13B,S23B:N23B,1) - flux*th33(S13B:N13B,S23B:N23B,1)
        IF (BC_3U .GT. 0) bc33(S13B:N13B,S23B:N23B,2) = bc33(S13B:N13B,S23B:N23B,2) - flux*th33(S13B:N13B,S23B:N23B,2)
        !-----------------------------------------------------------------------------------------------------
        
        
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L .GT. 0) rhs(S11B:N11B,1        ,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1)
        IF (BC_2U .GT. 0) rhs(S11B:N11B,N2       ,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,2)
        
        IF (BC_3L .GT. 0) rhs(S11B:N11B,S21B:N21B,1        ,1) = bc13(S11B:N11B,S21B:N21B,1)
        IF (BC_3U .GT. 0) rhs(S11B:N11B,S21B:N21B,N3       ,1) = bc13(S11B:N11B,S21B:N21B,2)
        
        IF (BC_1L .GT. 0) rhs(0        ,S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1)
        IF (BC_1U .GT. 0) rhs(N1       ,S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,2)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L .GT. 0) rhs(1        ,S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,1)
        IF (BC_1U .GT. 0) rhs(N1       ,S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2)
        
        IF (BC_3L .GT. 0) rhs(S12B:N12B,S22B:N22B,1        ,2) = bc23(S12B:N12B,S22B:N22B,1)
        IF (BC_3U .GT. 0) rhs(S12B:N12B,S22B:N22B,N3       ,2) = bc23(S12B:N12B,S22B:N22B,2)
        
        IF (BC_2L .GT. 0) rhs(S12B:N12B,0        ,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,1)
        IF (BC_2U .GT. 0) rhs(S12B:N12B,N2       ,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L .GT. 0) rhs(1        ,S23B:N23B,S33B:N33B,3) = bc31(S23B:N23B,S33B:N33B,1)
        IF (BC_1U .GT. 0) rhs(N1       ,S23B:N23B,S33B:N33B,3) = bc31(S23B:N23B,S33B:N33B,2)
        
        IF (BC_2L .GT. 0) rhs(S13B:N13B,1        ,S33B:N33B,3) = bc32(S13B:N13B,S33B:N33B,1)
        IF (BC_2U .GT. 0) rhs(S13B:N13B,N2       ,S33B:N33B,3) = bc32(S13B:N13B,S33B:N33B,2)
        
        IF (BC_3L .GT. 0) rhs(S13B:N13B,S23B:N23B,0        ,3) = bc33(S13B:N13B,S23B:N23B,1)
        IF (BC_3U .GT. 0) rhs(S13B:N13B,S23B:N23B,N3       ,3) = bc33(S13B:N13B,S23B:N23B,2)
        !-----------------------------------------------------------------------------------------------------
        
     END IF
     
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE rhs_vel
  
  
  
  
END MODULE mod_rhs
