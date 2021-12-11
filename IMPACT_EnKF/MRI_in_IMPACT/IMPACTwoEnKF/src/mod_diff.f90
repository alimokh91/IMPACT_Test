!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!*************************************************************************************************************

!> module containing subroutines for constructing difference operators and for interpolating values. It uses
!! modules mod_dims, mod_vars and mod_exchange.
MODULE mod_diff
  
  
  USE mod_dims
  USE mod_vars
  USE mod_exchange
  !USE mpi
  
  PRIVATE
  
  
  PUBLIC divergence, divergence2, divergence_transp
  PUBLIC gradient, gradient_transp
  PUBLIC Helmholtz, Helmholtz_explicit
  PUBLIC nonlinear, nonlinear_conc
  PUBLIC interpolate_vel
  PUBLIC outflow_bc
  PUBLIC bc_extrapolation, bc_extrapolation_transp
  PUBLIC interpolate_pre_vel , interpolate_vel_pre
  PUBLIC interpolate2_pre_vel, interpolate2_vel_pre
  PUBLIC first_pre_vel, first_vel_pre
  PUBLIC first_adv_pre, first_adv_vel
  PUBLIC Helmholtz_pre_explicit ! TEST!!!
  
  
  INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  ! TEST!!! Generell bei Randbehandlung bei ii=0,1 beginnen, bzw. bis ii=N1 rechnen!  
 
  !> subroutine that computes the divergece operator (on the pressure grid) in one spacial direction.
  SUBROUTINE divergence(m,phi,div)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)    ::  m                                           !< spacial dimension
  
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< ?
  REAL   , INTENT(inout) ::  div(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< divergence operator
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  CALL exchange(m,m,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 div(i,j,k) = cDu1(d1L,i)*phi(i+d1L,j,k)
!pgi$ unroll = n:8
                 DO ii = d1L+1, d1U
                    div(i,j,k) = div(i,j,k) + cDu1(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
!pgi$ unroll = n:8
                 DO jj = d2L, d2U
                    div(i,j,k) = div(i,j,k) + cDv2(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     !--------------------------------------------------------------------------------------------------------
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
!pgi$ unroll = n:8
                 DO kk = d3L, d3U
                    div(i,j,k) = div(i,j,k) + cDw3(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE divergence
  
  
  
  
  
  
  
  
  
  
  !> subroutine that computes the divergence (on the pressure grid) in all spacial dimensions 
  SUBROUTINE divergence2(phi,div)
  
  IMPLICIT NONE
  
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  REAL   , INTENT(  out) ::  div(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  CALL exchange(1,1,phi(b1L,b2L,b3L,1))
  CALL exchange(2,2,phi(b1L,b2L,b3L,2))
  CALL exchange(3,3,phi(b1L,b2L,b3L,3))
  
  
  !===========================================================================================================
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 div(i,j,k) = cDu1(d1L,i)*phi(i+d1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = d1L+1, d1U
                    div(i,j,k) = div(i,j,k) + cDu1(ii,i)*phi(i+ii,j,k,1)
                 END DO
!pgi$ unroll = n:8
                 DO jj = d2L, d2U
                    div(i,j,k) = div(i,j,k) + cDv2(jj,j)*phi(i,j+jj,k,2)
                 END DO
!pgi$ unroll = n:8
                 DO kk = d3L, d3U
                    div(i,j,k) = div(i,j,k) + cDw3(kk,k)*phi(i,j,k+kk,3)
                 END DO
              END DO
           END DO
        END DO
        
     !--------------------------------------------------------------------------------------------------------
     ELSE
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 div(i,j,k) = cDu1(d1L,i)*phi(i+d1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = d1L+1, d1U
                    div(i,j,k) = div(i,j,k) + cDu1(ii,i)*phi(i+ii,j,k,1)
                 END DO
!pgi$ unroll = n:8
                 DO jj = d2L, d2U
                    div(i,j,k) = div(i,j,k) + cDv2(jj,j)*phi(i,j+jj,k,2)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  END SUBROUTINE divergence2
  
  
  
  
  
  
  
  
  
  
  !> compute the transpose of the divergence operator 
  SUBROUTINE divergence_transp(m,phi,div)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  m                                          !< spacial direction
  
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< ?
  REAL   , INTENT(inout) ::  div(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< divergence (transpose) operator
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
        CALL exchange(m,0,phi)
        
        DO k = S31B, N31B ! "B"-Grenzen etwas konsistenter, aber nicht wirklich notwendig, da ohnehin mit bc_extrapolation nachmultipliziert wird ...
           DO j = S21B, N21B
              DO i = S11B, N11B
                 div(i,j,k) = cDu1T(g1L,i)*phi(i+g1L,j,k)
!pgi$ unroll = n:8
                 DO ii = g1L+1, g1U
                    div(i,j,k) = div(i,j,k) + cDu1T(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
        CALL exchange(m,0,phi)
        
        DO k = S32B, N32B
           DO j = S22B, N22B
              DO i = S12B, N12B
                 div(i,j,k) = cDv2T(g2L,j)*phi(i,j+g2L,k)
!pgi$ unroll = n:8
                 DO jj = g2L+1, g2U
                    div(i,j,k) = div(i,j,k) + cDv2T(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     !--------------------------------------------------------------------------------------------------------
        CALL exchange(m,0,phi)
        
        DO k = S33B, N33B
           DO j = S23B, N23B
              DO i = S13B, N13B
                 div(i,j,k) = cDw3T(g3L,k)*phi(i,j,k+g3L)
!pgi$ unroll = n:8
                 DO kk = g3L+1, g3U
                    div(i,j,k) = div(i,j,k) + cDw3T(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE divergence_transp
  
  
  
  
  
  
  
  
  
  
  !> subroutine that computes the gradient in direction m
  SUBROUTINE gradient(m,phi,grad)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  m                                           !< spacial direction
  
  REAL   , INTENT(inout) ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< ? 
  REAL   , INTENT(  out) ::  grad(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< gradient operator
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Randbedingungen könnten nur zum Teil in die Stencils eingebaut werden, so dass sich das   !
  !                vermutlich nicht wirklich lohnt.                                                          !
  !----------------------------------------------------------------------------------------------------------!
  
  
  CALL exchange(m,0,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
        DO k = S31, N31
           DO j = S21, N21
              DO i = S11, N11
                 grad(i,j,k) = cGp1(g1L,i)*phi(i+g1L,j,k)
!pgi$ unroll = n:8
                 DO ii = g1L+1, g1U
                    grad(i,j,k) = grad(i,j,k) + cGp1(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) grad(0 ,S21B:N21B,S31B:N31B) = 0.
     IF (BC_1U .GT. 0) grad(N1,S21B:N21B,S31B:N31B) = 0.
     IF (BC_2L .GT. 0) grad(S11B:N11B,1 ,S31B:N31B) = 0.
     IF (BC_2U .GT. 0) grad(S11B:N11B,N2,S31B:N31B) = 0.
     IF (BC_3L .GT. 0) grad(S11B:N11B,S21B:N21B,1 ) = 0.
     IF (BC_3U .GT. 0) grad(S11B:N11B,S21B:N21B,N3) = 0.
     
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
        DO k = S32, N32
           DO j = S22, N22
              DO i = S12, N12
                 grad(i,j,k) = cGp2(g2L,j)*phi(i,j+g2L,k)
!pgi$ unroll = n:8
                 DO jj = g2L+1, g2U
                    grad(i,j,k) = grad(i,j,k) + cGp2(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) grad(1 ,S22B:N22B,S32B:N32B) = 0.
     IF (BC_1U .GT. 0) grad(N1,S22B:N22B,S32B:N32B) = 0.
     IF (BC_2L .GT. 0) grad(S12B:N12B,0 ,S32B:N32B) = 0.
     IF (BC_2U .GT. 0) grad(S12B:N12B,N2,S32B:N32B) = 0.
     IF (BC_3L .GT. 0) grad(S12B:N12B,S22B:N22B,1 ) = 0.
     IF (BC_3U .GT. 0) grad(S12B:N12B,S22B:N22B,N3) = 0.
     
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     !--------------------------------------------------------------------------------------------------------
        DO k = S33, N33
           DO j = S23, N23
              DO i = S13, N13
                 grad(i,j,k) = cGp3(g3L,k)*phi(i,j,k+g3L)
!pgi$ unroll = n:8
                 DO kk = g3L+1, g3U
                    grad(i,j,k) = grad(i,j,k) + cGp3(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) grad(1 ,S23B:N23B,S33B:N33B) = 0.
     IF (BC_1U .GT. 0) grad(N1,S23B:N23B,S33B:N33B) = 0.
     IF (BC_2L .GT. 0) grad(S13B:N13B,1 ,S33B:N33B) = 0.
     IF (BC_2U .GT. 0) grad(S13B:N13B,N2,S33B:N33B) = 0.
     IF (BC_3L .GT. 0) grad(S13B:N13B,S23B:N23B,0 ) = 0.
     IF (BC_3U .GT. 0) grad(S13B:N13B,S23B:N23B,N3) = 0.
     
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE gradient
  
  
  
  
  
  
  
  
  
  
  !> computes the transpose of the gradient operator 
  SUBROUTINE gradient_transp(m,phi,grad)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  m                                           !< spacial direction
  
  REAL   , INTENT(inout) ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< ?
  REAL   , INTENT(  out) ::  grad(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< gradient (transpose) operator
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Umgekehrte Reihenfolge im Vergleich zu Subroutine "gradient".                             !
  !              - Randbedingungen könnten nur zum Teil in die Stencils eingebaut werden, so dass sich das   !
  !                vermutlich nicht wirklich lohnt.                                                          !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) phi(0 ,S21B:N21B,S31B:N31B) = 0.
     IF (BC_1U .GT. 0) phi(N1,S21B:N21B,S31B:N31B) = 0.
     IF (BC_2L .GT. 0) phi(S11B:N11B,1 ,S31B:N31B) = 0.
     IF (BC_2U .GT. 0) phi(S11B:N11B,N2,S31B:N31B) = 0.
     IF (BC_3L .GT. 0) phi(S11B:N11B,S21B:N21B,1 ) = 0.
     IF (BC_3U .GT. 0) phi(S11B:N11B,S21B:N21B,N3) = 0.
     
     !--------------------------------------------------------------------------------------------------------
        CALL exchange(m,m,phi) ! Muss nach Randbedingungen kommen!
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 grad(i,j,k) = cGp1T(d1L,i)*phi(i+d1L,j,k)
!pgi$ unroll = n:8
                 DO ii = d1L+1, d1U
                    grad(i,j,k) = grad(i,j,k) + cGp1T(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) phi(1 ,S22B:N22B,S32B:N32B) = 0.
     IF (BC_1U .GT. 0) phi(N1,S22B:N22B,S32B:N32B) = 0.
     IF (BC_2L .GT. 0) phi(S12B:N12B,0 ,S32B:N32B) = 0.
     IF (BC_2U .GT. 0) phi(S12B:N12B,N2,S32B:N32B) = 0.
     IF (BC_3L .GT. 0) phi(S12B:N12B,S22B:N22B,1 ) = 0.
     IF (BC_3U .GT. 0) phi(S12B:N12B,S22B:N22B,N3) = 0.
     
     !--------------------------------------------------------------------------------------------------------
        CALL exchange(m,m,phi) ! Muss nach Randbedingungen kommen!
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
!pgi$ unroll = n:8
                 DO jj = d2L, d2U
                    grad(i,j,k) = grad(i,j,k) + cGp2T(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) phi(1 ,S23B:N23B,S33B:N33B) = 0.
     IF (BC_1U .GT. 0) phi(N1,S23B:N23B,S33B:N33B) = 0.
     IF (BC_2L .GT. 0) phi(S13B:N13B,1 ,S33B:N33B) = 0.
     IF (BC_2U .GT. 0) phi(S13B:N13B,N2,S33B:N33B) = 0.
     IF (BC_3L .GT. 0) phi(S13B:N13B,S23B:N23B,0 ) = 0.
     IF (BC_3U .GT. 0) phi(S13B:N13B,S23B:N23B,N3) = 0.
     
     !--------------------------------------------------------------------------------------------------------
        CALL exchange(m,m,phi) ! Muss nach Randbedingungen kommen!
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
!pgi$ unroll = n:8
                 DO kk = d3L, d3U
                    grad(i,j,k) = grad(i,j,k) + cGp3T(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE gradient_transp
  
  
  
  
  
  
  
  
  
  
  !> subroutine that computes the Helmholtz operator.
  SUBROUTINE Helmholtz(m,exch_yes,phi,Lap)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  m                                           !< direction
  LOGICAL, INTENT(in   ) ::  exch_yes                                    !< whether to exchange between ghost cells 
  
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< ?
  REAL   , INTENT(  out) ::  Lap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< Laplace operator
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  REAL                   ::  dd1
  
  
  IF (exch_yes) THEN
     CALL exchange(1,m,phi)
     CALL exchange(2,m,phi)
     CALL exchange(3,m,phi)
  END IF
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
        IF (dimens == 3) THEN
           DO k = S31, N31
              DO j = S21, N21
                 DO i = S11, N11
                    dd1 = cu11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cu11(ii,i)*phi(i+ii,j,k)
                    END DO
!pgi$ unroll = n:8
                    DO jj = b2L, b2U
                       dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                    END DO
!pgi$ unroll = n:8
                    DO kk = b3L, b3U
                       dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
                    END DO
                    Lap(i,j,k) = phi(i,j,k) - multL*dd1
                 END DO
              END DO
           END DO
        ELSE
           DO k = S31, N31
              DO j = S21, N21
                 DO i = S11, N11
                    dd1 = cu11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cu11(ii,i)*phi(i+ii,j,k)
                    END DO
!pgi$ unroll = n:8
                    DO jj = b2L, b2U
                       dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                    END DO
                    Lap(i,j,k) = phi(i,j,k) - multL*dd1
                 END DO
              END DO
           END DO
        END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
        IF (dimens == 3) THEN
           DO k = S32, N32
              DO j = S22, N22
                 DO i = S12, N12
                    dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                    END DO
!pgi$ unroll = n:8
                    DO jj = b2L, b2U
                       dd1 = dd1 + cv22(jj,j)*phi(i,j+jj,k)
                    END DO
!pgi$ unroll = n:8
                    DO kk = b3L, b3U
                       dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
                    END DO
                    Lap(i,j,k) = phi(i,j,k) - multL*dd1
                 END DO
              END DO
           END DO
        ELSE
           DO k = S32, N32
              DO j = S22, N22
                 DO i = S12, N12
                    dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                    END DO
!pgi$ unroll = n:8
                    DO jj = b2L, b2U
                       dd1 = dd1 + cv22(jj,j)*phi(i,j+jj,k)
                    END DO
                    Lap(i,j,k) = phi(i,j,k) - multL*dd1
                 END DO
              END DO
           END DO
        END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 3 .AND. dimens == 3) THEN
     !--------------------------------------------------------------------------------------------------------
        DO k = S33, N33
           DO j = S23, N23
              DO i = S13, N13
                 dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                 END DO
!pgi$ unroll = n:8
                 DO kk = b3L, b3U
                    dd1 = dd1 + cw33(kk,k)*phi(i,j,k+kk)
                 END DO
                 Lap(i,j,k) = phi(i,j,k) - multL*dd1
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE Helmholtz
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE Helmholtz_explicit(exch_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  REAL                   ::  dd1
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Für den Aufbau der RHS bei expliziter Zeitintegration.                                    !
  !              - Randbedingungen müssen daher nicht berücksichtigt werden.                                 !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (exch_yes) CALL exchange_all_all(.TRUE.,vel)
  
  
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
        DO k = S31, N31
           DO j = S21, N21
              DO i = S11, N11
                 dd1 = cu11(b1L,i)*vel(i+b1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cu11(ii,i)*vel(i+ii,j,k,1)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*vel(i,j+jj,k,1)
                 END DO
!pgi$ unroll = n:8
                 DO kk = b3L, b3U
                    dd1 = dd1 + cp33(kk,k)*vel(i,j,k+kk,1)
                 END DO
                 nl(i,j,k,1) = nl(i,j,k,1) - multL*dd1
              END DO
           END DO
        END DO
     ELSE
        DO k = S31, N31
           DO j = S21, N21
              DO i = S11, N11
                 dd1 = cu11(b1L,i)*vel(i+b1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cu11(ii,i)*vel(i+ii,j,k,1)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*vel(i,j+jj,k,1)
                 END DO
                 nl(i,j,k,1) = nl(i,j,k,1) - multL*dd1
              END DO
           END DO
        END DO
     END IF
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
        DO k = S32, N32
           DO j = S22, N22
              DO i = S12, N12
                 dd1 = cp11(b1L,i)*vel(i+b1L,j,k,2)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*vel(i+ii,j,k,2)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cv22(jj,j)*vel(i,j+jj,k,2)
                 END DO
!pgi$ unroll = n:8
                 DO kk = b3L, b3U
                    dd1 = dd1 + cp33(kk,k)*vel(i,j,k+kk,2)
                 END DO
                 nl(i,j,k,2) = nl(i,j,k,2) - multL*dd1
              END DO
           END DO
        END DO
     ELSE
        DO k = S32, N32
           DO j = S22, N22
              DO i = S12, N12
                 dd1 = cp11(b1L,i)*vel(i+b1L,j,k,2)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*vel(i+ii,j,k,2)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cv22(jj,j)*vel(i,j+jj,k,2)
                 END DO
                 nl(i,j,k,2) = nl(i,j,k,2) - multL*dd1
              END DO
           END DO
        END DO
     END IF
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
        DO k = S33, N33
           DO j = S23, N23
              DO i = S13, N13
                 dd1 = cp11(b1L,i)*vel(i+b1L,j,k,3)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*vel(i+ii,j,k,3)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*vel(i,j+jj,k,3)
                 END DO
!pgi$ unroll = n:8
                 DO kk = b3L, b3U
                    dd1 = dd1 + cw33(kk,k)*vel(i,j,k+kk,3)
                 END DO
                 nl(i,j,k,3) = nl(i,j,k,3) - multL*dd1
              END DO
           END DO
        END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  END SUBROUTINE Helmholtz_explicit
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! relativ ungetestet ...
  SUBROUTINE Helmholtz_pre_explicit(exch_yes,phi,Lap)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(inout) ::  Lap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  REAL                   ::  dd1
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Für den Aufbau der RHS bei expliziter Zeitintegration.                                    !
  !              - Randbedingungen müssen daher nicht berücksichtigt werden.                                 !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (exch_yes) THEN
     CALL exchange(1,0,phi)
     CALL exchange(2,0,phi)
     CALL exchange(3,0,phi)
  END IF
  
  
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                 END DO
!pgi$ unroll = n:8
                 DO kk = b3L, b3U
                    dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
                 END DO
                 Lap(i,j,k) = Lap(i,j,k) - multL*dd1
              END DO
           END DO
        END DO
     ELSE
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                 END DO
                 Lap(i,j,k) = Lap(i,j,k) - multL*dd1
              END DO
           END DO
        END DO
     END IF
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  END SUBROUTINE Helmholtz_pre_explicit
  
  
  
  
  
  
  
  
  
  ! TEST!!! Teile davon (nur zentrale Operationen!) koennten jeweils durch interpolate2_pre_vel/interpolate2_vel_pre, first_adv_pre/first_adv_vel
  !         ersetzt werden (beachte aber Addition von nl!)
  ! TEST!!! umbenennen in advect... (?)
  !> computes advective (nonlinear) terms.
  SUBROUTINE nonlinear(exch_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  
  REAL                   ::  dd1
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - [Sij,Nij] ist immer eine Untermenge von [Sip,Nip]                                         !
  !              - Feld "res" wird mehrfach ausgetauscht, um mit einem statt drei skalaren Feldern arbeiten  !
  !                zu können. Im Prinzip könnte auch rhs(:,:,:,1:3) für die Zwischenspeicherung verwendet    !
  !                werden, jedoch sind dazu einige Umbaumassnahmen notwendig bei geringem Effizienzgewinn.   !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! worki muss bereits ausgetauscht sein!
  IF (exch_yes) CALL exchange_all_all(.TRUE.,vel)
  
  
  !===========================================================================================================
  !=== u*du/dx ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,1,work1,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S31, N31
        DO j = S21, N21
           DO i = S11, N11
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNu1U(n1L,i)*vel(i+n1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNu1U(ii,i)*vel(i+ii,j,k,1)
                 END DO
              ELSE
                 dd1 = cNu1D(n1L,i)*vel(i+n1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNu1D(ii,i)*vel(i+ii,j,k,1)
                 END DO
              END IF
              
              nl(i,j,k,1) = dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
        DO k = S31, N31
           DO j = S21, N21
              DO i = S11, N11
                 dd1 = cu1(b1L,i)*vel(i+b1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cu1(ii,i)*vel(i+ii,j,k,1)
                 END DO
                 
                 nl(i,j,k,1) = dd1*pp(i,j,k)
              END DO
           END DO
        END DO
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== v*du/dy ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,1,work2,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S31, N31
        DO j = S21, N21
           DO i = S11, N11
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNp2U(n2L,j)*vel(i,j+n2L,k,1)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNp2U(jj,j)*vel(i,j+jj,k,1)
                 END DO
              ELSE
                 dd1 = cNp2D(n2L,j)*vel(i,j+n2L,k,1)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNp2D(jj,j)*vel(i,j+jj,k,1)
                 END DO
              END IF
              
              nl(i,j,k,1) = nl(i,j,k,1) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
        DO k = S31, N31
           DO j = S21, N21
              DO i = S11, N11
                 dd1 = cp2(b2L,j)*vel(i,j+b2L,k,1)
!pgi$ unroll = n:8
                 DO jj = b2L+1, b2U
                    dd1 = dd1 + cp2(jj,j)*vel(i,j+jj,k,1)
                 END DO
                 
                 nl(i,j,k,1) = nl(i,j,k,1) + dd1*pp(i,j,k)
              END DO
           END DO
        END DO
  END IF
  !===========================================================================================================
  
  
  IF (dimens == 3) THEN
  !===========================================================================================================
  !=== w*du/dz ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,1,work3,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S31, N31
        DO j = S21, N21
           DO i = S11, N11
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNp3U(n3L,k)*vel(i,j,k+n3L,1)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNp3U(kk,k)*vel(i,j,k+kk,1)
                 END DO
              ELSE
                 dd1 = cNp3D(n3L,k)*vel(i,j,k+n3L,1)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNp3D(kk,k)*vel(i,j,k+kk,1)
                 END DO
              END IF
              
              nl(i,j,k,1) = nl(i,j,k,1) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
        DO k = S31, N31
           DO j = S21, N21
              DO i = S11, N11
                 dd1 = cp3(b3L,k)*vel(i,j,k+b3L,1)
!pgi$ unroll = n:8
                 DO kk = b3L+1, b3U
                    dd1 = dd1 + cp3(kk,k)*vel(i,j,k+kk,1)
                 END DO
                 
                 nl(i,j,k,1) = nl(i,j,k,1) + dd1*pp(i,j,k)
              END DO
           END DO
        END DO
  END IF
  !===========================================================================================================
  END IF
  
  
  
  
  
  
  !===========================================================================================================
  !=== u*dv/dx ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,2,work1,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S32, N32
        DO j = S22, N22
           DO i = S12, N12
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNp1U(n1L,i)*vel(i+n1L,j,k,2)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNp1U(ii,i)*vel(i+ii,j,k,2)
                 END DO
              ELSE
                 dd1 = cNp1D(n1L,i)*vel(i+n1L,j,k,2)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNp1D(ii,i)*vel(i+ii,j,k,2)
                 END DO
              END IF
              
              nl(i,j,k,2) = dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
        DO k = S32, N32
           DO j = S22, N22
              DO i = S12, N12
                 dd1 = cp1(b1L,i)*vel(i+b1L,j,k,2)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp1(ii,i)*vel(i+ii,j,k,2)
                 END DO
                 
                 nl(i,j,k,2) = dd1*pp(i,j,k)
              END DO
           END DO
        END DO
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== v*dv/dy ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,2,work2,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S32, N32
        DO j = S22, N22
           DO i = S12, N12
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNv2U(n2L,j)*vel(i,j+n2L,k,2)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNv2U(jj,j)*vel(i,j+jj,k,2)
                 END DO
              ELSE
                 dd1 = cNv2D(n2L,j)*vel(i,j+n2L,k,2)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNv2D(jj,j)*vel(i,j+jj,k,2)
                 END DO
              END IF
              
              nl(i,j,k,2) = nl(i,j,k,2) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
        DO k = S32, N32
           DO j = S22, N22
              DO i = S12, N12
                 dd1 = cv2(b2L,j)*vel(i,j+b2L,k,2)
!pgi$ unroll = n:8
                 DO jj = b2L+1, b2U
                    dd1 = dd1 + cv2(jj,j)*vel(i,j+jj,k,2)
                 END DO
                 
                 nl(i,j,k,2) = nl(i,j,k,2) + dd1*pp(i,j,k)
              END DO
           END DO
        END DO
  END IF
  !===========================================================================================================
  
  
  IF (dimens == 3) THEN
  !===========================================================================================================
  !=== w*dv/dz ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,2,work3,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S32, N32
        DO j = S22, N22
           DO i = S12, N12
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNp3U(n3L,k)*vel(i,j,k+n3L,2)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNp3U(kk,k)*vel(i,j,k+kk,2)
                 END DO
              ELSE
                 dd1 = cNp3D(n3L,k)*vel(i,j,k+n3L,2)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNp3D(kk,k)*vel(i,j,k+kk,2)
                 END DO
              END IF
              
              nl(i,j,k,2) = nl(i,j,k,2) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
        DO k = S32, N32
           DO j = S22, N22
              DO i = S12, N12
                 dd1 = cp3(b3L,k)*vel(i,j,k+b3L,2)
!pgi$ unroll = n:8
                 DO kk = b3L+1, b3U
                    dd1 = dd1 + cp3(kk,k)*vel(i,j,k+kk,2)
                 END DO
                 
                 nl(i,j,k,2) = nl(i,j,k,2) + dd1*pp(i,j,k)
              END DO
           END DO
        END DO
  END IF
  !===========================================================================================================
  END IF
  
  
  
  
  
  IF (dimens == 3) THEN
  !===========================================================================================================
  !=== u*dw/dx ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,3,work1,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S33, N33
        DO j = S23, N23
           DO i = S13, N13
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNp1U(n1L,i)*vel(i+n1L,j,k,3)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNp1U(ii,i)*vel(i+ii,j,k,3)
                 END DO
              ELSE
                 dd1 = cNp1D(n1L,i)*vel(i+n1L,j,k,3)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNp1D(ii,i)*vel(i+ii,j,k,3)
                 END DO
              END IF
              
              nl(i,j,k,3) = dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
        DO k = S33, N33
           DO j = S23, N23
              DO i = S13, N13
                 dd1 = cp1(b1L,i)*vel(i+b1L,j,k,3)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp1(ii,i)*vel(i+ii,j,k,3)
                 END DO
                 
                 nl(i,j,k,3) = dd1*pp(i,j,k)
              END DO
           END DO
        END DO
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== v*dw/dy ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,3,work2,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S33, N33
        DO j = S23, N23
           DO i = S13, N13
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNp2U(n2L,j)*vel(i,j+n2L,k,3)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNp2U(jj,j)*vel(i,j+jj,k,3)
                 END DO
              ELSE
                 dd1 = cNp2D(n2L,j)*vel(i,j+n2L,k,3)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNp2D(jj,j)*vel(i,j+jj,k,3)
                 END DO
              END IF
              
              nl(i,j,k,3) = nl(i,j,k,3) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
        DO k = S33, N33
           DO j = S23, N23
              DO i = S13, N13
                 dd1 = cp2(b2L,j)*vel(i,j+b2L,k,3)
!pgi$ unroll = n:8
                 DO jj = b2L+1, b2U
                    dd1 = dd1 + cp2(jj,j)*vel(i,j+jj,k,3)
                 END DO
                 
                 nl(i,j,k,3) = nl(i,j,k,3) + dd1*pp(i,j,k)
              END DO
           END DO
        END DO
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== w*dw/dz ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,3,work3,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S33, N33
        DO j = S23, N23
           DO i = S13, N13
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNw3U(n3L,k)*vel(i,j,k+n3L,3)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNw3U(kk,k)*vel(i,j,k+kk,3)
                 END DO
              ELSE
                 dd1 = cNw3D(n3L,k)*vel(i,j,k+n3L,3)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNw3D(kk,k)*vel(i,j,k+kk,3)
                 END DO
              END IF
              
              nl(i,j,k,3) = nl(i,j,k,3) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
        DO k = S33, N33
           DO j = S23, N23
              DO i = S13, N13
                 dd1 = cw3(b3L,k)*vel(i,j,k+b3L,3)
!pgi$ unroll = n:8
                 DO kk = b3L+1, b3U
                    dd1 = dd1 + cw3(kk,k)*vel(i,j,k+kk,3)
                 END DO
                 
                 nl(i,j,k,3) = nl(i,j,k,3) + dd1*pp(i,j,k)
              END DO
           END DO
        END DO
  END IF
  !===========================================================================================================
  END IF
  
  
  END SUBROUTINE nonlinear
  
  
  
  
  
  
  
  
  
  ! TEST!!! anderes Modul?
  ! TEST!!! basiert neu auf interpolate2_vel_pre ... ok??
  SUBROUTINE interpolate_vel(exch_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: -                                                                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
                   CALL interpolate2_vel_pre(exch_yes,1,vel(b1L,b2L,b3L,1),work1)
                   CALL interpolate2_vel_pre(exch_yes,2,vel(b1L,b2L,b3L,2),work2)
  IF (dimens == 3) CALL interpolate2_vel_pre(exch_yes,3,vel(b1L,b2L,b3L,3),work3)
  
  
  CALL exchange(1,0,work1)
  CALL exchange(2,0,work1)
  CALL exchange(3,0,work1)
  
  CALL exchange(1,0,work2)
  CALL exchange(2,0,work2)
  CALL exchange(3,0,work2)
  
  IF (dimens == 3) THEN
     CALL exchange(1,0,work3)
     CALL exchange(2,0,work3)
     CALL exchange(3,0,work3)
  END IF
  
  
  END SUBROUTINE interpolate_vel
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! noch relativ ungetestet!
  SUBROUTINE interpolate_pre_vel(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,inter)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  INTEGER, INTENT(in   ) ::  m
  INTEGER, INTENT(in   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  
  REAL   , INTENT(inout) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(inout) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !              - die Punkte auf (bzw. hinter) dem Rand der Wand-normalen Komponente werden auch bei        !
  !                kompakter Differenzierung im Feld immer explizit gerechnet, um nur eine Variante          !
  !                kompakter Differenzen abspeichern zu muessen (Extrapolation!).                            !
  !              - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (exch_yes) CALL exchange2(m,0,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
        DO k = SS3, NN3
           DO j = SS2, NN2
              DO i = S11B, N11B ! TEST!!! hier koennte man auch SS1, NN1 nehmen! gilt auch fuer andere Routinen!!
                 inter(i,j,k) = cIpu(g1L,i)*phi(i+g1L,j,k)
                 DO ii = g1L+1, g1U
                    inter(i,j,k) = inter(i,j,k) + cIpu(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN
        DO k = SS3, NN3
           DO j = S22B, N22B ! TEST!!! hier koennte man auch SS2, NN2 nehmen!
              DO i = SS1, NN1
                 inter(i,j,k) = cIpv(g2L,j)*phi(i,j+g2L,k)
                 DO jj = g2L+1, g2U
                    inter(i,j,k) = inter(i,j,k) + cIpv(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN
        DO k = S33B, N33B ! TEST!!! hier koennte man auch SS3, NN3 nehmen!
           DO j = SS2, NN2
              DO i = SS1, NN1
                 inter(i,j,k) = cIpw(g3L,k)*phi(i,j,k+g3L)
                 DO kk = g3L+1, g3U
                    inter(i,j,k) = inter(i,j,k) + cIpw(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE interpolate_pre_vel
  
  
  
  
  
  
  
  
  
  
  
  ! Wie interpolate_pre_vel, allerdings mit fixen Index-Limiten (ohne Rand)
  SUBROUTINE interpolate2_pre_vel(exch_yes,m,phi,inter)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  INTEGER, INTENT(in   ) ::  m
  
  REAL   , INTENT(inout) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(inout) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  IF (exch_yes) CALL exchange(m,0,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
        DO k = S31, N31
           DO j = S21, N21
              DO i = S11, N11
                 inter(i,j,k) = cIpu(g1L,i)*phi(i+g1L,j,k)
                 DO ii = g1L+1, g1U
                    inter(i,j,k) = inter(i,j,k) + cIpu(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN
        DO k = S32, N32
           DO j = S22, N22
              DO i = S12, N12
                 inter(i,j,k) = cIpv(g2L,j)*phi(i,j+g2L,k)
                 DO jj = g2L+1, g2U
                    inter(i,j,k) = inter(i,j,k) + cIpv(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN
        DO k = S33, N33
           DO j = S23, N23
              DO i = S13, N13
                 inter(i,j,k) = cIpw(g3L,k)*phi(i,j,k+g3L)
                 DO kk = g3L+1, g3U
                    inter(i,j,k) = inter(i,j,k) + cIpw(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE interpolate2_pre_vel
  
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! noch relativ ungetestet!
  SUBROUTINE interpolate_vel_pre(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,inter)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  INTEGER, INTENT(in)    ::  m
  INTEGER, INTENT(in   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  
  REAL   , INTENT(inout) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(inout) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (exch_yes) CALL exchange2(m,m,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
        DO k = SS3, NN3
           DO j = SS2, NN2
              DO i = S1p, N1p
                 inter(i,j,k) = cIup(d1L,i)*phi(i+d1L,j,k)
                 DO ii = d1L+1, d1U
                    inter(i,j,k) = inter(i,j,k) + cIup(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN
        DO k = SS3, NN3
           DO j = S2p, N2p
              DO i = SS1, NN1
                 inter(i,j,k) = cIvp(d2L,j)*phi(i,j+d2L,k)
                 DO jj = d2L+1, d2U
                    inter(i,j,k) = inter(i,j,k) + cIvp(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN
        DO k = S3p, N3p
           DO j = SS2, NN2
              DO i = SS1, NN1
                 inter(i,j,k) = cIwp(d3L,k)*phi(i,j,k+d3L)
                 DO kk = d3L+1, d3U
                    inter(i,j,k) = inter(i,j,k) + cIwp(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE interpolate_vel_pre
  
  
  
  
  
  
  
  
  
  
  
  ! Wie interpolate_vel_pre, allerdings mit fixen Index-Limiten (ohne Rand)
  SUBROUTINE interpolate2_vel_pre(exch_yes,m,phi,inter)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  INTEGER, INTENT(in)    ::  m
  
  REAL   , INTENT(inout) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(inout) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  IF (exch_yes) CALL exchange(m,m,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 inter(i,j,k) = cIup(d1L,i)*phi(i+d1L,j,k)
                 DO ii = d1L+1, d1U
                    inter(i,j,k) = inter(i,j,k) + cIup(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 inter(i,j,k) = cIvp(d2L,j)*phi(i,j+d2L,k)
                 DO jj = d2L+1, d2U
                    inter(i,j,k) = inter(i,j,k) + cIvp(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 inter(i,j,k) = cIwp(d3L,k)*phi(i,j,k+d3L)
                 DO kk = d3L+1, d3U
                    inter(i,j,k) = inter(i,j,k) + cIwp(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE interpolate2_vel_pre
  
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! noch relativ ungetestet!
  SUBROUTINE first_adv_pre(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,der,adv,upwind_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  INTEGER, INTENT(in   ) ::  m
  INTEGER, INTENT(in   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  LOGICAL, INTENT(in   ) ::  upwind_yes
  
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(inout) ::  der(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(in   ) ::  adv(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  REAL                   ::  dd1
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - die Punkte auf dem Rand der Wand-normalen Komponente werden auch bei kompakter            !
  !                Differenzierung im Feld immer explizit gerechnet, um nur eine Variante kompakter          !
  !                Differenzen abspeichern zu muessen.                                                       !
  !              - fuer die Punkte der Wand-normalen Komponente sind S12, N12 nur Platzhalter fuer S1p+1,    !
  !                N1p-1 usw.                                                                                !
  !              - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (exch_yes) CALL exchange2(m,0,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (upwind_yes) THEN
        DO k = SS3, NN3
           DO j = SS2, NN2
              DO i = S1p, N1p
                 IF (adv(i,j,k) .GE. 0.) THEN
                    dd1 = cNp1U(n1L,i)*phi(i+n1L,j,k)
                    DO ii = n1L+1, n1U
                       dd1 = dd1 + cNp1U(ii,i)*phi(i+ii,j,k)
                    END DO
                 ELSE
                    dd1 = cNp1D(n1L,i)*phi(i+n1L,j,k)
                    DO ii = n1L+1, n1U
                       dd1 = dd1 + cNp1D(ii,i)*phi(i+ii,j,k)
                    END DO
                 END IF
                 der(i,j,k) = dd1*adv(i,j,k)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
           DO k = SS3, NN3
              DO j = SS2, NN2
                 DO i = S1p, N1p
                    dd1 = cp1(b1L,i)*phi(i+b1L,j,k)
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cp1(ii,i)*phi(i+ii,j,k)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (upwind_yes) THEN
        DO k = SS3, NN3
           DO j = S2p, N2p
              DO i = SS1, NN1
                 IF (adv(i,j,k) .GE. 0.) THEN
                    dd1 = cNp2U(n2L,j)*phi(i,j+n2L,k)
                    DO jj = n2L+1, n2U
                       dd1 = dd1 + cNp2U(jj,j)*phi(i,j+jj,k)
                    END DO
                 ELSE
                    dd1 = cNp2D(n2L,j)*phi(i,j+n2L,k)
                    DO jj = n2L+1, n2U
                       dd1 = dd1 + cNp2D(jj,j)*phi(i,j+jj,k)
                    END DO
                 END IF
                 der(i,j,k) = dd1*adv(i,j,k)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
           DO k = SS3, NN3
              DO j = S2p, N2p
                 DO i = SS1, NN1
                    dd1 = cp2(b2L,j)*phi(i,j+b2L,k)
                    DO jj = b2L+1, b2U
                       dd1 = dd1 + cp2(jj,j)*phi(i,j+jj,k)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (upwind_yes) THEN
        DO k = S3p, N3p
           DO j = SS2, NN2
              DO i = SS1, NN1
                 IF (adv(i,j,k) .GE. 0.) THEN
                    dd1 = cNp3U(n3L,k)*phi(i,j,k+n3L)
                    DO kk = n3L+1, n3U
                       dd1 = dd1 + cNp3U(kk,k)*phi(i,j,k+kk)
                    END DO
                 ELSE
                    dd1 = cNp3D(n3L,k)*phi(i,j,k+n3L)
                    DO kk = n3L+1, n3U
                       dd1 = dd1 + cNp3D(kk,k)*phi(i,j,k+kk)
                    END DO
                 END IF
                 der(i,j,k) = dd1*adv(i,j,k)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
           DO k = S3p, N3p
              DO j = SS2, NN2
                 DO i = SS1, NN1
                    dd1 = cp3(b3L,k)*phi(i,j,k+b3L)
                    DO kk = b3L+1, b3U
                       dd1 = dd1 + cp3(kk,k)*phi(i,j,k+kk)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE first_adv_pre
  
  
  
  
  
  
  
  
  
  ! TEST!!! Routine kann Punkte der Rand-normalen Komponente auf dem Rand NICHT behandeln!! (will sie momentan aber auch gar nicht koennen)
  ! TEST!!! noch relativ ungetestet!
  SUBROUTINE first_adv_vel(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,der,adv,upwind_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  INTEGER, INTENT(in   ) ::  m
  INTEGER, INTENT(in   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  LOGICAL, INTENT(in   ) ::  upwind_yes
  
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(inout) ::  der(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(in   ) ::  adv(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  REAL                   ::  dd1
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (exch_yes) CALL exchange2(m,m,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (upwind_yes) THEN
        DO k = SS3, NN3
           DO j = SS2, NN2
              DO i = S11, N11
                 IF (adv(i,j,k) .GE. 0.) THEN
                    dd1 = cNu1U(n1L,i)*phi(i+n1L,j,k)
                    DO ii = n1L+1, n1U
                       dd1 = dd1 + cNu1U(ii,i)*phi(i+ii,j,k)
                    END DO
                 ELSE
                    dd1 = cNu1D(n1L,i)*phi(i+n1L,j,k)
                    DO ii = n1L+1, n1U
                       dd1 = dd1 + cNu1D(ii,i)*phi(i+ii,j,k)
                    END DO
                 END IF
                 der(i,j,k) = dd1*adv(i,j,k)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
           DO k = SS3, NN3
              DO j = SS2, NN2
                 DO i = S11, N11
                    dd1 = cu1(b1L,i)*phi(i+b1L,j,k)
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cu1(ii,i)*phi(i+ii,j,k)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (upwind_yes) THEN
        DO k = SS3, NN3
           DO j = S22, N22
              DO i = SS1, NN1
                 IF (adv(i,j,k) .GE. 0.) THEN
                    dd1 = cNv2U(n2L,j)*phi(i,j+n2L,k)
                    DO jj = n2L+1, n2U
                       dd1 = dd1 + cNv2U(jj,j)*phi(i,j+jj,k)
                    END DO
                 ELSE
                    dd1 = cNv2D(n2L,j)*phi(i,j+n2L,k)
                    DO jj = n2L+1, n2U
                       dd1 = dd1 + cNv2D(jj,j)*phi(i,j+jj,k)
                    END DO
                 END IF
                 der(i,j,k) = dd1*adv(i,j,k)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
           DO k = SS3, NN3
              DO j = S22, N22
                 DO i = SS1, NN1
                    dd1 = cv2(b2L,j)*phi(i,j+b2L,k)
                    DO jj = b2L+1, b2U
                       dd1 = dd1 + cv2(jj,j)*phi(i,j+jj,k)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (upwind_yes) THEN
        DO k = S33, N33
           DO j = SS2, NN2
              DO i = SS1, NN1
                 IF (adv(i,j,k) .GE. 0.) THEN
                    dd1 = cNw3U(n3L,k)*phi(i,j,k+n3L)
                    DO kk = n3L+1, n3U
                       dd1 = dd1 + cNw3U(kk,k)*phi(i,j,k+kk)
                    END DO
                 ELSE
                    dd1 = cNw3D(n3L,k)*phi(i,j,k+n3L)
                    DO kk = n3L+1, n3U
                       dd1 = dd1 + cNw3D(kk,k)*phi(i,j,k+kk)
                    END DO
                 END IF
                 der(i,j,k) = dd1*adv(i,j,k)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
           DO k = S33, N33
              DO j = SS2, NN2
                 DO i = SS1, NN1
                    dd1 = cw3(b3L,k)*phi(i,j,k+b3L)
                    DO kk = b3L+1, b3U
                       dd1 = dd1 + cw3(kk,k)*phi(i,j,k+kk)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE first_adv_vel
  
  
  
  
  
  
  
  
  
  
  !> subroutine that sets the boundary condition for advective boundaries 
  SUBROUTINE outflow_bc
  
  ! (revised on 06.08.2009)
  
  IMPLICIT NONE
  
  INTEGER                ::  m
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  REAL                   ::  dd1
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Es wird vorausgesetzt, dass vel(:,:,:,:) zuvor schon ausgetauscht, bzw. an den Raendern   !
  !                zu Null gesetzt wurde (was beides streng genommen aber nichtmal notwentwendig ist).       !
  !              - Kompakte Differenzen nicht notwendig, ist ohnehin nur ein Modell.                         !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== vel(:,:,:,1) ==========================================================================================
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  !--- Tangentialkomponenten ---------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(2,1,1) .AND. (BC_2L == 1 .OR. BC_2L == 3)) THEN ! TEST!!! Warum eigentlich BC_2L == 3? Wuerde generell nicht outlet(2,1,1) ausreichen?
     j = 1
     DO k = S31B, N31B
        DO i = S11B, N11B
           nlbc12(i,k,1) = cp2(b2L,j)*vel(i,j+b2L,k,1)
!pgi$ unroll = n:8
           DO jj = b2L+1, b2U
              nlbc12(i,k,1) = nlbc12(i,k,1) + cp2(jj,j)*vel(i,j+jj,k,1)
           END DO
           
           dd1 = cIpu(g1L,i)*drift2(i+g1L,k,1)
!pgi$ unroll = n:8
           DO ii = g1L+1, g1U
              dd1 = dd1 + cIpu(ii,i)*drift2(i+ii,k,1)
           END DO
           nlbc12(i,k,1) = dd1*nlbc12(i,k,1)
        END DO
     END DO
  END IF
  IF (outlet(2,2,1) .AND. (BC_2U == 1 .OR. BC_2U == 3)) THEN
     j = N2
     DO k = S31B, N31B
        DO i = S11B, N11B
           nlbc12(i,k,2) = cp2(b2L,j)*vel(i,j+b2L,k,1)
!pgi$ unroll = n:8
           DO jj = b2L+1, b2U
              nlbc12(i,k,2) = nlbc12(i,k,2) + cp2(jj,j)*vel(i,j+jj,k,1)
           END DO
           
           dd1 = cIpu(g1L,i)*drift2(i+g1L,k,2)
!pgi$ unroll = n:8
           DO ii = g1L+1, g1U
              dd1 = dd1 + cIpu(ii,i)*drift2(i+ii,k,2)
           END DO
           nlbc12(i,k,2) = dd1*nlbc12(i,k,2)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(3,1,1) .AND. (BC_3L == 1 .OR. BC_3L == 3)) THEN
     k = 1
     DO j = S21B, N21B
        DO i = S11B, N11B
           nlbc13(i,j,1) = cp3(b3L,k)*vel(i,j,k+b3L,1)
!pgi$ unroll = n:8
           DO kk = b3L+1, b3U
              nlbc13(i,j,1) = nlbc13(i,j,1) + cp3(kk,k)*vel(i,j,k+kk,1)
           END DO
           
           dd1 = cIpu(g1L,i)*drift3(i+g1L,j,1)
!pgi$ unroll = n:8
           DO ii = g1L+1, g1U
              dd1 = dd1 + cIpu(ii,i)*drift3(i+ii,j,1)
           END DO
           nlbc13(i,j,1) = dd1*nlbc13(i,j,1)
        END DO
     END DO
  END IF
  IF (outlet(3,2,1) .AND. (BC_3U == 1 .OR. BC_3U == 3)) THEN
     k = N3
     DO j = S21B, N21B
        DO i = S11B, N11B
           nlbc13(i,j,2) = cp3(b3L,k)*vel(i,j,k+b3L,1)
!pgi$ unroll = n:8
           DO kk = b3L+1, b3U
              nlbc13(i,j,2) = nlbc13(i,j,2) + cp3(kk,k)*vel(i,j,k+kk,1)
           END DO
           
           dd1 = cIpu(g1L,i)*drift3(i+g1L,j,2)
!pgi$ unroll = n:8
           DO ii = g1L+1, g1U
              dd1 = dd1 + cIpu(ii,i)*drift3(i+ii,j,2)
           END DO
           nlbc13(i,j,2) = dd1*nlbc13(i,j,2)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !--- Normalkomponente --------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(1,1,1) .AND. (BC_1L == 1 .OR. BC_1L == 2)) THEN
     i = 0
     DO k = S31B, N31B
        DO j = S21B, N21B
           nlbc11(j,k,1) = cDu1(d1L,1)*vel(1+d1L,j,k,1)
!pgi$ unroll = n:8
           DO ii = d1L+1, d1U
              nlbc11(j,k,1) = nlbc11(j,k,1) + cDu1(ii,1)*vel(1+ii,j,k,1)
           END DO
           nlbc11(j,k,1) = drift1(j,k,1)*nlbc11(j,k,1)
        END DO
     END DO
  END IF
  IF (outlet(1,2,1) .AND. (BC_1U == 1 .OR. BC_1U == 2)) THEN
     i = N1
     DO k = S31B, N31B
        DO j = S21B, N21B
           nlbc11(j,k,2) = cDu1(d1L,i)*vel(i+d1L,j,k,1)
!pgi$ unroll = n:8
           DO ii = d1L+1, d1U
              nlbc11(j,k,2) = nlbc11(j,k,2) + cDu1(ii,i)*vel(i+ii,j,k,1)
           END DO
           nlbc11(j,k,2) = drift1(j,k,2)*nlbc11(j,k,2)
        END DO
     END DO
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== vel(:,:,:,2) ==========================================================================================
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  !--- Tangentialkomponenten ---------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(1,1,2) .AND. (BC_1L == 1 .OR. BC_1L == 3)) THEN
     i = 1
     DO k = S32B, N32B
        DO j = S22B, N22B
           nlbc21(j,k,1) = cp1(b1L,i)*vel(i+b1L,j,k,2)
!pgi$ unroll = n:8
           DO ii = b1L+1, b1U
              nlbc21(j,k,1) = nlbc21(j,k,1) + cp1(ii,i)*vel(i+ii,j,k,2)
           END DO
           
           dd1 = cIpv(g2L,j)*drift1(j+g2L,k,1)
!pgi$ unroll = n:8
           DO jj = g2L+1, g2U
              dd1 = dd1 + cIpv(jj,j)*drift1(j+jj,k,1)
           END DO
           nlbc21(j,k,1) = dd1*nlbc21(j,k,1)
        END DO
     END DO
  END IF
  IF (outlet(1,2,2) .AND. (BC_1U == 1 .OR. BC_1U == 3)) THEN
     i = N1
     DO k = S32B, N32B
        DO j = S22B, N22B
           nlbc21(j,k,2) = cp1(b1L,i)*vel(i+b1L,j,k,2)
!pgi$ unroll = n:8
           DO ii = b1L+1, b1U
              nlbc21(j,k,2) = nlbc21(j,k,2) + cp1(ii,i)*vel(i+ii,j,k,2)
           END DO
           
           dd1 = cIpv(g2L,j)*drift1(j+g2L,k,2)
!pgi$ unroll = n:8
           DO jj = g2L+1, g2U
              dd1 = dd1 + cIpv(jj,j)*drift1(j+jj,k,2)
           END DO
           nlbc21(j,k,2) = dd1*nlbc21(j,k,2)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(3,1,2) .AND. (BC_3L == 1 .OR. BC_3L == 3)) THEN
     k = 1
     DO j = S22B, N22B
        DO i = S12B, N12B
           nlbc23(i,j,1) = cp3(b3L,k)*vel(i,j,k+b3L,2)
!pgi$ unroll = n:8
           DO kk = b3L+1, b3U
              nlbc23(i,j,1) = nlbc23(i,j,1) + cp3(kk,k)*vel(i,j,k+kk,2)
           END DO
           
           dd1 = cIpv(g2L,j)*drift3(i,j+g2L,1)
!pgi$ unroll = n:8
           DO jj = g2L+1, g2U
              dd1 = dd1 + cIpv(jj,j)*drift3(i,j+jj,1)
           END DO
           nlbc23(i,j,1) = dd1*nlbc23(i,j,1)
        END DO
     END DO
  END IF
  IF (outlet(3,2,2) .AND. (BC_3U == 1 .OR. BC_3U == 3)) THEN
     k = N3
     DO j = S22B, N22B
        DO i = S12B, N12B
           nlbc23(i,j,2) = cp3(b3L,k)*vel(i,j,k+b3L,2)
!pgi$ unroll = n:8
           DO kk = b3L+1, b3U
              nlbc23(i,j,2) = nlbc23(i,j,2) + cp3(kk,k)*vel(i,j,k+kk,2)
           END DO
           
           dd1 = cIpv(g2L,j)*drift3(i,j+g2L,2)
!pgi$ unroll = n:8
           DO jj = g2L+1, g2U
              dd1 = dd1 + cIpv(jj,j)*drift3(i,j+jj,2)
           END DO
           nlbc23(i,j,2) = dd1*nlbc23(i,j,2)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !--- Normalkomponente --------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(2,1,2) .AND. (BC_2L == 1 .OR. BC_2L == 2)) THEN
     j = 0
     DO k = S32B, N32B
        DO i = S12B, N12B
           nlbc22(i,k,1) = cDv2(d2L,1)*vel(i,1+d2L,k,2)
!pgi$ unroll = n:8
           DO jj = d2L+1, d2U
              nlbc22(i,k,1) = nlbc22(i,k,1) + cDv2(jj,1)*vel(i,1+jj,k,2)
           END DO
           nlbc22(i,k,1) = drift2(i,k,1)*nlbc22(i,k,1)
        END DO
     END DO
  END IF
  IF (outlet(2,2,2) .AND. (BC_2U == 1 .OR. BC_2U == 2)) THEN
     j = N2
     DO k = S32B, N32B
        DO i = S12B, N12B
           nlbc22(i,k,2) = cDv2(d2L,j)*vel(i,j+d2L,k,2)
!pgi$ unroll = n:8
           DO jj = d2L+1, d2U
              nlbc22(i,k,2) = nlbc22(i,k,2) + cDv2(jj,j)*vel(i,j+jj,k,2)
           END DO
           nlbc22(i,k,2) = drift2(i,k,2)*nlbc22(i,k,2)
        END DO
     END DO
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== vel(:,:,:,3) ==========================================================================================
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  !--- Tangentialkomponenten ---------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(1,1,3) .AND. (BC_1L == 1 .OR. BC_1L == 3)) THEN
     i = 1
     DO k = S33B, N33B
        DO j = S23B, N23B
           nlbc31(j,k,1) = cp1(b1L,i)*vel(i+b1L,j,k,3)
!pgi$ unroll = n:8
           DO ii = b1L+1, b1U
              nlbc31(j,k,1) = nlbc31(j,k,1) + cp1(ii,i)*vel(i+ii,j,k,3)
           END DO
           
           dd1 = cIpw(g3L,k)*drift1(j,k+g3L,1)
!pgi$ unroll = n:8
           DO kk = g3L+1, g3U
              dd1 = dd1 + cIpw(kk,k)*drift1(j,k+kk,1)
           END DO
           nlbc31(j,k,1) = dd1*nlbc31(j,k,1)
        END DO
     END DO
  END IF
  IF (outlet(1,2,3) .AND. (BC_1U == 1 .OR. BC_1U == 3)) THEN
     i = N1
     DO k = S33B, N33B
        DO j = S23B, N23B
           nlbc31(j,k,2) = cp1(b1L,i)*vel(i+b1L,j,k,3)
!pgi$ unroll = n:8
           DO ii = b1L+1, b1U
              nlbc31(j,k,2) = nlbc31(j,k,2) + cp1(ii,i)*vel(i+ii,j,k,3)
           END DO
           
           dd1 = cIpw(g3L,k)*drift1(j,k+g3L,2)
!pgi$ unroll = n:8
           DO kk = g3L+1, g3U
              dd1 = dd1 + cIpw(kk,k)*drift1(j,k+kk,2)
           END DO
           nlbc31(j,k,2) = dd1*nlbc31(j,k,2)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(2,1,3) .AND. (BC_2L == 1 .OR. BC_2L == 3)) THEN
     j = 1
     DO k = S33B, N33B
        DO i = S13B, N13B
           nlbc32(i,k,1) = cp2(b2L,j)*vel(i,j+b2L,k,3)
!pgi$ unroll = n:8
           DO jj = b2L+1, b2U
              nlbc32(i,k,1) = nlbc32(i,k,1) + cp2(jj,j)*vel(i,j+jj,k,3)
           END DO
           
           dd1 = cIpw(g3L,k)*drift2(i,k+g3L,1)
!pgi$ unroll = n:8
           DO kk = g3L+1, g3U
              dd1 = dd1 + cIpw(kk,k)*drift2(i,k+kk,1)
           END DO
           nlbc32(i,k,1) = dd1*nlbc32(i,k,1)
        END DO
     END DO
  END IF
  IF (outlet(2,2,3) .AND. (BC_2U == 1 .OR. BC_2U == 3)) THEN
     j = N2
     DO k = S33B, N33B
        DO i = S13B, N13B
           nlbc32(i,k,2) = cp2(b2L,j)*vel(i,j+b2L,k,3)
!pgi$ unroll = n:8
           DO jj = b2L+1, b2U
              nlbc32(i,k,2) = nlbc32(i,k,2) + cp2(jj,j)*vel(i,j+jj,k,3)
           END DO
           
           dd1 = cIpw(g3L,k)*drift2(i,k+g3L,2)
!pgi$ unroll = n:8
           DO kk = g3L+1, g3U
              dd1 = dd1 + cIpw(kk,k)*drift2(i,k+kk,2)
           END DO
           nlbc32(i,k,2) = dd1*nlbc32(i,k,2)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !--- Normalkomponente --------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(3,1,3) .AND. (BC_3L == 1 .OR. BC_3L == 2)) THEN
     k = 0
     DO j = S23B, N23B
        DO i = S13B, N13B
           nlbc33(i,j,1) = cDw3(d3L,1)*vel(i,j,1+d3L,3)
!pgi$ unroll = n:8
           DO kk = d3L+1, d3U
              nlbc33(i,j,1) = nlbc33(i,j,1) + cDw3(kk,1)*vel(i,j,1+kk,3)
           END DO
           nlbc33(i,j,1) = drift3(i,j,1)*nlbc33(i,j,1)
        END DO
     END DO
  END IF
  IF (outlet(3,2,3) .AND. (BC_3U == 1 .OR. BC_3U == 2)) THEN
     k = N3
     DO j = S23B, N23B
        DO i = S13B, N13B
           nlbc33(i,j,2) = cDw3(d3L,k)*vel(i,j,k+d3L,3)
!pgi$ unroll = n:8
           DO kk = d3L+1, d3U
              nlbc33(i,j,2) = nlbc33(i,j,2) + cDw3(kk,k)*vel(i,j,k+kk,3)
           END DO
           nlbc33(i,j,2) = drift3(i,j,2)*nlbc33(i,j,2)
        END DO
     END DO
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE outflow_bc
  
  
  
  
  
  
  
  
  
  
  !> subroutine that extrapolates boundary conditions for wall-normal velocities onto the boundaries.
  SUBROUTINE bc_extrapolation(m,phi)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  m                                           !< direction
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< ?
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Achtung: - Wird in "product_div_grad" und "explicit" verwendet.                                          !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! TEST!!! bislang nur Dirichlet-RB eingebaut!
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 1) THEN
     IF (BC_1L .GT. 0) THEN
        i = 0
        DO k = S31B, N31B
           DO j = S21B, N21B
!pgi$ unroll = n:8
              DO ii = 0, d1U
                 phi(i,j,k) = phi(i,j,k) - cIup(ii,1)*phi(1+ii,j,k)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIup(-1,1)
           END DO
        END DO
     END IF
     IF (BC_1U .GT. 0) THEN
        i = N1
        DO k = S31B, N31B
           DO j = S21B, N21B
!pgi$ unroll = n:8
              DO ii = d1L, -1
                 phi(i,j,k) = phi(i,j,k) - cIup(ii,i)*phi(i+ii,j,k)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIup(0,i)
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN
     IF (BC_2L .GT. 0) THEN
        j = 0
        DO k = S32B, N32B
           DO i = S12B, N12B
!pgi$ unroll = n:8
              DO jj = 0, d2U
                 phi(i,j,k) = phi(i,j,k) - cIvp(jj,1)*phi(i,1+jj,k)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIvp(-1,1)
           END DO
        END DO
     END IF
     IF (BC_2U .GT. 0) THEN
        j = N2
        DO k = S32B, N32B
           DO i = S12B, N12B
!pgi$ unroll = n:8
              DO jj = d2L, -1
                 phi(i,j,k) = phi(i,j,k) - cIvp(jj,j)*phi(i,j+jj,k)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIvp(0,j)
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN
     IF (BC_3L .GT. 0) THEN
        k = 0
        DO j = S23B, N23B
           DO i = S13B, N13B
!pgi$ unroll = n:8
              DO kk = 0, d3U
                 phi(i,j,k) = phi(i,j,k) - cIwp(kk,1)*phi(i,j,1+kk)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIwp(-1,1)
           END DO
        END DO
     END IF
     IF (BC_3U .GT. 0) THEN
        k = N3
        DO j = S23B, N23B
           DO i = S13B, N13B
!pgi$ unroll = n:8
              DO kk = d3L, -1
                 phi(i,j,k) = phi(i,j,k) - cIwp(kk,k)*phi(i,j,k+kk)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIwp(0,k)
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  
  
  END SUBROUTINE bc_extrapolation
  
  
  
  
  
  
  
  
  
  
  !> extrapolates the boundary conditions for wall-normal velocities to the boundaries for the transpose 
  !! Laplacian
  !! @note: only Dirichlet BCs have been implemented.
  SUBROUTINE bc_extrapolation_transp(m,phi)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  m
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Achtung: - Wird in "product_div_grad" und "explicit" verwendet.                                          !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! TEST!!! bislang nur Dirichlet-RB eingebaut!
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 1) THEN
     IF (BC_1L .GT. 0) THEN
        i = 0
        DO k = S31B, N31B
           DO j = S21B, N21B
!pgi$ unroll = n:8
              DO ii = 0, d1U
                 phi(i+ii+1,j,k) = phi(i+ii+1,j,k) - phi(i,j,k)*cIup(ii,1)/cIup(-1,1)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIup(-1,1)
           END DO
        END DO
     END IF
     IF (BC_1U .GT. 0) THEN
        i = N1
        DO k = S31B, N31B
           DO j = S21B, N21B
!pgi$ unroll = n:8
              DO ii = d1L, -1
                 phi(i+ii,j,k) = phi(i+ii,j,k) - phi(i,j,k)*cIup(ii,i)/cIup(0,i)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIup(0,i)
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN
     IF (BC_2L .GT. 0) THEN
        j = 0
        DO k = S32B, N32B
           DO i = S12B, N12B
!pgi$ unroll = n:8
              DO jj = 0, d2U
                 phi(i,j+jj+1,k) = phi(i,j+jj+1,k) - phi(i,j,k)*cIvp(jj,1)/cIvp(-1,1)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIvp(-1,1)
           END DO
        END DO
     END IF
     IF (BC_2U .GT. 0) THEN
        j = N2
        DO k = S32B, N32B
           DO i = S12B, N12B
!pgi$ unroll = n:8
              DO jj = d2L, -1
                 phi(i,j+jj,k) = phi(i,j+jj,k) - phi(i,j,k)*cIvp(jj,j)/cIvp(0,j)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIvp(0,j)
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN
     IF (BC_3L .GT. 0) THEN
        k = 0
        DO j = S23B, N23B
           DO i = S13B, N13B
!pgi$ unroll = n:8
              DO kk = 0, d3U
                 phi(i,j,k+kk+1) = phi(i,j,k+kk+1) - phi(i,j,k)*cIwp(kk,1)/cIwp(-1,1)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIwp(-1,1)
           END DO
        END DO
     END IF
     IF (BC_3U .GT. 0) THEN
        k = N3
        DO j = S23B, N23B
           DO i = S13B, N13B
!pgi$ unroll = n:8
              DO kk = d3L, -1
                 phi(i,j,k+kk) = phi(i,j,k+kk) - phi(i,j,k)*cIwp(kk,k)/cIwp(0,k)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIwp(0,k)
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  
  
  END SUBROUTINE bc_extrapolation_transp
  
  
  
END MODULE mod_diff
