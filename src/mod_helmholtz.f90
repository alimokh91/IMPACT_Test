!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!*************************************************************************************************************

!> module containing routines for the Helmholtz problem. It uses modules mod_dims, mod_vars, mod_diff and 
!! mod_exchange.
MODULE mod_helmholtz
  
  
  USE mod_dims
  USE mod_vars
  USE mod_diff
  USE mod_exchange
  
  
  PRIVATE
  
  PUBLIC product_Helmholtz, product_Helmholtz_precond
  PUBLIC product_Helmholtz_relax, product_Helmholtz_relax_coarse
  PUBLIC relaxation_Helmholtz, relaxation_Helmholtz_coarse
  PUBLIC handle_corner_conc
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  ! TEST!!! Generell wäre es sinnvoll, für alle Geschwindigkeiten und den Druck eigene Stencils zu speichern,
  !         um die Programmierung übersichtlicher zu halten und die IF-Abfragen zu vermeiden.
  
  !> subroutine that constructs the Laplace operator. 
  SUBROUTINE product_Helmholtz(phi,Hel)
  
  IMPLICIT NONE
  
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< ?
  REAL   , INTENT(out  ) ::  Hel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< Laplace operator
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Routine "Helmholtz" (ohne RB) ist ausgelagert, da sie auch aus RHS heraus aufgerufen wird,!
  !                wo keine RB notwendig sind. Sie werden erst anschliessend aufgeprägt.                     !
  !              - Austausch und Null-Setzen am Rand wird in Helmholtz-Routine erledigt.                     !
  !              - Interpolation in Wand-normaler Richting muss bei i=0 etc. um 1 versetzt werden!           !
  !              - Randbedingungen müssen auch in Ecken und Kanten gerechnet werden, daher die Intervalle    !
  !                S11B:N11B & Co..                                                                          !
  !              - Reihenfolge der Randbedingungen muss jeweils angepasst werden, da sich Ecken und Kanten   !
  !                überschneiden! Für die lid-driven cavity sind die tangentialen Richtungen zuletzt aufzu-  !
  !                zuprägen.                                                                                 !
  !              - Randbedingungen können NICHT in die Stencils mit eingebaut werden, da hier der Laplace-   !
  !                Operator gebildet wird (anders als z.B. bei der Divergenz).                               !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (thetaL == 0.) THEN
     IF (direction == 1) Hel(S11:N11,S21:N21,S31:N31) = phi(S11:N11,S21:N21,S31:N31)
     IF (direction == 2) Hel(S12:N12,S22:N22,S32:N32) = phi(S12:N12,S22:N22,S32:N32)
     IF (direction == 3) Hel(S13:N13,S23:N23,S33:N33) = phi(S13:N13,S23:N23,S33:N33)

  ELSE
     CALL Helmholtz(direction,.TRUE.,phi,Hel)
  END IF
  
  !===========================================================================================================
  IF (direction == 1) THEN
     
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L == 1 .OR. BC_2L == 3) Hel(S11B:N11B,1 ,S31B:N31B) = phi(S11B:N11B,1 ,S31B:N31B)
     IF (BC_2U == 1 .OR. BC_2U == 3) Hel(S11B:N11B,N2,S31B:N31B) = phi(S11B:N11B,N2,S31B:N31B)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L == 2 .OR. BC_2L == 4) THEN
        j = 1
        DO k = S31B, N31B
           DO i = S11B, N11B
              Hel(i,j,k) = cp2(b2L,j)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
              DO jj = b2L+1, b2U
                 Hel(i,j,k) = Hel(i,j,k) + cp2(jj,j)*phi(i,j+jj,k)
              END DO
           END DO
        END DO
     END IF
     IF (BC_2U == 2 .OR. BC_2U == 4) THEN
        j = N2
        DO k = S31B, N31B
           DO i = S11B, N11B
              Hel(i,j,k) = cp2(b2L,j)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
              DO jj = b2L+1, b2U
                 Hel(i,j,k) = Hel(i,j,k) + cp2(jj,j)*phi(i,j+jj,k)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3L == 1 .OR. BC_3L == 3) Hel(S11B:N11B,S21B:N21B,1 ) = phi(S11B:N11B,S21B:N21B,1 )
     IF (BC_3U == 1 .OR. BC_3U == 3) Hel(S11B:N11B,S21B:N21B,N3) = phi(S11B:N11B,S21B:N21B,N3)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3L == 2 .OR. BC_3L == 4) THEN
        k = 1
        DO j = S21B, N21B
           DO i = S11B, N11B
              Hel(i,j,k) = cp3(b3L,k)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
              DO kk = b3L+1, b3U
                 Hel(i,j,k) = Hel(i,j,k) + cp3(kk,k)*phi(i,j,k+kk)
              END DO
           END DO
        END DO
     END IF
     IF (BC_3U == 2 .OR. BC_3U == 4) THEN
        k = N3
        DO j = S21B, N21B
           DO i = S11B, N11B
              Hel(i,j,k) = cp3(b3L,k)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
              DO kk = b3L+1, b3U
                 Hel(i,j,k) = Hel(i,j,k) + cp3(kk,k)*phi(i,j,k+kk)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L == 1 .OR. BC_1L == 2) THEN
        i = 0
        DO k = S31B, N31B
           DO j = S21B, N21B
              Hel(i,j,k) = cIup(d1L,1)*phi(1+d1L,j,k)
!pgi$ unroll = n:8
              DO ii = d1L+1, d1U
                 Hel(i,j,k) = Hel(i,j,k) + cIup(ii,1)*phi(1+ii,j,k)
              END DO
           END DO
        END DO
     END IF
     IF (BC_1U == 1 .OR. BC_1U == 2) THEN
        i = N1
        DO k = S31B, N31B
           DO j = S21B, N21B
              Hel(i,j,k) = cIup(d1L,i)*phi(i+d1L,j,k)
!pgi$ unroll = n:8
              DO ii = d1L+1, d1U
                 Hel(i,j,k) = Hel(i,j,k) + cIup(ii,i)*phi(i+ii,j,k)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L == 3 .OR. BC_1L == 4) THEN
        i = 0
        DO k = S31B, N31B
           DO j = S21B, N21B
              Hel(i,j,k) = cDu1(d1L,i)*phi(1+d1L,j,k)
!pgi$ unroll = n:8
              DO ii = d1L+1, d1U
                 Hel(i,j,k) = Hel(i,j,k) + cDu1(ii,i)*phi(1+ii,j,k)
              END DO
           END DO
        END DO
     END IF
     IF (BC_1U == 3 .OR. BC_1U == 4) THEN
        i = N1
        DO k = S31B, N31B
           DO j = S21B, N21B
              Hel(i,j,k) = cDu1(d1L,i)*phi(i+d1L,j,k)
!pgi$ unroll = n:8
              DO ii = d1L+1, d1U
                 Hel(i,j,k) = Hel(i,j,k) + cDu1(ii,i)*phi(i+ii,j,k)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  IF (direction == 2) THEN
     
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L == 1 .OR. BC_1L == 3) Hel(1 ,S22B:N22B,S32B:N32B) = phi(1 ,S22B:N22B,S32B:N32B)
     IF (BC_1U == 1 .OR. BC_1U == 3) Hel(N1,S22B:N22B,S32B:N32B) = phi(N1,S22B:N22B,S32B:N32B)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L == 2 .OR. BC_1L == 4) THEN
        i = 1
        DO k = S32B, N32B
           DO j = S22B, N22B
              Hel(i,j,k) = cp1(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
              DO ii = b1L+1, b1U
                 Hel(i,j,k) = Hel(i,j,k) + cp1(ii,i)*phi(i+ii,j,k)
              END DO
           END DO
        END DO
     END IF
     IF (BC_1U == 2 .OR. BC_1U == 4) THEN
        i = N1
        DO k = S32B, N32B
           DO j = S22B, N22B
              Hel(i,j,k) = cp1(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
              DO ii = b1L+1, b1U
                 Hel(i,j,k) = Hel(i,j,k) + cp1(ii,i)*phi(i+ii,j,k)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3L == 1 .OR. BC_3L == 3) Hel(S12B:N12B,S22B:N22B,1 ) = phi(S12B:N12B,S22B:N22B,1 )
     IF (BC_3U == 1 .OR. BC_3U == 3) Hel(S12B:N12B,S22B:N22B,N3) = phi(S12B:N12B,S22B:N22B,N3)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3L == 2 .OR. BC_3L == 4) THEN
        k = 1
        DO j = S22B, N22B
           DO i = S12B, N12B
              Hel(i,j,k) = cp3(b3L,k)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
              DO kk = b3L+1, b3U
                 Hel(i,j,k) = Hel(i,j,k) + cp3(kk,k)*phi(i,j,k+kk)
              END DO
           END DO
        END DO
     END IF
     IF (BC_3U == 2 .OR. BC_3U == 4) THEN
        k = N3
        DO j = S22B, N22B
           DO i = S12B, N12B
              Hel(i,j,k) = cp3(b3L,k)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
              DO kk = b3L+1, b3U
                 Hel(i,j,k) = Hel(i,j,k) + cp3(kk,k)*phi(i,j,k+kk)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L == 1 .OR. BC_2L == 2) THEN
        j = 0
        DO k = S32B, N32B
           DO i = S12B, N12B
              Hel(i,j,k) = cIvp(d2L,1)*phi(i,1+d2L,k)
!pgi$ unroll = n:8
              DO jj = d2L+1, d2U
                 Hel(i,j,k) = Hel(i,j,k) + cIvp(jj,1)*phi(i,1+jj,k)
              END DO
           END DO
        END DO
     END IF
     IF (BC_2U == 1 .OR. BC_2U == 2) THEN
        j = N2
        DO k = S32B, N32B
           DO i = S12B, N12B
              Hel(i,j,k) = cIvp(d2L,j)*phi(i,j+d2L,k)
!pgi$ unroll = n:8
              DO jj = d2L+1, d2U
                 Hel(i,j,k) = Hel(i,j,k) + cIvp(jj,j)*phi(i,j+jj,k)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L == 3 .OR. BC_2L == 4) THEN
        j = 0
        DO k = S32B, N32B
           DO i = S12B, N12B
              Hel(i,j,k) = cDv2(d2L,1)*phi(i,1+d2L,k)
!pgi$ unroll = n:8
              DO jj = d2L+1, d2U
                 Hel(i,j,k) = Hel(i,j,k) + cDv2(jj,1)*phi(i,1+jj,k)
              END DO
           END DO
        END DO
     END IF
     IF (BC_2U == 3 .OR. BC_2U == 4) THEN
        j = N2
        DO k = S32B, N32B
           DO i = S12B, N12B
              Hel(i,j,k) = cDv2(d2L,j)*phi(i,j+d2L,k)
!pgi$ unroll = n:8
              DO jj = d2L+1, d2U
                 Hel(i,j,k) = Hel(i,j,k) + cDv2(jj,j)*phi(i,j+jj,k)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  IF (direction == 3) THEN
     
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L == 1 .OR. BC_1L == 3) Hel(1 ,S23B:N23B,S33B:N33B) = phi(1 ,S23B:N23B,S33B:N33B)
     IF (BC_1U == 1 .OR. BC_1U == 3) Hel(N1,S23B:N23B,S33B:N33B) = phi(N1,S23B:N23B,S33B:N33B)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L == 2 .OR. BC_1L == 4) THEN
        i = 1
        DO k = S33B, N33B
           DO j = S23B, N23B
              Hel(i,j,k) = cp1(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
              DO ii = b1L+1, b1U
                 Hel(i,j,k) = Hel(i,j,k) + cp1(ii,i)*phi(i+ii,j,k)
              END DO
           END DO
        END DO
     END IF
     IF (BC_1U == 2 .OR. BC_1U == 4) THEN
        i = N1
        DO k = S33B, N33B
           DO j = S23B, N23B
              Hel(i,j,k) = cp1(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
              DO ii = b1L+1, b1U
                 Hel(i,j,k) = Hel(i,j,k) + cp1(ii,i)*phi(i+ii,j,k)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L == 1 .OR. BC_2L == 3) Hel(S13B:N13B,1 ,S33B:N33B) = phi(S13B:N13B,1 ,S33B:N33B)
     IF (BC_2U == 1 .OR. BC_2U == 3) Hel(S13B:N13B,N2,S33B:N33B) = phi(S13B:N13B,N2,S33B:N33B)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L == 2 .OR. BC_2L == 4) THEN
        j = 1
        DO k = S33B, N33B
           DO i = S13B, N13B
              Hel(i,j,k) = cp2(b2L,j)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
              DO jj = b2L+1, b2U
                 Hel(i,j,k) = Hel(i,j,k) + cp2(jj,j)*phi(i,j+jj,k)
              END DO
           END DO
        END DO
     END IF
     IF (BC_2U == 2 .OR. BC_2U == 4) THEN
        j = N2
        DO k = S33B, N33B
           DO i = S13B, N13B
              Hel(i,j,k) = cp2(b2L,j)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
              DO jj = b2L+1, b2U
                 Hel(i,j,k) = Hel(i,j,k) + cp2(jj,j)*phi(i,j+jj,k)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3L == 1 .OR. BC_3L == 2) THEN
        k = 0
        DO j = S23B, N23B
           DO i = S13B, N13B
              Hel(i,j,k) = cIwp(d3L,1)*phi(i,j,1+d3L)
!pgi$ unroll = n:8
              DO kk = d3L+1, d3U
                 Hel(i,j,k) = Hel(i,j,k) + cIwp(kk,1)*phi(i,j,1+kk)
              END DO
           END DO
        END DO
     END IF
     IF (BC_3U == 1 .OR. BC_3U == 2) THEN
        k = N3
        DO j = S23B, N23B
           DO i = S13B, N13B
              Hel(i,j,k) = cIwp(d3L,k)*phi(i,j,k+d3L)
!pgi$ unroll = n:8
              DO kk = d3L+1, d3U
                 Hel(i,j,k) = Hel(i,j,k) + cIwp(kk,k)*phi(i,j,k+kk)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3L == 3 .OR. BC_3L == 4) THEN
        k = 0
        DO j = S23B, N23B
           DO i = S13B, N13B
              Hel(i,j,k) = cDw3(d3L,1)*phi(i,j,1+d3L)
!pgi$ unroll = n:8
              DO kk = d3L+1, d3U
                 Hel(i,j,k) = Hel(i,j,k) + cDw3(kk,1)*phi(i,j,1+kk)
              END DO
           END DO
        END DO
     END IF
     IF (BC_3U == 3 .OR. BC_3U == 4) THEN
        k = N3
        DO j = S23B, N23B
           DO i = S13B, N13B
              Hel(i,j,k) = cDw3(d3L,k)*phi(i,j,k+d3L)
!pgi$ unroll = n:8
              DO kk = d3L+1, d3U
                 Hel(i,j,k) = Hel(i,j,k) + cDw3(kk,k)*phi(i,j,k+kk)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE product_Helmholtz






  
  
  
  
  
  SUBROUTINE product_Helmholtz_relax(phi,Hel) ! TEST!!! 2D-Variante fehlt noch ...
  
  IMPLICIT NONE
  
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(out  ) ::  Hel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, j, k, g
  
  
  !-----------------------------------------------------------------------------------------------------!
  ! Anmerkung: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gewählt sind!        !
  !            - Reihenfolge der Randbedingungen analog zu product_Helmholtz!                           !
  !-----------------------------------------------------------------------------------------------------!
  
  ! Alternative: hier direkt auf Druckgitter differenzieren und gleiche Restriktion wie für Druck verwenden!
  
  g = 1
  
  CALL exchange_relax(g,0,0,0,direction,.TRUE.,phi)
  
  !===========================================================================================================
  IF (direction == 1) THEN
     
     !--------------------------------------------------------------------------------------------------------
     !--- Feld -----------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
        DO k = S31, N31
           DO j = S21, N21
!pgi$ unroll = n:8
              DO i = S11, N11
                 Hel(i,j,k) =  phi(i,j,k) - multL*(                                     &
                         &     cu11R(-1,i,g)*phi(i-1,j,k) + cu11R(1,i,g)*phi(i+1,j,k)   &
                         &  +  cp22R(-1,j,g)*phi(i,j-1,k) + cp22R(1,j,g)*phi(i,j+1,k)   &
                         &  +  cp33R(-1,k,g)*phi(i,j,k-1) + cp33R(1,k,g)*phi(i,j,k+1)   &
                         &  + (cu11R(0,i,g) + cp22R(0,j,g) + cp33R(0,k,g))*phi(i,j,k))
              END DO
           END DO
        END DO
     ELSE
        DO k = S31, N31
           DO j = S21, N21
!pgi$ unroll = n:8
              DO i = S11, N11
                 Hel(i,j,k) =  phi(i,j,k) - multL*(                                     &
                         &     cu11R(-1,i,g)*phi(i-1,j,k) + cu11R(1,i,g)*phi(i+1,j,k)   &
                         &  +  cp22R(-1,j,g)*phi(i,j-1,k) + cp22R(1,j,g)*phi(i,j+1,k)   &
                         &  + (cu11R(0,i,g) + cp22R(0,j,g))*phi(i,j,k))
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L == 1 .OR. BC_2L == 3) Hel(S11B:N11B,1 ,S31B:N31B) = phi(S11B:N11B,1 ,S31B:N31B)
     IF (BC_2U == 1 .OR. BC_2U == 3) Hel(S11B:N11B,N2,S31B:N31B) = phi(S11B:N11B,N2,S31B:N31B)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L == 2 .OR. BC_2L == 4) THEN
        j = 1
        DO k = S31B, N31B
!pgi$ unroll = n:8
           DO i = S11B, N11B
              Hel(i,j,k) = cp22R(0,j,g)*phi(i,j,k) + cp22R( 1,j,g)*phi(i,j+1,k)
           END DO
        END DO
     END IF
     IF (BC_2U == 2 .OR. BC_2U == 4) THEN
        j = N2
        DO k = S31B, N31B
!pgi$ unroll = n:8
           DO i = S11B, N11B
              Hel(i,j,k) = cp22R(0,j,g)*phi(i,j,k) + cp22R(-1,j,g)*phi(i,j-1,k)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3L == 1 .OR. BC_3L == 3) Hel(S11B:N11B,S21B:N21B,1 ) = phi(S11B:N11B,S21B:N21B,1 )
     IF (BC_3U == 1 .OR. BC_3U == 3) Hel(S11B:N11B,S21B:N21B,N3) = phi(S11B:N11B,S21B:N21B,N3)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3L == 2 .OR. BC_3L == 4) THEN
        k = 1
        DO j = S21B, N21B
!pgi$ unroll = n:8
           DO i = S11B, N11B
              Hel(i,j,k) = cp33R(0,k,g)*phi(i,j,k) + cp33R( 1,k,g)*phi(i,j,k+1)
           END DO
        END DO
     END IF
     IF (BC_3U == 2 .OR. BC_3U == 4) THEN
        k = N3
        DO j = S21B, N21B
!pgi$ unroll = n:8
           DO i = S11B, N11B
              Hel(i,j,k) = cp33R(0,k,g)*phi(i,j,k) + cp33R(-1,k,g)*phi(i,j,k-1)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) THEN
        i = 0
        DO k = S31B, N31B
!pgi$ unroll = n:8
           DO j = S21B, N21B
              Hel(i,j,k) = cu11R(0,i,g)*phi(i,j,k) + cu11R( 1,i,g)*phi(i+1,j,k)
           END DO
        END DO
     END IF
     IF (BC_1U .GT. 0) THEN
        i = N1
        DO k = S31B, N31B
!pgi$ unroll = n:8
           DO j = S21B, N21B
              Hel(i,j,k) = cu11R(0,i,g)*phi(i,j,k) + cu11R(-1,i,g)*phi(i-1,j,k)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  IF (direction == 2) THEN
     
     !--------------------------------------------------------------------------------------------------------
     !--- Feld -----------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
        DO k = S32, N32
           DO j = S22, N22
!pgi$ unroll = n:8
              DO i = S12, N12
                 Hel(i,j,k) =  phi(i,j,k) - multL*(                                     &
                         &     cp11R(-1,i,g)*phi(i-1,j,k) + cp11R(1,i,g)*phi(i+1,j,k)   &
                         &  +  cv22R(-1,j,g)*phi(i,j-1,k) + cv22R(1,j,g)*phi(i,j+1,k)   &
                         &  +  cp33R(-1,k,g)*phi(i,j,k-1) + cp33R(1,k,g)*phi(i,j,k+1)   &
                         &  + (cp11R(0,i,g) + cv22R(0,j,g) + cp33R(0,k,g))*phi(i,j,k))
              END DO
           END DO
        END DO
     ELSE
        DO k = S32, N32
           DO j = S22, N22
!pgi$ unroll = n:8
              DO i = S12, N12
                 Hel(i,j,k) =  phi(i,j,k) - multL*(                                     &
                         &     cp11R(-1,i,g)*phi(i-1,j,k) + cp11R(1,i,g)*phi(i+1,j,k)   &
                         &  +  cv22R(-1,j,g)*phi(i,j-1,k) + cv22R(1,j,g)*phi(i,j+1,k)   &
                         &  + (cp11R(0,i,g) + cv22R(0,j,g))*phi(i,j,k))
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L == 1 .OR. BC_1L == 3) Hel(1 ,S22B:N22B,S32B:N32B) = phi(1 ,S22B:N22B,S32B:N32B)
     IF (BC_1U == 1 .OR. BC_1U == 3) Hel(N1,S22B:N22B,S32B:N32B) = phi(N1,S22B:N22B,S32B:N32B)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L == 2 .OR. BC_1L == 4) THEN
        i = 1
        DO k = S32B, N32B
!pgi$ unroll = n:8
           DO j = S22B, N22B
              Hel(i,j,k) = cp11R(0,i,g)*phi(i,j,k) + cp11R( 1,i,g)*phi(i+1,j,k)
           END DO
        END DO
     END IF
     IF (BC_1U == 2 .OR. BC_1U == 4) THEN
        i = N1
        DO k = S32B, N32B
!pgi$ unroll = n:8
           DO j = S22B, N22B
              Hel(i,j,k) = cp11R(0,i,g)*phi(i,j,k) + cp11R(-1,i,g)*phi(i-1,j,k)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3L == 1 .OR. BC_3L == 3) Hel(S12B:N12B,S22B:N22B,1 ) = phi(S12B:N12B,S22B:N22B,1 )
     IF (BC_3U == 1 .OR. BC_3U == 3) Hel(S12B:N12B,S22B:N22B,N3) = phi(S12B:N12B,S22B:N22B,N3)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3L == 2 .OR. BC_3L == 4) THEN
        k = 1
        DO j = S22B, N22B
!pgi$ unroll = n:8
           DO i = S12B, N12B
              Hel(i,j,k) = cp33R(0,k,g)*phi(i,j,k) + cp33R( 1,k,g)*phi(i,j,k+1)
           END DO
        END DO
     END IF
     IF (BC_3U == 2 .OR. BC_3U == 4) THEN
        k = N3
        DO j = S22B, N22B
!pgi$ unroll = n:8
           DO i = S12B, N12B
              Hel(i,j,k) = cp33R(0,k,g)*phi(i,j,k) + cp33R(-1,k,g)*phi(i,j,k-1)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L .GT. 0) THEN
        j = 0
        DO k = S32B, N32B
!pgi$ unroll = n:8
           DO i = S12B, N12B
              Hel(i,j,k) = cv22R(0,j,g)*phi(i,j,k) + cv22R( 1,j,g)*phi(i,j+1,k)
           END DO
        END DO
     END IF
     IF (BC_2U .GT. 0) THEN
        j = N2
        DO k = S32B, N32B
!pgi$ unroll = n:8
           DO i = S12B, N12B
              Hel(i,j,k) = cv22R(0,j,g)*phi(i,j,k) + cv22R(-1,j,g)*phi(i,j-1,k)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  IF (direction == 3) THEN
     
     !--------------------------------------------------------------------------------------------------------
     !--- Feld -----------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     ! dimens /= 3 durch direction /= 3 trivial erfüllt.
     DO k = S33, N33
        DO j = S23, N23
!pgi$ unroll = n:8
           DO i = S13, N13
              Hel(i,j,k) =  phi(i,j,k) - multL*(                                     &
                      &     cp11R(-1,i,g)*phi(i-1,j,k) + cp11R(1,i,g)*phi(i+1,j,k)   &
                      &  +  cp22R(-1,j,g)*phi(i,j-1,k) + cp22R(1,j,g)*phi(i,j+1,k)   &
                      &  +  cw33R(-1,k,g)*phi(i,j,k-1) + cw33R(1,k,g)*phi(i,j,k+1)   &
                      &  + (cp11R(0,i,g) + cp22R(0,j,g) + cw33R(0,k,g))*phi(i,j,k))
           END DO
        END DO
     END DO
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L == 1 .OR. BC_1L == 3) Hel(1 ,S23B:N23B,S33B:N33B) = phi(1 ,S23B:N23B,S33B:N33B)
     IF (BC_1U == 1 .OR. BC_1U == 3) Hel(N1,S23B:N23B,S33B:N33B) = phi(N1,S23B:N23B,S33B:N33B)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L == 2 .OR. BC_1L == 4) THEN
        i = 1
        DO k = S33B, N33B
!pgi$ unroll = n:8
           DO j = S23B, N23B
              Hel(i,j,k) = cp11R(0,i,g)*phi(i,j,k) + cp11R( 1,i,g)*phi(i+1,j,k)
           END DO
        END DO
     END IF
     IF (BC_1U == 2 .OR. BC_1U == 4) THEN
        i = N1
        DO k = S33B, N33B
!pgi$ unroll = n:8
           DO j = S23B, N23B
              Hel(i,j,k) = cp11R(0,i,g)*phi(i,j,k) + cp11R(-1,i,g)*phi(i-1,j,k)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L == 1 .OR. BC_2L == 3) Hel(S13B:N13B,1 ,S33B:N33B) = phi(S13B:N13B,1 ,S33B:N33B)
     IF (BC_2U == 1 .OR. BC_2U == 3) Hel(S13B:N13B,N2,S33B:N33B) = phi(S13B:N13B,N2,S33B:N33B)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L == 2 .OR. BC_2L == 4) THEN
        j = 1
        DO k = S33B, N33B
!pgi$ unroll = n:8
           DO i = S13B, N13B
              Hel(i,j,k) = cp22R(0,j,g)*phi(i,j,k) + cp22R( 1,j,g)*phi(i,j+1,k)
           END DO
        END DO
     END IF
     IF (BC_2U == 2 .OR. BC_2U == 4) THEN
        j = N2
        DO k = S33B, N33B
!pgi$ unroll = n:8
           DO i = S13B, N13B
              Hel(i,j,k) = cp22R(0,j,g)*phi(i,j,k) + cp22R(-1,j,g)*phi(i,j-1,k)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3L .GT. 0) THEN
        k = 0
        DO j = S23B, N23B
!pgi$ unroll = n:8
           DO i = S13B, N13B
              Hel(i,j,k) = cw33R(0,k,g)*phi(i,j,k) + cw33R( 1,k,g)*phi(i,j,k+1)
           END DO
        END DO
     END IF
     IF (BC_3U .GT. 0) THEN
        k = N3
        DO j = S23B, N23B
!pgi$ unroll = n:8
           DO i = S13B, N13B
              Hel(i,j,k) = cw33R(0,k,g)*phi(i,j,k) + cw33R(-1,k,g)*phi(i,j,k-1)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE product_Helmholtz_relax
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE product_Helmholtz_relax_coarse(g,phi,Hel) ! TEST!!! 2D-Variante fehlt noch ...
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  g
  
  REAL   , INTENT(inout) ::  phi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL   , INTENT(out  ) ::  Hel(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  INTEGER                ::  i, N1, N1R, N11R
  INTEGER                ::  j, N2, N2R, N22R
  INTEGER                ::  k, N3, N3R, N33R
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkung: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gewählt sind!             !
  !----------------------------------------------------------------------------------------------------------!
  
  N1 = NN(1,g)
  N2 = NN(2,g)
  N3 = NN(3,g)
  
  N1R  = N1 + d1R ! TEST!!! Substituieren ...
  N2R  = N2 + d2R
  N3R  = N3 + d3R
  
  N11R = N1 + d11R
  N22R = N2 + d22R
  N33R = N3 + d33R
  
  
  CALL exchange_relax(g,0,0,0,0,.TRUE.,phi)
  
  
  !===========================================================================================================
  IF (direction == 1) THEN
     
     !--------------------------------------------------------------------------------------------------------
     !--- Feld -----------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
        DO k = S33R, N33R
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 Hel(i,j,k) = phi(i,j,k) - multL*(                                       &
                         &     cu11R(-1,i,g)*phi(i-1,j,k) + cu11R(1,i,g)*phi(i+1,j,k)    &
                         &   + cp22R(-1,j,g)*phi(i,j-1,k) + cp22R(1,j,g)*phi(i,j+1,k)    &
                         &   + cp33R(-1,k,g)*phi(i,j,k-1) + cp33R(1,k,g)*phi(i,j,k+1)    &
                         &   + phi(i,j,k)*(cu11R(0,i,g) + cp22R(0,j,g) + cp33R(0,k,g)))
              END DO
           END DO
        END DO
     ELSE
        DO k = S33R, N33R
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 Hel(i,j,k) = phi(i,j,k) - multL*(                                       &
                         &     cu11R(-1,i,g)*phi(i-1,j,k) + cu11R(1,i,g)*phi(i+1,j,k)    &
                         &   + cp22R(-1,j,g)*phi(i,j-1,k) + cp22R(1,j,g)*phi(i,j+1,k)    &
                         &   + phi(i,j,k)*(cu11R(0,i,g) + cp22R(0,j,g)))
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L == 1 .OR. BC_2L == 3) Hel(S1R:N1R,1 ,S3R:N3R) = phi(S1R:N1R,1 ,S3R:N3R)
     IF (BC_2U == 1 .OR. BC_2U == 3) Hel(S1R:N1R,N2,S3R:N3R) = phi(S1R:N1R,N2,S3R:N3R)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L == 2 .OR. BC_2L == 4) THEN
        j = 1
        DO k = S3R, N3R
!pgi$ unroll = n:8
           DO i = S1R, N1R
              Hel(i,j,k) = cp22R(0,j,g)*phi(i,j,k) + cp22R( 1,j,g)*phi(i,j+1,k)
           END DO
        END DO
     END IF
     IF (BC_2U == 2 .OR. BC_2U == 4) THEN
        j = N2
        DO k = S3R, N3R
!pgi$ unroll = n:8
           DO i = S1R, N1R
              Hel(i,j,k) = cp22R(0,j,g)*phi(i,j,k) + cp22R(-1,j,g)*phi(i,j-1,k)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3L == 1 .OR. BC_3L == 3) Hel(S1R:N1R,S2R:N2R,1 ) = phi(S1R:N1R,S2R:N2R,1 )
     IF (BC_3U == 1 .OR. BC_3U == 3) Hel(S1R:N1R,S2R:N2R,N3) = phi(S1R:N1R,S2R:N2R,N3)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3L == 2 .OR. BC_3L == 4) THEN
        k = 1
        DO j = S2R, N2R
!pgi$ unroll = n:8
           DO i = S1R, N1R
              Hel(i,j,k) = cp33R(0,k,g)*phi(i,j,k) + cp33R( 1,k,g)*phi(i,j,k+1)
           END DO
        END DO
     END IF
     IF (BC_3U == 2 .OR. BC_3U == 4) THEN
        k = N3
        DO j = S2R, N2R
!pgi$ unroll = n:8
           DO i = S1R, N1R
              Hel(i,j,k) = cp33R(0,k,g)*phi(i,j,k) + cp33R(-1,k,g)*phi(i,j,k-1)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) THEN
        i = 1
        DO k = S3R, N3R
!pgi$ unroll = n:8
           DO j = S2R, N2R
              Hel(i,j,k) = cu11R(0,i,g)*phi(i,j,k) + cu11R( 1,i,g)*phi(i+1,j,k)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1U .GT. 0) THEN
        i = N1
        DO k = S3R, N3R
!pgi$ unroll = n:8
           DO j = S2R, N2R
              Hel(i,j,k) = cu11R(0,i,g)*phi(i,j,k) + cu11R(-1,i,g)*phi(i-1,j,k)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  IF (direction == 2) THEN
     
     !--------------------------------------------------------------------------------------------------------
     !--- Feld -----------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
        DO k = S33R, N33R
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 Hel(i,j,k) = phi(i,j,k) - multL*(                                       &
                         &     cp11R(-1,i,g)*phi(i-1,j,k) + cp11R(1,i,g)*phi(i+1,j,k)    &
                         &   + cv22R(-1,j,g)*phi(i,j-1,k) + cv22R(1,j,g)*phi(i,j+1,k)    &
                         &   + cp33R(-1,k,g)*phi(i,j,k-1) + cp33R(1,k,g)*phi(i,j,k+1)    &
                         &   + phi(i,j,k)*(cp11R(0,i,g) + cv22R(0,j,g) + cp33R(0,k,g)))
              END DO
           END DO
        END DO
     ELSE
        DO k = S33R, N33R
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 Hel(i,j,k) = phi(i,j,k) - multL*(                                       &
                         &     cp11R(-1,i,g)*phi(i-1,j,k) + cp11R(1,i,g)*phi(i+1,j,k)    &
                         &   + cv22R(-1,j,g)*phi(i,j-1,k) + cv22R(1,j,g)*phi(i,j+1,k)    &
                         &   + phi(i,j,k)*(cp11R(0,i,g) + cv22R(0,j,g)))
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L == 1 .OR. BC_1L == 3) Hel(1 ,S2R:N2R,S3R:N3R) = phi(1 ,S2R:N2R,S3R:N3R)
     IF (BC_1U == 1 .OR. BC_1U == 3) Hel(N1,S2R:N2R,S3R:N3R) = phi(N1,S2R:N2R,S3R:N3R)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L == 2 .OR. BC_1L == 4) THEN
        i = 1
        DO k = S3R, N3R
!pgi$ unroll = n:8
           DO j = S2R, N2R
              Hel(i,j,k) = cp11R(0,i,g)*phi(i,j,k) + cp11R( 1,i,g)*phi(i+1,j,k)
           END DO
        END DO
     END IF
     IF (BC_1U == 2 .OR. BC_1U == 4) THEN
        i = N1
        DO k = S3R, N3R
!pgi$ unroll = n:8
           DO j = S2R, N2R
              Hel(i,j,k) = cp11R(0,i,g)*phi(i,j,k) + cp11R(-1,i,g)*phi(i-1,j,k)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3L == 1 .OR. BC_3L == 3) Hel(S1R:N1R,S2R:N2R,1 ) = phi(S1R:N1R,S2R:N2R,1 )
     IF (BC_3U == 1 .OR. BC_3U == 3) Hel(S1R:N1R,S2R:N2R,N3) = phi(S1R:N1R,S2R:N2R,N3)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3L == 2 .OR. BC_3L == 4) THEN
        k = 1
        DO j = S2R, N2R
!pgi$ unroll = n:8
           DO i = S1R, N1R
              Hel(i,j,k) = cp33R(0,k,g)*phi(i,j,k) + cp33R( 1,k,g)*phi(i,j,k+1)
           END DO
        END DO
     END IF
     IF (BC_3U == 2 .OR. BC_3U == 4) THEN
        k = N3
        DO j = S2R, N2R
!pgi$ unroll = n:8
           DO i = S1R, N1R
              Hel(i,j,k) = cp33R(0,k,g)*phi(i,j,k) + cp33R(-1,k,g)*phi(i,j,k-1)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L .GT. 0) THEN
        j = 1
        DO k = S3R, N3R
!pgi$ unroll = n:8
           DO i = S1R, N1R
              Hel(i,j,k) = cv22R(0,j,g)*phi(i,j,k) + cv22R( 1,j,g)*phi(i,j+1,k)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2U .GT. 0) THEN
        j = N2
        DO k = S3R, N3R
!pgi$ unroll = n:8
           DO i = S1R, N1R
              Hel(i,j,k) = cv22R(0,j,g)*phi(i,j,k) + cv22R(-1,j,g)*phi(i,j-1,k)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  ! dimens /= 3 durch direction /= 3 trivial erfüllt.
  IF (direction == 3) THEN
     DO k = S33R, N33R
        DO j = S22R, N22R
!pgi$ unroll = n:8
           DO i = S11R, N11R
              Hel(i,j,k) = phi(i,j,k) - multL*(                                       &
                      &     cp11R(-1,i,g)*phi(i-1,j,k) + cp11R(1,i,g)*phi(i+1,j,k)    &
                      &   + cp22R(-1,j,g)*phi(i,j-1,k) + cp22R(1,j,g)*phi(i,j+1,k)    &
                      &   + cw33R(-1,k,g)*phi(i,j,k-1) + cw33R(1,k,g)*phi(i,j,k+1)    &
                      &   + phi(i,j,k)*(cp11R(0,i,g) + cp22R(0,j,g) + cw33R(0,k,g)))
           END DO
        END DO
     END DO
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L == 1 .OR. BC_1L == 3) Hel(1 ,S2R:N2R,S3R:N3R) = phi(1 ,S2R:N2R,S3R:N3R)
     IF (BC_1U == 1 .OR. BC_1U == 3) Hel(N1,S2R:N2R,S3R:N3R) = phi(N1,S2R:N2R,S3R:N3R)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_1L == 2 .OR. BC_1L == 4) THEN
        i = 1
        DO k = S3R, N3R
!pgi$ unroll = n:8
           DO j = S2R, N2R
              Hel(i,j,k) = cp11R(0,i,g)*phi(i,j,k) + cp11R( 1,i,g)*phi(i+1,j,k)
           END DO
        END DO
     END IF
     IF (BC_1U == 2 .OR. BC_1U == 4) THEN
        i = N1
        DO k = S3R, N3R
!pgi$ unroll = n:8
           DO j = S2R, N2R
              Hel(i,j,k) = cp11R(0,i,g)*phi(i,j,k) + cp11R(-1,i,g)*phi(i-1,j,k)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L == 1 .OR. BC_2L == 3) Hel(S1R:N1R,1 ,S3R:N3R) = phi(S1R:N1R,1 ,S3R:N3R)
     IF (BC_2U == 1 .OR. BC_2U == 3) Hel(S1R:N1R,N2,S3R:N3R) = phi(S1R:N1R,N2,S3R:N3R)
     !--------------------------------------------------------------------------------------------------------
     IF (BC_2L == 2 .OR. BC_2L == 4) THEN
        j = 1
        DO k = S3R, N3R
!pgi$ unroll = n:8
           DO i = S1R, N1R
              Hel(i,j,k) = cp22R(0,j,g)*phi(i,j,k) + cp22R( 1,j,g)*phi(i,j+1,k)
           END DO
        END DO
     END IF
     IF (BC_2U == 2 .OR. BC_2U == 4) THEN
        j = N2
        DO k = S3R, N3R
!pgi$ unroll = n:8
           DO i = S1R, N1R
              Hel(i,j,k) = cp22R(0,j,g)*phi(i,j,k) + cp22R(-1,j,g)*phi(i,j-1,k)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3L .GT. 0) THEN
        k = 1
        DO j = S2R, N2R
!pgi$ unroll = n:8
           DO i = S1R, N1R
              Hel(i,j,k) = cw33R(0,k,g)*phi(i,j,k) + cw33R( 1,k,g)*phi(i,j,k+1)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (BC_3U .GT. 0) THEN
        k = N3
        DO j = S2R, N2R
!pgi$ unroll = n:8
           DO i = S1R, N1R
              Hel(i,j,k) = cw33R(0,k,g)*phi(i,j,k) + cw33R(-1,k,g)*phi(i,j,k-1)
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE product_Helmholtz_relax_coarse
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE relaxation_Helmholtz(init_yes,n_relax,bb,rel) ! TEST!!! 2D-Variante fehlt noch ...
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  init_yes
  INTEGER, INTENT(in   ) ::  n_relax
  
  REAL   , INTENT(in   ) ::  bb (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(inout) ::  rel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, j, k, g, r
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gewählt sind!           !
  !              - Randbedingungen müssen auch in Ecken und Kanten gerechnet werden, daher die Intervalle    !
  !                S11B:N11B & Co.                                                                           !
  !              - Initialisierung bezieht sich nur auf Feldbereich, daher die Intervalle S11:N11 usw..      !
  !              - Reihenfolge der Randbedingungen analog zu product_Helmholtz!                              !
  !              - Dirichlet-Randbedingungen sind zur Optimierung der Konvergenzrate z.T. VOR der Relaxation !
  !                eingefügt!                                                                                !
  !----------------------------------------------------------------------------------------------------------!
  
  g = 1
  
  ! Nur zum testen, falls Initialisierung unten fehlerhaft erscheint:
  !IF (init_yes) rel = 0.
  
  
  DO r = 1, n_relax
     
     IF (.NOT. (r == 1 .AND. init_yes)) CALL exchange_relax(g,0,0,0,direction,.TRUE.,rel)
     
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     IF (direction == 1) THEN
        
        !=====================================================================================================
        IF (r == 1 .AND. init_yes) THEN
           IF (BC_1L .LE. 0) rel(S11-1  ,S21:N21,S31:N31) = 0.
           IF (BC_2L .LE. 0) rel(S11:N11,S21-1  ,S31:N31) = 0.
           IF (BC_3L .LE. 0) rel(S11:N11,S21:N21,S31-1  ) = 0.
        END IF
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L == 1 .OR. BC_2L == 3) rel(S11B:N11B,1 ,S31B:N31B) = bb(S11B:N11B,1 ,S31B:N31B)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L == 2 .OR. BC_2L == 4) THEN
           j = 1
           IF (r == 1 .AND. init_yes) THEN
              DO k = S31B, N31B
!pgi$ unroll = n:8
                 DO i = S11B, N11B
                    rel(i,j,k) = bb(i,j,k) / cp22R(0,j,g)
                 END DO
              END DO
           ELSE
              DO k = S31B, N31B
!pgi$ unroll = n:8
                 DO i = S11B, N11B
                    rel(i,j,k) = (bb(i,j,k) - cp22R( 1,j,g)*rel(i,j+1,k)) / cp22R(0,j,g)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L == 1 .OR. BC_3L == 3) rel(S11B:N11B,S21B:N21B,1 ) = bb(S11B:N11B,S21B:N21B,1 )
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L == 2 .OR. BC_3L == 4) THEN
           k = 1
           IF (r == 1 .AND. init_yes) THEN
              DO j = S21B, N21B
!pgi$ unroll = n:8
                 DO i = S11B, N11B
                    rel(i,j,k) = bb(i,j,k) / cp33R(0,k,g)
                 END DO
              END DO
           ELSE
              DO j = S21B, N21B
!pgi$ unroll = n:8
                 DO i = S11B, N11B
                    rel(i,j,k) = (bb(i,j,k) - cp33R( 1,k,g)*rel(i,j,k+1)) / cp33R(0,k,g)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L .GT. 0) THEN
           i = 0
           IF (r == 1 .AND. init_yes) THEN
              DO k = S31B, N31B
!pgi$ unroll = n:8
                 DO j = S21B, N21B
                    rel(i,j,k) = bb(i,j,k) / cu11R(0,i,g)
                 END DO
              END DO
           ELSE
              DO k = S31B, N31B
!pgi$ unroll = n:8
                 DO j = S21B, N21B
                    rel(i,j,k) = (bb(i,j,k) - cu11R( 1,i,g)*rel(i+1,j,k)) / cu11R(0,i,g)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Feld --------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (r == 1 .AND. init_yes) THEN
           
           DO k = S31, N31
              DO j = S21, N21
!pgi$ unroll = n:8
                 DO i = S11, N11
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                             &
                       &    cu11R(-1,i,g)*rel(i-1,j,k) + cp22R(-1,j,g)*rel(i,j-1,k) + cp33R(-1,k,g)*rel(i,j,k-1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cu11R(0,i,g) + cp22R(0,j,g) + cp33R(0,k,g)))
                 END DO
              END DO
           END DO
           
        ELSE
           
           DO k = S31, N31
              DO j = S21, N21
!pgi$ unroll = n:8
                 DO i = S11, N11
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                &
                                &        cu11R(-1,i,g)*rel(i-1,j,k) + cu11R(1,i,g)*rel(i+1,j,k)     &
                                &      + cp22R(-1,j,g)*rel(i,j-1,k) + cp22R(1,j,g)*rel(i,j+1,k)     &
                                &      + cp33R(-1,k,g)*rel(i,j,k-1) + cp33R(1,k,g)*rel(i,j,k+1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cu11R(0,i,g) + cp22R(0,j,g) + cp33R(0,k,g)))
                 END DO
              END DO
           END DO
           
        END IF
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2U == 1 .OR. BC_2U == 3) rel(S11B:N11B,N2,S31B:N31B) = bb(S11B:N11B,N2,S31B:N31B)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2U == 2 .OR. BC_2U == 4) THEN
           j = N2
           DO k = S31B, N31B
!pgi$ unroll = n:8
              DO i = S11B, N11B
                 rel(i,j,k) = (bb(i,j,k) - cp22R(-1,j,g)*rel(i,j-1,k)) / cp22R(0,j,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U == 1 .OR. BC_3U == 3) rel(S11B:N11B,S21B:N21B,N3) = bb(S11B:N11B,S21B:N21B,N3)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U == 2 .OR. BC_3U == 4) THEN
           k = N3
           DO j = S21B, N21B
!pgi$ unroll = n:8
              DO i = S11B, N11B
                 rel(i,j,k) = (bb(i,j,k) - cp33R(-1,k,g)*rel(i,j,k-1)) / cp33R(0,k,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1U .GT. 0) THEN
           i = N1
           DO k = S31B, N31B
!pgi$ unroll = n:8
              DO j = S21B, N21B
                 rel(i,j,k) = (bb(i,j,k) - cu11R(-1,i,g)*rel(i-1,j,k)) / cu11R(0,i,g)
              END DO
           END DO
        END IF
        !=====================================================================================================
        
     END IF
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     IF (direction == 2) THEN
        
        !=====================================================================================================
        IF (r == 1 .AND. init_yes) THEN
           IF (BC_1L .LE. 0) rel(S12-1  ,S22:N22,S32:N32) = 0.
           IF (BC_2L .LE. 0) rel(S12:N12,S22-1  ,S32:N32) = 0.
           IF (BC_3L .LE. 0) rel(S12:N12,S22:N22,S32-1  ) = 0.
        END IF
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L == 1 .OR. BC_1L == 3) rel(1 ,S22B:N22B,S32B:N32B) = bb(1 ,S22B:N22B,S32B:N32B)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L == 2 .OR. BC_1L == 4) THEN
           i = 1
           IF (r == 1 .AND. init_yes) THEN
              DO k = S32B, N32B
!pgi$ unroll = n:8
                 DO j = S22B, N22B
                    rel(i,j,k) = bb(i,j,k) / cp11R(0,i,g)
                 END DO
              END DO
           ELSE
              DO k = S32B, N32B
!pgi$ unroll = n:8
                 DO j = S22B, N22B
                    rel(i,j,k) = (bb(i,j,k) - cp11R( 1,i,g)*rel(i+1,j,k)) / cp11R(0,i,g)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L == 1 .OR. BC_3L == 3) rel(S12B:N12B,S22B:N22B,1 ) = bb(S12B:N12B,S22B:N22B,1 )
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L == 2 .OR. BC_3L == 4) THEN
           k = 1
           IF (r == 1 .AND. init_yes) THEN
              DO j = S22B, N22B
!pgi$ unroll = n:8
                 DO i = S12B, N12B
                    rel(i,j,k) = bb(i,j,k) / cp33R(0,k,g)
                 END DO
              END DO
           ELSE
              DO j = S22B, N22B
!pgi$ unroll = n:8
                 DO i = S12B, N12B
                    rel(i,j,k) = (bb(i,j,k) - cp33R( 1,k,g)*rel(i,j,k+1)) / cp33R(0,k,g)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L .GT. 0) THEN
           j = 0
           IF (r == 1 .AND. init_yes) THEN
              DO k = S32B, N32B
!pgi$ unroll = n:8
                 DO i = S12B, N12B
                    rel(i,j,k) = bb(i,j,k) / cv22R(0,j,g)
                 END DO
              END DO
           ELSE
              DO k = S32B, N32B
!pgi$ unroll = n:8
                 DO i = S12B, N12B
                    rel(i,j,k) = (bb(i,j,k) - cv22R( 1,j,g)*rel(i,j+1,k)) / cv22R(0,j,g)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Feld --------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (r == 1 .AND. init_yes) THEN
           
           DO k = S32, N32
              DO j = S22, N22
!pgi$ unroll = n:8
                 DO i = S12, N12
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                             &
                       &    cp11R(-1,i,g)*rel(i-1,j,k) + cv22R(-1,j,g)*rel(i,j-1,k) + cp33R(-1,k,g)*rel(i,j,k-1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cp11R(0,i,g) + cv22R(0,j,g) + cp33R(0,k,g)))
                 END DO
              END DO
           END DO
           
        ELSE
           
           DO k = S32, N32
              DO j = S22, N22
!pgi$ unroll = n:8
                 DO i = S12, N12
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                &
                                &        cp11R(-1,i,g)*rel(i-1,j,k) + cp11R(1,i,g)*rel(i+1,j,k)     &
                                &      + cv22R(-1,j,g)*rel(i,j-1,k) + cv22R(1,j,g)*rel(i,j+1,k)     &
                                &      + cp33R(-1,k,g)*rel(i,j,k-1) + cp33R(1,k,g)*rel(i,j,k+1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cp11R(0,i,g) + cv22R(0,j,g) + cp33R(0,k,g)))
                 END DO
              END DO
           END DO
           
        END IF
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1U == 1 .OR. BC_1U == 3) rel(N1,S22B:N22B,S32B:N32B) = bb(N1,S22B:N22B,S32B:N32B)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1U == 2 .OR. BC_1U == 4) THEN
           i = N1
           DO k = S32B, N32B
!pgi$ unroll = n:8
              DO j = S22B, N22B
                 rel(i,j,k) = (bb(i,j,k) - cp11R(-1,i,g)*rel(i-1,j,k)) / cp11R(0,i,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U == 1 .OR. BC_3U == 3) rel(S12B:N12B,S22B:N22B,N3) = bb(S12B:N12B,S22B:N22B,N3)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U == 2 .OR. BC_3U == 4) THEN
           k = N3
           DO j = S22B, N22B
!pgi$ unroll = n:8
              DO i = S12B, N12B
                 rel(i,j,k) = (bb(i,j,k) - cp33R(-1,k,g)*rel(i,j,k-1)) / cp33R(0,k,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2U .GT. 0) THEN
           j = N2
           DO k = S32, N32
!pgi$ unroll = n:8
              DO i = S12, N12
                 rel(i,j,k) = (bb(i,j,k) - cv22R(-1,j,g)*rel(i,j-1,k)) / cv22R(0,j,g)
              END DO
           END DO
        END IF
        !=====================================================================================================
        
     END IF
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     IF (direction == 3) THEN
        
        !=====================================================================================================
        IF (r == 1 .AND. init_yes) THEN
           IF (BC_1L .LE. 0) rel(S13-1  ,S23:N23,S33:N33) = 0.
           IF (BC_2L .LE. 0) rel(S13:N13,S23-1  ,S33:N33) = 0.
           IF (BC_3L .LE. 0) rel(S13:N13,S23:N23,S33-1  ) = 0.
        END IF
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L == 1 .OR. BC_1L == 3) rel(1 ,S23B:N23B,S33B:N33B) = bb(1 ,S23B:N23B,S33B:N33B)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L == 2 .OR. BC_1L == 4) THEN
           i = 1
           IF (r == 1 .AND. init_yes) THEN
              DO k = S33B, N33B
!pgi$ unroll = n:8
                 DO j = S23B, N23B
                    rel(i,j,k) = bb(i,j,k) / cp11R(0,i,g)
                 END DO
              END DO
           ELSE
              DO k = S33B, N33B
!pgi$ unroll = n:8
                 DO j = S23B, N23B
                    rel(i,j,k) = (bb(i,j,k) - cp11R( 1,i,g)*rel(i+1,j,k)) / cp11R(0,i,g)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L == 1 .OR. BC_2L == 3) rel(S13B:N13B,1 ,S33B:N33B) = bb(S13B:N13B,1 ,S33B:N33B)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L == 2 .OR. BC_2L == 4) THEN
           j = 1
           IF (r == 1 .AND. init_yes) THEN
              DO k = S33B, N33B
!pgi$ unroll = n:8
                 DO i = S13B, N13B
                    rel(i,j,k) = bb(i,j,k) / cp22R(0,j,g)
                 END DO
              END DO
           ELSE
              DO k = S33B, N33B
!pgi$ unroll = n:8
                 DO i = S13B, N13B
                    rel(i,j,k) = (bb(i,j,k) - cp22R( 1,j,g)*rel(i,j+1,k)) / cp22R(0,j,g)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L .GT. 0) THEN
           k = 0
           IF (r == 1 .AND. init_yes) THEN
              DO j = S23B, N23B
!pgi$ unroll = n:8
                 DO i = S13B, N13B
                    rel(i,j,k) = bb(i,j,k) / cw33R(0,k,g)
                 END DO
              END DO
           ELSE
              DO j = S23B, N23B
!pgi$ unroll = n:8
                 DO i = S13B, N13B
                    rel(i,j,k) = (bb(i,j,k) - cw33R( 1,k,g)*rel(i,j,k+1)) / cw33R(0,k,g)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Feld --------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (r == 1 .AND. init_yes) THEN
           
           DO k = S33, N33
              DO j = S23, N23
!pgi$ unroll = n:8
                 DO i = S13, N13
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                             &
                       &    cp11R(-1,i,g)*rel(i-1,j,k) + cp22R(-1,j,g)*rel(i,j-1,k) + cw33R(-1,k,g)*rel(i,j,k-1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cp11R(0,i,g) + cp22R(0,j,g) + cw33R(0,k,g)))
                 END DO
              END DO
           END DO
           
        ELSE
           
           DO k = S33, N33
              DO j = S23, N23
!pgi$ unroll = n:8
                 DO i = S13, N13
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                &
                                &        cp11R(-1,i,g)*rel(i-1,j,k) + cp11R(1,i,g)*rel(i+1,j,k)     &
                                &      + cp22R(-1,j,g)*rel(i,j-1,k) + cp22R(1,j,g)*rel(i,j+1,k)     &
                                &      + cw33R(-1,k,g)*rel(i,j,k-1) + cw33R(1,k,g)*rel(i,j,k+1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cp11R(0,i,g) + cp22R(0,j,g) + cw33R(0,k,g)))
                 END DO
              END DO
           END DO
           
        END IF
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1U == 1 .OR. BC_1U == 3) rel(N1,S23B:N23B,S33B:N33B) = bb(N1,S23B:N23B,S33B:N33B)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1U == 2 .OR. BC_1U == 4) THEN
           i = N1
           DO k = S33B, N33B
!pgi$ unroll = n:8
              DO j = S23B, N23B
                 rel(i,j,k) = (bb(i,j,k) - cp11R(-1,i,g)*rel(i-1,j,k)) / cp11R(0,i,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2U == 1 .OR. BC_2U == 3) rel(S13B:N13B,N2,S33B:N33B) = bb(S13B:N13B,N2,S33B:N33B)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2U == 2 .OR. BC_2U == 4) THEN
           j = N2
           DO k = S33B, N33B
!pgi$ unroll = n:8
              DO i = S13B, N13B
                 rel(i,j,k) = (bb(i,j,k) - cp22R(-1,j,g)*rel(i,j-1,k)) / cp22R(0,j,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U .GT. 0) THEN
           k = N3
           DO j = S23B, N23B
!pgi$ unroll = n:8
              DO i = S13B, N13B
                 rel(i,j,k) = (bb(i,j,k) - cw33R(-1,k,g)*rel(i,j,k-1)) / cw33R(0,k,g)
              END DO
           END DO
        END IF
        !=====================================================================================================
        
     END IF
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     
  END DO
  
  
  END SUBROUTINE relaxation_Helmholtz
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE relaxation_Helmholtz_coarse(init_yes,n_relax,g,bb,rel) ! TEST!!! 2D-Variante fehlt noch ...
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  init_yes
  INTEGER, INTENT(in   ) ::  n_relax
  
  INTEGER, INTENT(in   ) ::  g
  
  REAL   , INTENT(in   ) ::  bb (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL   , INTENT(inout) ::  rel(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  INTEGER                ::  i, N1, N1R, N11R
  INTEGER                ::  j, N2, N2R, N22R
  INTEGER                ::  k, N3, N3R, N33R
  
  INTEGER                ::  r
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Dirichlet-Randbedingungen sind der Konvergenzrate wegen überall VOR der Relaxation        !
  !                eingefügt!                                                                                !
  !----------------------------------------------------------------------------------------------------------!
  
  N1 = NN(1,g)
  N2 = NN(2,g)
  N3 = NN(3,g)
  
  N1R  = N1 + d1R
  N2R  = N2 + d2R
  N3R  = N3 + d3R
  
  N11R = N1 + d11R
  N22R = N2 + d22R
  N33R = N3 + d33R
  
  
  ! Nur zum testen, falls Initialisierung unten fehlerhaft erscheint:
  !IF (init_yes) rel = 0.
  
  
  IF (init_yes) THEN
     IF (BC_1L .LE. 0) rel(S11R-1   ,S22R:N22R,S33R:N33R) = 0.
     IF (BC_2L .LE. 0) rel(S11R:N11R,S22R-1   ,S33R:N33R) = 0.
     IF (BC_3L .LE. 0) rel(S11R:N11R,S22R:N22R,S33R-1   ) = 0.
  END IF
  
  
  DO r = 1, n_relax
     
     IF (.NOT. (r == 1 .AND. init_yes)) CALL exchange_relax(g,0,0,0,0,.TRUE.,rel)
     
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     IF (direction == 1) THEN
        
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L == 1 .OR. BC_2L == 3) rel(S1R:N1R,1 ,S3R:N3R) = bb(S1R:N1R,1 ,S3R:N3R)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L == 2 .OR. BC_2L == 4) THEN
           j = 1
           IF (r == 1 .AND. init_yes) THEN
              DO k = S3R, N3R
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    rel(i,j,k) = bb(i,j,k) / cp22R(0,j,g)
                 END DO
              END DO
           ELSE
              DO k = S3R, N3R
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    rel(i,j,k) = (bb(i,j,k) - cp22R( 1,j,g)*rel(i,j+1,k)) / cp22R(0,j,g)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L == 1 .OR. BC_3L == 3) rel(S1R:N1R,S2R:N2R,1 ) = bb(S1R:N1R,S2R:N2R,1 )
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L == 2 .OR. BC_3L == 4) THEN
           k = 1
           IF (r == 1 .AND. init_yes) THEN
              DO j = S2R, N2R
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    rel(i,j,k) = bb(i,j,k) / cp33R(0,k,g)
                 END DO
              END DO
           ELSE
              DO j = S2R, N2R
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    rel(i,j,k) = (bb(i,j,k) - cp33R( 1,k,g)*rel(i,j,k+1)) / cp33R(0,k,g)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L .GT. 0) THEN
           i = 1
           IF (r == 1 .AND. init_yes) THEN
              DO k = S3R, N3R
!pgi$ unroll = n:8
                 DO j = S2R, N2R
                    rel(i,j,k) = bb(i,j,k) / cu11R(0,i,g)
                 END DO
              END DO
           ELSE
              DO k = S3R, N3R
!pgi$ unroll = n:8
                 DO j = S2R, N2R
                    rel(i,j,k) = (bb(i,j,k) - cu11R( 1,i,g)*rel(i+1,j,k)) / cu11R(0,i,g)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Feld --------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (r == 1 .AND. init_yes) THEN
           DO k = S33R, N33R
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                             &
                       &    cu11R(-1,i,g)*rel(i-1,j,k) + cp22R(-1,j,g)*rel(i,j-1,k) + cp33R(-1,k,g)*rel(i,j,k-1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cu11R(0,i,g) + cp22R(0,j,g) + cp33R(0,k,g)))
                 END DO
              END DO
           END DO
        ELSE
           DO k = S33R, N33R
             DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                &
                                &        cu11R(-1,i,g)*rel(i-1,j,k) + cu11R(1,i,g)*rel(i+1,j,k)     &
                                &      + cp22R(-1,j,g)*rel(i,j-1,k) + cp22R(1,j,g)*rel(i,j+1,k)     &
                                &      + cp33R(-1,k,g)*rel(i,j,k-1) + cp33R(1,k,g)*rel(i,j,k+1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cu11R(0,i,g) + cp22R(0,j,g) + cp33R(0,k,g)))
                 END DO
              END DO
           END DO
        END IF
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2U == 1 .OR. BC_2U == 3) rel(S1R:N1R,N2,S3R:N3R) = bb(S1R:N1R,N2,S3R:N3R)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2U == 2 .OR. BC_2U == 4) THEN
           j = N2
           DO k = S3R, N3R
!pgi$ unroll = n:8
              DO i = S1R, N1R
                 rel(i,j,k) = (bb(i,j,k) - cp22R(-1,j,g)*rel(i,j-1,k)) / cp22R(0,j,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U == 1 .OR. BC_3U == 3) rel(S1R:N1R,S2R:N2R,N3) = bb(S1R:N1R,S2R:N2R,N3)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U == 2 .OR. BC_3U == 4) THEN
           k = N3
           DO j = S2R, N2R
!pgi$ unroll = n:8
              DO i = S1R, N1R
                 rel(i,j,k) = (bb(i,j,k) - cp33R(-1,k,g)*rel(i,j,k-1)) / cp33R(0,k,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1U .GT. 0) THEN
           i = N1
           DO k = S3R, N3R
!pgi$ unroll = n:8
              DO j = S2R, N2R
                 rel(i,j,k) = (bb(i,j,k) - cu11R(-1,i,g)*rel(i-1,j,k)) / cu11R(0,i,g)
              END DO
           END DO
        END IF
        !=====================================================================================================
        
     END IF
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     IF (direction == 2) THEN
        
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L == 1 .OR. BC_1L == 3) rel(1 ,S2R:N2R,S3R:N3R) = bb(1 ,S2R:N2R,S3R:N3R)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L == 2 .OR. BC_1L == 4) THEN
           i = 1
           IF (r == 1 .AND. init_yes) THEN
              DO k = S3R, N3R
!pgi$ unroll = n:8
                 DO j = S2R, N2R
                    rel(i,j,k) = bb(i,j,k) / cp11R(0,i,g)
                 END DO
              END DO
           ELSE
              DO k = S3R, N3R
!pgi$ unroll = n:8
                 DO j = S2R, N2R
                    rel(i,j,k) = (bb(i,j,k) - cp11R( 1,i,g)*rel(i+1,j,k)) / cp11R(0,i,g)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L == 1 .OR. BC_3L == 3) rel(S1R:N1R,S2R:N2R,1 ) = bb(S1R:N1R,S2R:N2R,1 )
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L == 2 .OR. BC_3L == 4) THEN
           k = 1
           IF (r == 1 .AND. init_yes) THEN
              DO j = S2R, N2R
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    rel(i,j,k) = bb(i,j,k) / cp33R(0,k,g)
                 END DO
              END DO
           ELSE
              DO j = S2R, N2R
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    rel(i,j,k) = (bb(i,j,k) - cp33R( 1,k,g)*rel(i,j,k+1)) / cp33R(0,k,g)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L .GT. 0) THEN
           j = 1
           IF (r == 1 .AND. init_yes) THEN
              DO k = S3R, N3R
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    rel(i,j,k) = bb(i,j,k) / cv22R(0,j,g)
                 END DO
              END DO
           ELSE
              DO k = S3R, N3R
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    rel(i,j,k) = (bb(i,j,k) - cv22R( 1,j,g)*rel(i,j+1,k)) / cv22R(0,j,g)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Feld --------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (r == 1 .AND. init_yes) THEN
           DO k = S33R, N33R
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                             &
                       &    cp11R(-1,i,g)*rel(i-1,j,k) + cv22R(-1,j,g)*rel(i,j-1,k) + cp33R(-1,k,g)*rel(i,j,k-1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cp11R(0,i,g) + cv22R(0,j,g) + cp33R(0,k,g)))
                 END DO
              END DO
           END DO
        ELSE
           DO k = S33R, N33R
             DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                &
                                &        cp11R(-1,i,g)*rel(i-1,j,k) + cp11R(1,i,g)*rel(i+1,j,k)     &
                                &      + cv22R(-1,j,g)*rel(i,j-1,k) + cv22R(1,j,g)*rel(i,j+1,k)     &
                                &      + cp33R(-1,k,g)*rel(i,j,k-1) + cp33R(1,k,g)*rel(i,j,k+1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cp11R(0,i,g) + cv22R(0,j,g) + cp33R(0,k,g)))
                 END DO
              END DO
           END DO
        END IF
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1U == 1 .OR. BC_1U == 3) rel(N1,S2R:N2R,S3R:N3R) = bb(N1,S2R:N2R,S3R:N3R)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1U == 2 .OR. BC_1U == 4) THEN
           i = N1
           DO k = S3R, N3R
!pgi$ unroll = n:8
              DO j = S2R, N2R
                 rel(i,j,k) = (bb(i,j,k) - cp11R(-1,i,g)*rel(i-1,j,k)) / cp11R(0,i,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U == 1 .OR. BC_3U == 3) rel(S1R:N1R,S2R:N2R,N3) = bb(S1R:N1R,S2R:N2R,N3)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U == 2 .OR. BC_3U == 4) THEN
           k = N3
           DO j = S2R, N2R
!pgi$ unroll = n:8
              DO i = S1R, N1R
                 rel(i,j,k) = (bb(i,j,k) - cp33R(-1,k,g)*rel(i,j,k-1)) / cp33R(0,k,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2U .GT. 0) THEN
           j = N2
           DO k = S3R, N3R
!pgi$ unroll = n:8
              DO i = S1R, N1R
                 rel(i,j,k) = (bb(i,j,k) - cv22R(-1,j,g)*rel(i,j-1,k)) / cv22R(0,j,g)
              END DO
           END DO
        END IF
        !=====================================================================================================
        
     END IF
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     IF (direction == 3) THEN
        
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L == 1 .OR. BC_1L == 3) rel(1 ,S2R:N2R,S3R:N3R) = bb(1 ,S2R:N2R,S3R:N3R)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L == 2 .OR. BC_1L == 4) THEN
           i = 1
           IF (r == 1 .AND. init_yes) THEN
              DO k = S3R, N3R
!pgi$ unroll = n:8
                 DO j = S2R, N2R
                    rel(i,j,k) = bb(i,j,k) / cp11R(0,i,g)
                 END DO
              END DO
           ELSE
              DO k = S3R, N3R
!pgi$ unroll = n:8
                 DO j = S2R, N2R
                    rel(i,j,k) = (bb(i,j,k) - cp11R( 1,i,g)*rel(i+1,j,k)) / cp11R(0,i,g)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L == 1 .OR. BC_2L == 3) rel(S1R:N1R,1 ,S3R:N3R) = bb(S1R:N1R,1 ,S3R:N3R)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L == 2 .OR. BC_2L == 4) THEN
           j = 1
           IF (r == 1 .AND. init_yes) THEN
              DO k = S3R, N3R
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    rel(i,j,k) = bb(i,j,k) / cp22R(0,j,g)
                 END DO
              END DO
           ELSE
              DO k = S3R, N3R
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    rel(i,j,k) = (bb(i,j,k) - cp22R( 1,j,g)*rel(i,j+1,k)) / cp22R(0,j,g)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L .GT. 0) THEN
           k = 1
           IF (r == 1 .AND. init_yes) THEN
              DO j = S2R, N2R
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    rel(i,j,k) = bb(i,j,k) / cw33R(0,k,g)
                 END DO
              END DO
           ELSE
              DO j = S2R, N2R
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    rel(i,j,k) = (bb(i,j,k) - cw33R( 1,k,g)*rel(i,j,k+1)) / cw33R(0,k,g)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Feld --------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (r == 1 .AND. init_yes) THEN
           DO k = S33R, N33R
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                             &
                       &    cp11R(-1,i,g)*rel(i-1,j,k) + cp22R(-1,j,g)*rel(i,j-1,k) + cw33R(-1,k,g)*rel(i,j,k-1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cp11R(0,i,g) + cp22R(0,j,g) + cw33R(0,k,g)))
                 END DO
              END DO
           END DO
        ELSE
           DO k = S33R, N33R
             DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                &
                                &        cp11R(-1,i,g)*rel(i-1,j,k) + cp11R(1,i,g)*rel(i+1,j,k)     &
                                &      + cp22R(-1,j,g)*rel(i,j-1,k) + cp22R(1,j,g)*rel(i,j+1,k)     &
                                &      + cw33R(-1,k,g)*rel(i,j,k-1) + cw33R(1,k,g)*rel(i,j,k+1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cp11R(0,i,g) + cp22R(0,j,g) + cw33R(0,k,g)))
                 END DO
              END DO
           END DO
        END IF
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1U == 1 .OR. BC_1U == 3) rel(N1,S2R:N2R,S3R:N3R) = bb(N1,S2R:N2R,S3R:N3R)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1U == 2 .OR. BC_1U == 4) THEN
           i = N1
           DO k = S3R, N3R
!pgi$ unroll = n:8
              DO j = S2R, N2R
                 rel(i,j,k) = (bb(i,j,k) - cp11R(-1,i,g)*rel(i-1,j,k)) / cp11R(0,i,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2U == 1 .OR. BC_2U == 3) rel(S1R:N1R,N2,S3R:N3R) = bb(S1R:N1R,N2,S3R:N3R)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2U == 2 .OR. BC_2U == 4) THEN
           j = N2
           DO k = S3R, N3R
!pgi$ unroll = n:8
              DO i = S1R, N1R
                 rel(i,j,k) = (bb(i,j,k) - cp22R(-1,j,g)*rel(i,j-1,k)) / cp22R(0,j,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U .GT. 0) THEN
           k = N3
           DO j = S2R, N2R
!pgi$ unroll = n:8
              DO i = S1R, N1R
                 rel(i,j,k) = (bb(i,j,k) - cw33R(-1,k,g)*rel(i,j,k-1)) / cw33R(0,k,g)
              END DO
           END DO
        END IF
        !=====================================================================================================
        
     END IF
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     
  END DO
  
  
  END SUBROUTINE relaxation_Helmholtz_coarse
  
  
  
END MODULE mod_helmholtz
