!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!*************************************************************************************************************

!> module containing subroutines for computing the Laplacian operator. It uses modules mod_dims, mod-vars,
!! mod_diff and mod_exchange.
MODULE mod_laplace


  USE mod_dims
  USE mod_vars
  USE mod_diff
  USE mod_exchange
  
  
  PRIVATE
  
  PUBLIC product_div_grad, product_div_grad_transp
  PUBLIC product_div_grad_relax, relaxation_div_grad, relaxation_div_grad_inv
  PUBLIC handle_corner_Lap  
  
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> subroutine that computes the Laplacian operator 
  SUBROUTINE product_div_grad(phi,Lap)
  
  IMPLICIT NONE
  
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !> ?
  REAL   , INTENT(out  ) ::  Lap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !> Laplacian
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  INTEGER                ::  m
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Achtung: - Ein direktes Aufstellen des Laplace-Operators funktioniert i.A. nur sehr schlecht als Vor-    !
  !            konditionierer.                                                                               !
  !          - Das Feld "grad" kann nur eliminiert werden, wenn der div-grad-Stencil komplett ausmultipli-   !
  !            ziert wird (aktuell zu aufwendig).                                                            !
  !----------------------------------------------------------------------------------------------------------!
  
  
  DO m = 1, dimens
     
     CALL gradient        (m,phi,dig)
     CALL bc_extrapolation(m,dig    )
     CALL divergence      (m,dig,Lap)
     
  END DO
  
  IF (corner_yes) CALL handle_corner_Lap(1,Lap)
  
  
  END SUBROUTINE product_div_grad
  
  
  
  
  
  
  
  
  
  
  !> subroutine that computes the transpose of the Laplacian 
  SUBROUTINE product_div_grad_transp(phi,Lap)
  
  IMPLICIT NONE
  
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(out  ) ::  Lap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  INTEGER                ::  m
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Achtung: - Alle Operatoren werden transponiert und in umgekehrter Reihenfolge angewendet im Vergleich zu !
  !            "product_div_grad".                                                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (corner_yes) CALL handle_corner_Lap(1,phi)
  
  DO m = 1, dimens
     
     CALL divergence_transp      (m,phi,dig)
     CALL bc_extrapolation_transp(m,dig    )
     CALL gradient_transp        (m,dig,Lap)
     
  END DO
  
  
  END SUBROUTINE product_div_grad_transp
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE product_div_grad_relax(g,phi,Lap) ! TEST!!! aufraeumen und Variablen substituieren ...
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  g
  
  REAL   , INTENT(inout) ::  phi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL   , INTENT(out  ) ::  Lap(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  INTEGER                ::  i, N1, S1R, N1R, S11R, N11R
  INTEGER                ::  j, N2, S2R, N2R, S22R, N22R
  INTEGER                ::  k, N3, S3R, N3R, S33R, N33R
  
  INTEGER                ::  BC_1L, BC_1U, BC_2L, BC_2U, BC_3L, BC_3U
  
  !*****************************************************************************************
  N1 = NN(1,g)
  N2 = NN(2,g)
  N3 = NN(3,g)
  
  BC_1L = BC(1,1,g)
  BC_1U = BC(2,1,g)
  BC_2L = BC(1,2,g)
  BC_2U = BC(2,2,g)
  BC_3L = BC(1,3,g)
  BC_3U = BC(2,3,g)
  
  S1R  = SNB(1,1,g)
  S2R  = SNB(1,2,g)
  S3R  = SNB(1,3,g)
  
  N1R  = SNB(2,1,g)
  N2R  = SNB(2,2,g)
  N3R  = SNB(2,3,g)
  
  S11R = SNF(1,1,g)
  S22R = SNF(1,2,g)
  S33R = SNF(1,3,g)
  
  N11R = SNF(2,1,g)
  N22R = SNF(2,2,g)
  N33R = SNF(2,3,g)
  !*****************************************************************************************
  
  
  !-----------------------------------------------------------------------------------------------------!
  ! Anmerkung: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gewählt sind!        !
  !            - Direktes Verbinden mit Restriktion erscheint nicht sehr attraktiv, da Lap auf dem      !
  !              jeweils feineren Gitter ohnehin benötigt wird und coarse nur per MOD( , ) oder IF-     !
  !              Abfrage belegt werden könnte. Zum anderen Wird das redundante feinere Feld auch bei    !
  !              der Interpolation verwendet (ebenfalls zur Vereinfachung der Programmierung) und ist   !
  !              somit ohnehin vorhanden. Geschwindigkeitsmässig ist vermutlich auch nicht mehr viel zu !
  !              holen.                                                                                 !        
  !-----------------------------------------------------------------------------------------------------!
  
  
  CALL exchange_relax(g,0,0,0,0,.TRUE.,phi)
  
  
  !===========================================================================================================
  IF (dimens == 3) THEN
     
     DO k = S33R, N33R
        DO j = S22R, N22R
!pgi$ unroll = n:8
           DO i = S11R, N11R
              Lap(i,j,k) =  cdg1(-1,i,g)*phi(i-1,j,k) + cdg1(1,i,g)*phi(i+1,j,k)   &
                      &  +  cdg2(-1,j,g)*phi(i,j-1,k) + cdg2(1,j,g)*phi(i,j+1,k)   &
                      &  +  cdg3(-1,k,g)*phi(i,j,k-1) + cdg3(1,k,g)*phi(i,j,k+1)   &
                      &  + (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))*phi(i,j,k)
           END DO
        END DO
     END DO
     
  !===========================================================================================================
  ELSE
     
     DO k = S33R, N33R
        DO j = S22R, N22R
!pgi$ unroll = n:8
           DO i = S11R, N11R
              Lap(i,j,k) =  cdg1(-1,i,g)*phi(i-1,j,k) + cdg1(1,i,g)*phi(i+1,j,k)   &
                      &  +  cdg2(-1,j,g)*phi(i,j-1,k) + cdg2(1,j,g)*phi(i,j+1,k)   &
                      &  + (cdg1(0,i,g) + cdg2(0,j,g))*phi(i,j,k)
           END DO
        END DO
     END DO
     
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Randbedingungen =======================================================================================
  !===========================================================================================================
  
  IF (BC_1L .GT. 0) THEN
     i = 1
     DO k = S3R, N3R
!pgi$ unroll = n:8
        DO j = S2R, N2R
           Lap(i,j,k) = cdg1(0,i,g)*phi(i,j,k) + cdg1(1,i,g)*phi(i+1,j,k)
        END DO
     END DO
  END IF
  
  IF (BC_1U .GT. 0) THEN
     i = N1
     DO k = S3R, N3R
!pgi$ unroll = n:8
        DO j = S2R, N2R
           Lap(i,j,k) = cdg1(-1,i,g)*phi(i-1,j,k) + cdg1(0,i,g)*phi(i,j,k)
        END DO
     END DO
  END IF
  
  !===========================================================================================================
  
  IF (BC_2L .GT. 0) THEN
     j = 1
     DO k = S3R, N3R
!pgi$ unroll = n:8
        DO i = S1R, N1R
           Lap(i,j,k) = cdg2(0,j,g)*phi(i,j,k) + cdg2(1,j,g)*phi(i,j+1,k)
        END DO
     END DO
  END IF
  
  IF (BC_2U .GT. 0) THEN
     j = N2
     DO k = S3R, N3R
!pgi$ unroll = n:8
        DO i = S1R, N1R
           Lap(i,j,k) = cdg2(-1,j,g)*phi(i,j-1,k) + cdg2(0,j,g)*phi(i,j,k)
        END DO
     END DO
  END IF
  
  !===========================================================================================================
  
  IF (BC_3L .GT. 0) THEN
     k = 1
     DO j = S2R, N2R
!pgi$ unroll = n:8
        DO i = S1R, N1R
           Lap(i,j,k) = cdg3(0,k,g)*phi(i,j,k) + cdg3(1,k,g)*phi(i,j,k+1)
        END DO
     END DO
  END IF
  
  IF (BC_3U .GT. 0) THEN
     k = N3
     DO j = S2R, N2R
!pgi$ unroll = n:8
        DO i = S1R, N1R
           Lap(i,j,k) = cdg3(-1,k,g)*phi(i,j,k-1) + cdg3(0,k,g)*phi(i,j,k)
        END DO
     END DO
  END IF
  
  !===========================================================================================================
  
  
  IF (corner_yes) CALL handle_corner_Lap(g,Lap)
  
  
  END SUBROUTINE product_div_grad_relax
  
  
  
  
  
  
  
  
  
  ! Habe RB-Linienrelaxation herausgeschmissen, weil die Gewinne in der Geschwindigkeit zu gering waren.
  ! Gleichzeitig war der Code viel zu lang und unuebersichtlich.
  SUBROUTINE relaxation_div_grad(init_yes,n_relax,g,bb,rel) ! TEST!!! reine 2D-Variante fehlt noch ... ! TEST!!! aufraeumen und Variablen substituieren ...
  
  IMPLICIT NONE
  
  !*************************************************************************************************
  INTEGER                ::  BC_1L, BC_1U, BC_2L, BC_2U, BC_3L, BC_3U
  !*************************************************************************************************
  
  LOGICAL, INTENT(in   ) ::  init_yes
  INTEGER, INTENT(in   ) ::  n_relax
  
  INTEGER, INTENT(in   ) ::  g
  
  REAL   , INTENT(in   ) ::  bb  (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL   , INTENT(inout) ::  rel (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL                   ::  comp(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U)) ! FELD!!!
  
  INTEGER                ::  i, N1, S1R, N1R, S11R, N11R
  INTEGER                ::  j, N2, S2R, N2R, S22R, N22R
  INTEGER                ::  k, N3, S3R, N3R, S33R, N33R
  
  INTEGER                ::  r, ss
  REAL                   ::  mult
  
  ! TEST!!! herausziehen?
  ! RBGS_mode =< 0 : naturally ordered Gauss-Seidel (slow on CPU, cannot be vectorized)
  ! RBGS_mode == 1 : 3D-Red-Black Gauss-Seidel (faster on CPU due to vectorization, parallelization-independent)
  ! RBGS_mode >= 2 : 2D-Red-Black Gauss-Seidel (normally the fastest variant on CPU due to better cache utilization)
  INTEGER, PARAMETER     ::  RBGS_mode = 2
  LOGICAL, PARAMETER     ::  SOR_yes = .FALSE.
  !REAL   , PARAMETER     ::  omega   = 0.8! 1.2 ! 1.27 !1.27  ! omega = 0.9
  REAL                   ::  omega
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Optimization annotations:                                                                                !
  !    - line relaxation runs on Rosa (Cray XT5, Istanbul hex-core processors) in single core mode with      !
  !      about 5.6-5.8% peak performance in direction 1, 4.2-4.6% in direction 2, 3.9-4.0% in direction 3.   !
  !    - line relaxation in direction 3 is about 50% more expensive than in direction 1; direction 2 is      !
  !      somewhere in between.                                                                               !
  !    - line relaxation in direction 1 is typically 30% of the total execution time(!).                     !
  !    - RB-GS converges slightly worse than standard GS.                                                    !
  !    - RB-GS allows full vectorization (execpt some operations in direction 1), provided that the loops    !
  !      are reordered to let direction 1 be the innermost (speedup was not as good as hoped (<5%) for       !
  !      direction 1 and 2, direction 3 was not tested).                                                     !
  !    - unrolling the non-vectorizable loops 8 times gives roughly the best results.                        !
  !    - IF statements are automatically pulled out of the loops by the compiler.                            !
  !----------------------------------------------------------------------------------------------------------!
  
  ! Allgemein, pgf90: "Invariant if transformation" <==> IF-Statements werden aus dem Loop herauskopiert, muss man nicht von Hand machen!
  
  
  !*****************************************************************************************
  N1 = NN(1,g)
  N2 = NN(2,g)
  N3 = NN(3,g)
  
  BC_1L = BC(1,1,g)
  BC_1U = BC(2,1,g)
  BC_2L = BC(1,2,g)
  BC_2U = BC(2,2,g)
  BC_3L = BC(1,3,g)
  BC_3U = BC(2,3,g)
  
  S1R  = SNB(1,1,g)
  S2R  = SNB(1,2,g)
  S3R  = SNB(1,3,g)
  
  N1R  = SNB(2,1,g)
  N2R  = SNB(2,2,g)
  N3R  = SNB(2,3,g)
  
  S11R = SNF(1,1,g)
  S22R = SNF(1,2,g)
  S33R = SNF(1,3,g)
  
  N11R = SNF(2,1,g)
  N22R = SNF(2,2,g)
  N33R = SNF(2,3,g)
  !*****************************************************************************************
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gewählt sind!           !
  !              - Bei der Initialisierung müssen die Intervallgrenzen SiR:NiR anstelle von SiiR:NiiR        !
  !                gewählt werden, da beim  Aufbau der RHS ("vec") innerhalb der Linienrelaxation IMMER auch !
  !                die Randbereiche SiR und NiR aufgebaut aber ggf. mit einer Korrektur wieder überschrieben !
  !                werden. Ansonsten wäre eine Initialisierung allein im Feldbereich SiiR:NiiR ausreichend,  !
  !                kann aber z.B. zu Floating Point Exceptions führen (die für die eigentliche Rechnung      !
  !                allerdings irrelevant wären)!                                                             !
  !              - obere Stirnflächen werden bei der Initialisierung ebenfalls berücksichtigt, da ggf.       !
  !                verschiedene Richtungen nacheinander bearbeitet werden und dies beim Aufbauen der rechten !
  !                Seite sonst berücksichtigt werden müsste. ==> geringerer Mehraufwand beim Rechnen,        !
  !                weniger Programmierdetails.                                                               !
  !              - LU-Zerlegung (d.h. band, mult) kann gespeichert werden, wenn mindestens eine Richtung     !
  !                äquidistant ist. Der Lösungsaufwand würde sich etwa halbieren!!                           !
  !              - "r == 1 .AND. init_yes" sollte idealerweise aus den Schleifen herausgezogen werden, was   !
  !                hier aber aus Gründen der Übersicht bisher nicht ausgeführt wurde.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  omega = 1.
  
  
  IF (Jacobi_yes) THEN
  
  DO r = 1, n_relax
     
     IF (r == 1 .AND. init_yes) THEN
        !=====================================================================================================
        IF (BC_1L .GT. 0) THEN
           i = 1
           DO k = S33R, N33R
!pgi$ unroll = n:8
              DO j = S22R, N22R
                 rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L .GT. 0) THEN
           j = 1
           DO k = S33R, N33R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L .GT. 0 .AND. dimens == 3) THEN
           k = 1
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
              END DO
           END DO
        END IF
        !=====================================================================================================
        DO k = S33R, N33R
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 rel(i,j,k) = omega*bb(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
              END DO
           END DO
        END DO
        !=====================================================================================================
        IF (BC_1U .GT. 0) THEN
           i = N1
           DO k = S33R, N33R
!pgi$ unroll = n:8
              DO j = S22R, N22R
                 rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2U .GT. 0) THEN
           j = N2
           DO k = S33R, N33R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U .GT. 0 .AND. dimens == 3) THEN
           k = N3
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
              END DO
           END DO
        END IF
        !=====================================================================================================
        
     ELSE
        
        CALL exchange_relax(g,0,0,0,0,.TRUE.,rel)
        
        !=====================================================================================================
        IF (BC_1L .GT. 0) THEN
           i = 1
           DO k = S33R, N33R
!pgi$ unroll = n:8
              DO j = S22R, N22R
                 comp(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L .GT. 0) THEN
           j = 1
           DO k = S33R, N33R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 comp(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L .GT. 0 .AND. dimens == 3) THEN
           k = 1
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 comp(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
              END DO
           END DO
        END IF
        !=====================================================================================================
        DO k = S33R, N33R
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 comp(i,j,k) = omega*(bb(i,j,k)                                                &
                             &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                             &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                             &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                             &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
              END DO
           END DO
        END DO
        !=====================================================================================================
        IF (BC_1U .GT. 0) THEN
           i = N1
           DO k = S33R, N33R
!pgi$ unroll = n:8
              DO j = S22R, N22R
                 comp(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2U .GT. 0) THEN
           j = N2
           DO k = S33R, N33R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 comp(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U .GT. 0 .AND. dimens == 3) THEN
           k = N3
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 comp(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
              END DO
           END DO
        END IF
        !=====================================================================================================
        
        rel(1:N1,1:N2,1:N3) = comp(1:N1,1:N2,1:N3)
        
     END IF
     
  END DO
  
  
  !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************
  ELSE
  
  
  IF (init_yes) THEN
     IF (BC_1L .LE. 0) rel(S11R-1 ,S2R:N2R,S3R:N3R) = 0.
     IF (BC_1U .LE. 0) rel(N11R+1 ,S2R:N2R,S3R:N3R) = 0.
     IF (BC_2L .LE. 0) rel(S1R:N1R,S22R-1 ,S3R:N3R) = 0.
     IF (BC_2U .LE. 0) rel(S1R:N1R,N22R+1 ,S3R:N3R) = 0.
     IF (BC_3L .LE. 0) rel(S1R:N1R,S2R:N2R,S33R-1 ) = 0.
     IF (BC_3U .LE. 0) rel(S1R:N1R,S2R:N2R,N33R+1 ) = 0.
     
     !rel = 0. ! TEST!!!  Nur zum testen, falls Initialisierung unten fehlerhaft erscheint:
  END IF
  
  
  DO r = 1, n_relax
     
     IF (.NOT. (r == 1 .AND. init_yes)) CALL exchange_relax(g,0,0,0,0,.TRUE.,rel)
     
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     IF (impl_dir(1) == 0 .AND. impl_dir(2) == 0 .AND. impl_dir(3) == 0) THEN
        
        IF (RBGS_mode .GT. 0) THEN
           !-------------------------------------------------------------------------------------------------!
           ! Hinweis:                                                                                        !
           ! Die Red-Black-Sortierung sollte für eine gute Performance global sein, d.h. ueber die Blöcke    !
           ! hinweg implementiert sein. Die aktuelle Umsetzung macht sich zunutze, dass im Parallelbetrieb   !
           ! die Anzahl der Gitterpunkte in jeder Raumrichtung automatisch eine gerade Zahl ist (die Tat-    !
           ! sache, dass in Randblöcken ggf. eine ungerade Anzahl Gitterpunkte abgearbeitet wird hat damit   !
           ! nichts zu tun!). Daher kann sich der Shift "ss" in jedem Block einfach auf den lokalen Gitter-  !
           ! punkt rel(1,1,1) beziehen.                                                                      !
           !                                                                                                 !
           ! Die aktuelle Schreibweise in den Schleifen ist äquivalent zu:                                   !
           ! DO k = S33R, N33R                                                                               !
           !    DO j = S22R, N22R                                                                            !
           !       DO i = S11R, N11R                                                                         !
           !          IF (MOD(i+j+k) == 1) ... ! (red)                                                       !
           !          ! bzw.                                                                                 !
           !          IF (MOD(i+j+k) == 0) ... ! (black)                                                     !
           !       END DO                                                                                    !
           !    END DO                                                                                       !
           ! END DO                                                                                          !
           !-------------------------------------------------------------------------------------------------!
           
           !==================================================================================================
           !==================================================================================================
           IF (BC_1L .GT. 0) THEN
              i = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
                    ss = MOD(i+k+S22R+1,2)
!pgi$ unroll = n:8
                    DO j = ss+S22R, N22R, 2
                       rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
                    ss = MOD(i+k+S22R+1,2)
!pgi$ unroll = n:8
                    DO j = ss+S22R, N22R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_1U .GT. 0) THEN
              i = N1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
                    ss = MOD(i+k+S22R+1,2)
!pgi$ unroll = n:8
                    DO j = ss+S22R, N22R, 2
                       rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
                    ss = MOD(i+k+S22R+1,2)
!pgi$ unroll = n:8
                    DO j = ss+S22R, N22R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !==================================================================================================
           IF (BC_2L .GT. 0) THEN
              j = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_2U .GT. 0) THEN
              j = N2
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !==================================================================================================
           IF (RBGS_mode == 1) THEN
           !--------------------------------------------------------------------------------------------------
           IF (BC_3L .GT. 0) THEN
              k = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                    END DO
                 END DO
              ELSE
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_3U .GT. 0) THEN
              k = N3
              IF (r == 1 .AND. init_yes) THEN
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                    END DO
                 END DO
              ELSE
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           ELSE
           !--------------------------------------------------------------------------------------------------
           IF (BC_3L .GT. 0) THEN
              k = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                    END DO
                 END DO
              ELSE
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           END IF
           !==================================================================================================
           IF (RBGS_mode == 1) THEN
           !--------------------------------------------------------------------------------------------------
           IF (r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R+1,2) ! "+k" ist offenbar elementar für grosse Konvergenzrate. Symmetrie zu relaxation_div_grad_inv ist hier interessanterweise kontraproduktiv!
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*bb(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       ELSE
                          rel(i,j,k) =       bb(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
              END DO
              !CALL exchange_relax(g,0,0,0,0,.TRUE.,rel) ! TEST!!! Austausch bringt praktisch nichts, macht aber die Glättung unabhängig von der Parallelisierung (sonst teilweise Jacobi-Iteration)!
              DO k = S33R, N33R
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                       ELSE
                          rel(i,j,k) =       (bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                       ELSE
                          rel(i,j,k) =       (bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
              END DO
              !CALL exchange_relax(g,0,0,0,0,.TRUE.,rel) ! TEST!!! Austausch bringt praktisch nichts, macht aber die Glättung unabhängig von der Parallelisierung (sonst teilweise Jacobi-Iteration)!
              DO k = S33R, N33R
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                       ELSE
                          rel(i,j,k) =       (bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
              END DO
           END IF
           !--------------------------------------------------------------------------------------------------
           ELSE
           !--------------------------------------------------------------------------------------------------
           IF (r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R+1,2) ! "+k" ist offenbar elementar für grosse Konvergenzrate. Symmetrie zu relaxation_div_grad_inv ist hier interessanterweise kontraproduktiv!
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       ELSE
                          rel(i,j,k) =       (bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1)                           )    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                       ELSE
                          rel(i,j,k) =       (bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1)                           )    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                       ELSE
                          rel(i,j,k) =       (bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                       ELSE
                          rel(i,j,k) =       (bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
              END DO
           END IF
           !--------------------------------------------------------------------------------------------------
           END IF
           !==================================================================================================
           IF (BC_1L .GT. 0) THEN
              i = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
                    ss = MOD(i+k+S22R,2)
!pgi$ unroll = n:8
                    DO j = ss+S22R, N22R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
                    ss = MOD(i+k+S22R,2)
!pgi$ unroll = n:8
                    DO j = ss+S22R, N22R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_1U .GT. 0) THEN
              i = N1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
                    ss = MOD(i+k+S22R,2)
!pgi$ unroll = n:8
                    DO j = ss+S22R, N22R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
                    ss = MOD(i+k+S22R,2)
!pgi$ unroll = n:8
                    DO j = ss+S22R, N22R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !==================================================================================================
           IF (BC_2L .GT. 0) THEN
              j = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_2U .GT. 0) THEN
              j = N2
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !==================================================================================================
           IF (RBGS_mode == 1) THEN
           !--------------------------------------------------------------------------------------------------
           IF (BC_3L .GT. 0) THEN
              k = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g)
                    END DO
                 END DO
              ELSE
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_3U .GT. 0) THEN
              k = N3
              IF (r == 1 .AND. init_yes) THEN
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g)
                    END DO
                 END DO
              ELSE
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           ELSE
           !--------------------------------------------------------------------------------------------------
           IF (BC_3U .GT. 0) THEN
              k = N3
              IF (r == 1 .AND. init_yes) THEN
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g)
                    END DO
                 END DO
              ELSE
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           END IF
           !==================================================================================================
           !==================================================================================================
           
        ELSE
           
           !==================================================================================================
           !==================================================================================================
           IF (BC_1L .GT. 0) THEN
              i = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
!pgi$ unroll = n:8
                    DO j = S22R, N22R
                       rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
!pgi$ unroll = n:8
                    DO j = S22R, N22R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_2L .GT. 0) THEN
              j = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_3L .GT. 0) THEN
              k = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                    END DO
                 END DO
              ELSE
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !==================================================================================================
           IF (r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       IF (SOR_yes) THEN
                          rel(i,j,k) = bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k) - cdg2(-1,j,g)*rel(i,j-1,k) - cdg3(-1,k,g)*rel(i,j,k-1)
                          rel(i,j,k) = omega*rel(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       ELSE
                          rel(i,j,k) = bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k) - cdg2(-1,j,g)*rel(i,j-1,k) - cdg3(-1,k,g)*rel(i,j,k-1)
                          rel(i,j,k) = rel(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                       ELSE
                          rel(i,j,k) =       (bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
              END DO
           END IF
           !==================================================================================================
           IF (BC_1U .GT. 0) THEN
              i = N1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
!pgi$ unroll = n:8
                    DO j = S22R, N22R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
!pgi$ unroll = n:8
                    DO j = S22R, N22R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_2U .GT. 0) THEN
              j = N2
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_3U .GT. 0) THEN
              k = N3
              IF (r == 1 .AND. init_yes) THEN
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g)
                    END DO
                 END DO
              ELSE
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !==================================================================================================
           !==================================================================================================
        END IF
        
     END IF
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     IF (impl_dir(1) == 1) THEN
        
        !=====================================================================================================
        IF (BC_2L .GT. 0) THEN
           j = 1
           IF (r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L .GT. 0) THEN
           k = 1
           IF (r == 1 .AND. init_yes) THEN
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                 END DO
              END DO
           ELSE
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
        DO k = S33R, N33R
           DO j = S22R, N22R
              
              !--- Diagonalelement / rechte Seite aufbauen ---------------------------------------------------
              IF (r == 1 .AND. init_yes) THEN
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    band1(1,i) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                    band1(2,i) = bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k) - cdg3(-1,k,g)*rel(i,j,k-1)
                 END DO
              ELSE
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    band1(1,i) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                    band1(2,i) = bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k) - cdg3(-1,k,g)*rel(i,j,k-1)   &
                                    &      - cdg2( 1,j,g)*rel(i,j+1,k) - cdg3( 1,k,g)*rel(i,j,k+1)
                 END DO
              END IF
              
              !--- SOR (1) -----------------------------------------------------------------------------------
              IF (SOR_yes .AND. .NOT. (r == 1 .AND. init_yes)) THEN
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    SOR1(i) = rel(i,j,k)
                 END DO
              END IF
              
              
              !--- mit RB überschreiben ----------------------------------------------------------------------
              IF (BC_1L .GT. 0) THEN
                 band1(1,S1R) = cdg1(0,S1R,g)
                 band1(2,S1R) = bb  (S1R,j,k)
              ELSE
                 band1(2,S1R) = band1(2,S1R) - cdg1(-1,S1R,g)*rel(S1R-1,j,k)
              END IF
              
              IF (BC_1U .GT. 0) THEN
                 band1(1,N1R) = cdg1(0,N1R,g)
                 band1(2,N1R) = bb  (N1R,j,k)
              ELSE
                 band1(2,N1R) = band1(2,N1R) - cdg1( 1,N1R,g)*rel(N1R+1,j,k)
              END IF
              
              
              !--- Gauss-Elimination hoch --------------------------------------------------------------------
!pgi$ unroll = n:8
              DO i = S1R+1, N1R
                 mult       = cdg1(-1,i,g) / band1(1,i-1)
                 band1(1,i) = band1(1,i) - mult * cdg1(1,i-1,g)
                 band1(2,i) = band1(2,i) - mult * band1(2,i-1)
              END DO
              
              
              !--- Gauss-Elimination runter ------------------------------------------------------------------
              rel(N1R,j,k) = band1(2,N1R) / band1(1,N1R)
!pgi$ unroll = n:8
              DO i = N1R-1, S1R, -1
                 rel(i,j,k) = (band1(2,i) - cdg1(1,i,g)*rel(i+1,j,k)) / band1(1,i)
              END DO
              
              
              !--- SOR (2) -----------------------------------------------------------------------------------
              IF (SOR_yes) THEN ! TEST!!! Das koennte man auch besser mit der vorherigen Schleife verbinden ...
                 IF (r == 1 .AND. init_yes) THEN
!pgi$ unroll = n:8
                    DO i = S1R, N1R
                       rel(i,j,k) = omega*rel(i,j,k)
                    END DO
                 ELSE
!pgi$ unroll = n:8
                    DO i = S1R, N1R
                       rel(i,j,k) = omega*rel(i,j,k) + (1.-omega)*SOR1(i)
                    END DO
                 END IF
              END IF
              
           END DO
        END DO
        !=====================================================================================================
        IF (BC_2U .GT. 0) THEN
           j = N2
           IF (r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g)
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U .GT. 0) THEN
           k = N3
           IF (r == 1 .AND. init_yes) THEN
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g)
                 END DO
              END DO
           ELSE
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
     END IF
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     IF (impl_dir(2) == 1) THEN
        
        IF (.NOT. (impl_dir(1) == 0 .AND. r == 1 .AND. init_yes)) CALL exchange_relax(g,0,0,0,0,.TRUE.,rel) ! TEST!!! Zwischenaktualisierung auch bei Helmholtz-Routinen einbauen!
        !=====================================================================================================
        IF (BC_1L .GT. 0) THEN
           i = 1
           IF (impl_dir(1) == 0 .AND. r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO j = S22R, N22R
                    rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO j = S22R, N22R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L .GT. 0) THEN
           k = 1
           IF (impl_dir(1) == 0 .AND. r == 1 .AND. init_yes) THEN
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                 END DO
              END DO
           ELSE
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
        DO k = S33R, N33R
           DO i = S11R, N11R
              
              !--- Diagonalelement / rechte Seite aufbauen ---------------------------------------------------
              IF (impl_dir(1) == 0 .AND. r == 1 .AND. init_yes) THEN
!pgi$ unroll = n:8
                 DO j = S2R, N2R
                    band2(1,j) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                    band2(2,j) = bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)
                 END DO
              ELSE
!pgi$ unroll = n:8
                 DO j = S2R, N2R
                    band2(1,j) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                    band2(2,j) = bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)   &
                                    &      - cdg1( 1,i,g)*rel(i+1,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)
                 END DO
              END IF
              
              !--- SOR (1) -----------------------------------------------------------------------------------
              IF (SOR_yes .AND. .NOT. (impl_dir(1) == 0 .AND. r == 1 .AND. init_yes)) THEN
!pgi$ unroll = n:8
                 DO j = S2R, N2R
                    SOR2(j) = rel(i,j,k)
                 END DO
              END IF
              
              
              !--- mit RB überschreiben ----------------------------------------------------------------------
              IF (BC_2L .GT. 0) THEN
                 band2(1,S2R) = cdg2(0,S2R,g)
                 band2(2,S2R) = bb  (i,S2R,k)
              ELSE
                 band2(2,S2R) = band2(2,S2R) - cdg2(-1,S2R,g)*rel(i,S2R-1,k)
              END IF
              
              IF (BC_2U .GT. 0) THEN
                 band2(1,N2R) = cdg2(0,N2R,g)
                 band2(2,N2R) = bb  (i,N2R,k)
              ELSE
                 band2(2,N2R) = band2(2,N2R) - cdg2( 1,N2R,g)*rel(i,N2R+1,k)
              END IF
              
              
              !--- Gauss-Elimination hoch --------------------------------------------------------------------
!pgi$ unroll = n:8
              DO j = S2R+1, N2R
                 mult       = cdg2(-1,j,g) / band2(1,j-1)
                 band2(1,j) = band2(1,j) - mult * cdg2(1,j-1,g)
                 band2(2,j) = band2(2,j) - mult * band2(2,j-1)
              END DO
              
              
              !--- Gauss-Elimination runter ------------------------------------------------------------------
              rel(i,N2R,k) = band2(2,N2R) / band2(1,N2R)
!pgi$ unroll = n:8
              DO j = N2R-1, S2R, -1
                 rel(i,j,k) = (band2(2,j) - cdg2(1,j,g)*rel(i,j+1,k)) / band2(1,j)
              END DO
              
              
              !--- SOR (2) -----------------------------------------------------------------------------------
              IF (SOR_yes) THEN
                 IF (impl_dir(1) == 0 .AND. r == 1 .AND. init_yes) THEN
!pgi$ unroll = n:8
                    DO j = S2R, N2R
                       rel(i,j,k) = omega*rel(i,j,k)
                    END DO
                 ELSE
!pgi$ unroll = n:8
                    DO j = S2R, N2R
                       rel(i,j,k) = omega*rel(i,j,k) + (1.-omega)*SOR2(j)
                    END DO
                 END IF
              END IF
              
           END DO
        END DO
        !=====================================================================================================
        IF (BC_1U .GT. 0) THEN
           i = N1
           IF (impl_dir(1) == 0 .AND. r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO j = S22R, N22R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g)
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO j = S22R, N22R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U .GT. 0) THEN
           k = N3
           IF (impl_dir(1) == 0 .AND. r == 1 .AND. init_yes) THEN
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g)
                 END DO
              END DO
           ELSE
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
        
     END IF
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     IF (impl_dir(3) == 1) THEN
        
        IF (.NOT. (impl_dir(1) == 0 .AND. impl_dir(2) == 0 .AND. r == 1 .AND. init_yes)) CALL exchange_relax(g,0,0,0,0,.TRUE.,rel)
        !=====================================================================================================
        IF (BC_1L .GT. 0) THEN
           i = 1
           IF (impl_dir(1) == 0 .AND. impl_dir(2) == 0 .AND. r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO j = S22R, N22R
                    rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO j = S22R, N22R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L .GT. 0) THEN
           j = 1
           IF (impl_dir(1) == 0 .AND. impl_dir(2) == 0 .AND. r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
        DO j = S22R, N22R
           DO i = S11R, N11R
              
              !--- Diagonalelement / rechte Seite aufbauen ---------------------------------------------------
              IF (impl_dir(1) == 0 .AND. impl_dir(2) == 0 .AND. r == 1 .AND. init_yes) THEN
!pgi$ unroll = n:8
                 DO k = S3R, N3R
                    band3(1,k) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                    band3(2,k) = bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)
                 END DO
              ELSE
!pgi$ unroll = n:8
                 DO k = S3R, N3R
                    band3(1,k) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                    band3(2,k) = bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)   &
                                    &      - cdg1( 1,i,g)*rel(i+1,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)
                 END DO
              END IF
              
              !--- SOR (1) -----------------------------------------------------------------------------------
              IF (SOR_yes .AND. .NOT. (impl_dir(1) == 0 .AND. impl_dir(2) == 0 .AND. r == 1 .AND. init_yes)) THEN
!pgi$ unroll = n:8
                 DO k = S3R, N3R
                    SOR3(k) = rel(i,j,k)
                 END DO
              END IF
              
              
              !--- mit RB überschreiben ----------------------------------------------------------------------
              IF (BC_3L .GT. 0) THEN
                 band3(1,S3R) = cdg3(0,S3R,g)
                 band3(2,S3R) = bb  (i,j,S3R)
              ELSE
                 band3(2,S3R) = band3(2,S3R) - cdg3(-1,S3R,g)*rel(i,j,S3R-1)
              END IF
              
              IF (BC_3U .GT. 0) THEN
                 band3(1,N3R) = cdg3(0,N3R,g)
                 band3(2,N3R) = bb  (i,j,N3R)
              ELSE
                 band3(2,N3R) = band3(2,N3R) - cdg3( 1,N3R,g)*rel(i,j,N3R+1)
              END IF
              
              
              !--- Gauss-Elimination hoch --------------------------------------------------------------------
!pgi$ unroll = n:8
              DO k = S3R+1, N3R
                 mult       = cdg3(-1,k,g) / band3(1,k-1)
                 band3(1,k) = band3(1,k) - mult * cdg3(1,k-1,g)
                 band3(2,k) = band3(2,k) - mult * band3(2,k-1)
              END DO
              
              
              !--- Gauss-Elimination runter ------------------------------------------------------------------
              rel(i,j,N3R) = band3(2,N3R) / band3(1,N3R)
!pgi$ unroll = n:8
              DO k = N3R-1, S3R, -1
                 rel(i,j,k) = (band3(2,k) - cdg3(1,k,g)*rel(i,j,k+1)) / band3(1,k)
              END DO
              
              
              !--- SOR (2) -----------------------------------------------------------------------------------
              IF (SOR_yes) THEN
                 IF (impl_dir(1) == 0 .AND. impl_dir(2) == 0 .AND. r == 1 .AND. init_yes) THEN
!pgi$ unroll = n:8
                    DO k = S3R, N3R
                       rel(i,j,k) = omega*rel(i,j,k)
                    END DO
                 ELSE
!pgi$ unroll = n:8
                    DO k = S3R, N3R
                       rel(i,j,k) = omega*rel(i,j,k) + (1.-omega)*SOR3(k)
                    END DO
                 END IF
              END IF
              
           END DO
        END DO
        !=====================================================================================================
        IF (BC_1U .GT. 0) THEN
           i = N1
           IF (impl_dir(1) == 0 .AND. impl_dir(2) == 0 .AND. r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO j = S22R, N22R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g)
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO j = S22R, N22R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2U .GT. 0) THEN
           j = N2
           IF (impl_dir(1) == 0 .AND. impl_dir(2) == 0 .AND. r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g)
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
        
     END IF
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     
  END DO
  !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************
  
  END IF
  
  
  IF (corner_yes) CALL handle_corner_Lap(g,rel)
  
  
  END SUBROUTINE relaxation_div_grad
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE relaxation_div_grad_inv(init_yes,n_relax,g,bb,rel) ! TEST!!! 2D-Variante fehlt noch ... ! TEST!!! aufraeumen und Variablen substituieren ...
  
  IMPLICIT NONE
  
  !*************************************************************************************************
  INTEGER                ::  BC_1L, BC_1U, BC_2L, BC_2U, BC_3L, BC_3U
  !*************************************************************************************************
  
  LOGICAL, INTENT(in   ) ::  init_yes
  INTEGER, INTENT(in   ) ::  n_relax
  
  INTEGER, INTENT(in   ) ::  g
  
  REAL   , INTENT(in   ) ::  bb  (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL   , INTENT(inout) ::  rel (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL                   ::  comp(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U)) ! FELD!!!
  
  INTEGER                ::  i, N1, S1R, N1R, S11R, N11R
  INTEGER                ::  j, N2, S2R, N2R, S22R, N22R
  INTEGER                ::  k, N3, S3R, N3R, S33R, N33R
  
  INTEGER                ::  r, ss
  REAL                   ::  mult
  
  ! TEST!!! herausziehen?
  ! RBGS_mode =< 0 : naturally ordered Gauss-Seidel (slow on CPU, cannot be vectorized)
  ! RBGS_mode == 1 : 3D-Red-Black Gauss-Seidel (faster on CPU due to vectorization, parallelization-independent)
  ! RBGS_mode >= 2 : 2D-Red-Black Gauss-Seidel (normally the fastest variant on CPU due to better cache utilization)
  INTEGER, PARAMETER     ::  RBGS_mode = 2
  LOGICAL, PARAMETER     ::  SOR_yes = .FALSE.
  !REAL   , PARAMETER     ::  omega   = 0.8! 1.2 ! 1.27 !1.27  ! omega = 0.9
  REAL                   ::  omega
  
  ! Allgemein, pgf90: "Invariant if transformation" <==> IF-Statements werden aus dem Loop herauskopiert, muss man nicht von Hand machen!
  
  
  !*****************************************************************************************
  N1 = NN(1,g)
  N2 = NN(2,g)
  N3 = NN(3,g)
  
  BC_1L = BC(1,1,g)
  BC_1U = BC(2,1,g)
  BC_2L = BC(1,2,g)
  BC_2U = BC(2,2,g)
  BC_3L = BC(1,3,g)
  BC_3U = BC(2,3,g)
  
  S1R  = SNB(1,1,g)
  S2R  = SNB(1,2,g)
  S3R  = SNB(1,3,g)
  
  N1R  = SNB(2,1,g)
  N2R  = SNB(2,2,g)
  N3R  = SNB(2,3,g)
  
  S11R = SNF(1,1,g)
  S22R = SNF(1,2,g)
  S33R = SNF(1,3,g)
  
  N11R = SNF(2,1,g)
  N22R = SNF(2,2,g)
  N33R = SNF(2,3,g)
  !*****************************************************************************************
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gewählt sind!           !
  !              - Bei der Initialisierung müssen die Intervallgrenzen SiR:NiR anstelle von SiiR:NiiR        !
  !                gewählt werden, da beim  Aufbau der RHS ("vec") innerhalb der Linienrelaxation IMMER auch !
  !                die Randbereiche SiR und NiR aufgebaut aber ggf. mit einer Korrektur wieder überschrieben !
  !                werden. Ansonsten wäre eine Initialisierung allein im Feldbereich SiiR:NiiR ausreichend,  !
  !                kann aber z.B. zu Floating Point Exceptions führen (die für die eigentliche Rechnung      !
  !                allerdings irrelevant wären)!                                                             !
  !              - obere Stirnflächen werden bei der Initialisierung ebenfalls berücksichtigt, da ggf.       !
  !                verschiedene Richtungen nacheinander bearbeitet werden und dies beim Aufbauen der rechten !
  !                Seite sonst berücksichtigt werden müsste. ==> geringerer Mehraufwand beim Rechnen,        !
  !                weniger Programmierdetails.                                                               !
  !              - LU-Zerlegung (d.h. band, mult) kann gespeichert werden, wenn mindestens eine Richtung     !
  !                äquidistant ist. Der Lösungsaufwand würde sich etwa halbieren!!                           !
  !              - "r == 1 .AND. init_yes" sollte idealerweise aus den Schleifen herausgezogen werden, was   !
  !                hier aber aus Gründen der Übersicht bisher nicht ausgeführt wurde.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  omega = 1.
  
  
  IF (Jacobi_yes) THEN
  
  DO r = 1, n_relax
     
     IF (r == 1 .AND. init_yes) THEN
        !=====================================================================================================
        IF (BC_1L .GT. 0) THEN
           i = 1
           DO k = S33R, N33R
!pgi$ unroll = n:8
              DO j = S22R, N22R
                 rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L .GT. 0) THEN
           j = 1
           DO k = S33R, N33R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L .GT. 0 .AND. dimens == 3) THEN
           k = 1
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
              END DO
           END DO
        END IF
        !=====================================================================================================
        DO k = S33R, N33R
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 rel(i,j,k) = omega*bb(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
              END DO
           END DO
        END DO
        !=====================================================================================================
        IF (BC_1U .GT. 0) THEN
           i = N1
           DO k = S33R, N33R
!pgi$ unroll = n:8
              DO j = S22R, N22R
                 rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2U .GT. 0) THEN
           j = N2
           DO k = S33R, N33R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U .GT. 0 .AND. dimens == 3) THEN
           k = N3
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
              END DO
           END DO
        END IF
        !=====================================================================================================
        
     ELSE
        
        CALL exchange_relax(g,0,0,0,0,.TRUE.,rel)
        
        !=====================================================================================================
        IF (BC_1L .GT. 0) THEN
           i = 1
           DO k = S33R, N33R
!pgi$ unroll = n:8
              DO j = S22R, N22R
                 comp(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L .GT. 0) THEN
           j = 1
           DO k = S33R, N33R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 comp(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L .GT. 0 .AND. dimens == 3) THEN
           k = 1
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 comp(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
              END DO
           END DO
        END IF
        !=====================================================================================================
        DO k = S33R, N33R
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 comp(i,j,k) = omega*(bb(i,j,k)                                                &
                             &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                             &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                             &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                             &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
              END DO
           END DO
        END DO
        !=====================================================================================================
        IF (BC_1U .GT. 0) THEN
           i = N1
           DO k = S33R, N33R
!pgi$ unroll = n:8
              DO j = S22R, N22R
                 comp(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2U .GT. 0) THEN
           j = N2
           DO k = S33R, N33R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 comp(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
              END DO
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U .GT. 0 .AND. dimens == 3) THEN
           k = N3
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 comp(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
              END DO
           END DO
        END IF
        !=====================================================================================================
        
        rel(1:N1,1:N2,1:N3) = comp(1:N1,1:N2,1:N3)
        
     END IF
     
  END DO
  
  
  !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************
  ELSE
  
  
  IF (init_yes) THEN
     IF (BC_1L .LE. 0) rel(S11R-1 ,S2R:N2R,S3R:N3R) = 0.
     IF (BC_1U .LE. 0) rel(N11R+1 ,S2R:N2R,S3R:N3R) = 0.
     IF (BC_2L .LE. 0) rel(S1R:N1R,S22R-1 ,S3R:N3R) = 0.
     IF (BC_2U .LE. 0) rel(S1R:N1R,N22R+1 ,S3R:N3R) = 0.
     IF (BC_3L .LE. 0) rel(S1R:N1R,S2R:N2R,S33R-1 ) = 0.
     IF (BC_3U .LE. 0) rel(S1R:N1R,S2R:N2R,N33R+1 ) = 0.
     
     !rel = 0. ! TEST!!!  Nur zum testen, falls Initialisierung unten fehlerhaft erscheint:
  END IF
  
  
  DO r = 1, n_relax
     
     IF (.NOT. (r == 1 .AND. init_yes)) CALL exchange_relax(g,0,0,0,0,.TRUE.,rel)
     
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     IF (impl_dir(1) == 0 .AND. impl_dir(2) == 0 .AND. impl_dir(3) == 0) THEN
        
        IF (RBGS_mode .GT. 0) THEN
           !-------------------------------------------------------------------------------------------------!
           ! Hinweis:                                                                                        !
           ! Die Red-Black-Sortierung sollte für eine gute Performance global sein, d.h. ueber die Blöcke    !
           ! hinweg implementiert sein. Die aktuelle Umsetzung macht sich zunutze, dass im Parallelbetrieb   !
           ! die Anzahl der Gitterpunkte in jeder Raumrichtung automatisch eine gerade Zahl ist (die Tat-    !
           ! sache, dass in Randblöcken ggf. eine ungerade Anzahl Gitterpunkte abgearbeitet wird hat damit   !
           ! nichts zu tun!). Daher kann sich der Shift "ss" in jedem Block einfach auf den lokalen Gitter-  !
           ! punkt rel(1,1,1) beziehen.                                                                      !
           !                                                                                                 !
           ! Die aktuelle Schreibweise in den Schleifen ist äquivalent zu:                                   !
           ! DO k = S33R, N33R                                                                               !
           !    DO j = S22R, N22R                                                                            !
           !       DO i = S11R, N11R                                                                         !
           !          IF (MOD(i+j+k) == 1) ... ! (red)                                                       !
           !          ! bzw.                                                                                 !
           !          IF (MOD(i+j+k) == 0) ... ! (black)                                                     !
           !       END DO                                                                                    !
           !    END DO                                                                                       !
           ! END DO                                                                                          !
           !-------------------------------------------------------------------------------------------------!
           
           !==================================================================================================
           !==================================================================================================
           IF (BC_1L .GT. 0) THEN
              i = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
                    ss = MOD(i+k+S22R,2)
!pgi$ unroll = n:8
                    DO j = ss+S22R, N22R, 2
                       rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
                    ss = MOD(i+k+S22R,2)
!pgi$ unroll = n:8
                    DO j = ss+S22R, N22R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_1U .GT. 0) THEN
              i = N1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
                    ss = MOD(i+k+S22R,2)
!pgi$ unroll = n:8
                    DO j = ss+S22R, N22R, 2
                       rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
                    ss = MOD(i+k+S22R,2)
!pgi$ unroll = n:8
                    DO j = ss+S22R, N22R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !==================================================================================================
           IF (BC_2L .GT. 0) THEN
              j = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_2U .GT. 0) THEN
              j = N2
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !==================================================================================================
           IF (RBGS_mode == 1) THEN
           !--------------------------------------------------------------------------------------------------
           IF (BC_3L .GT. 0) THEN
              k = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                    END DO
                 END DO
              ELSE
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_3U .GT. 0) THEN
              k = N3
              IF (r == 1 .AND. init_yes) THEN
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                    END DO
                 END DO
              ELSE
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           ELSE
           !--------------------------------------------------------------------------------------------------
           IF (BC_3U .GT. 0) THEN
              k = N3
              IF (r == 1 .AND. init_yes) THEN
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                    END DO
                 END DO
              ELSE
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           END IF
           !==================================================================================================
           IF (RBGS_mode == 1) THEN
           !--------------------------------------------------------------------------------------------------
           IF (r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R,2) ! "+k" ist offenbar elementar für grosse Konvergenzrate. Symmetrie zu relaxation_div_grad_inv ist hier interessanterweise kontraproduktiv!
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*bb(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       ELSE
                          rel(i,j,k) =       bb(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
              END DO
              !CALL exchange_relax(g,0,0,0,0,.TRUE.,rel) ! TEST!!! Austausch bringt praktisch nichts, macht aber die Glättung unabhängig von der Parallelisierung (sonst teilweise Jacobi-Iteration)!
              DO k = S33R, N33R
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                       ELSE
                          rel(i,j,k) =       (bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                       ELSE
                          rel(i,j,k) =       (bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
              END DO
              !CALL exchange_relax(g,0,0,0,0,.TRUE.,rel) ! TEST!!! Austausch bringt praktisch nichts, macht aber die Glättung unabhängig von der Parallelisierung (sonst teilweise Jacobi-Iteration)!
              DO k = S33R, N33R
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                       ELSE
                          rel(i,j,k) =       (bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
              END DO
           END IF
           !--------------------------------------------------------------------------------------------------
           ELSE
           !--------------------------------------------------------------------------------------------------
           IF (r == 1 .AND. init_yes) THEN
              DO k = N33R, S33R, -1
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R,2) ! "+k" ist offenbar elementar für grosse Konvergenzrate. Symmetrie zu relaxation_div_grad_inv ist hier interessanterweise kontraproduktiv!
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       ELSE
                          rel(i,j,k) =       (bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &                                  - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                       ELSE
                          rel(i,j,k) =       (bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &                                  - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
              END DO
           ELSE
              DO k = N33R, S33R, -1
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                       ELSE
                          rel(i,j,k) =       (bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                       ELSE
                          rel(i,j,k) =       (bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
              END DO
           END IF
           !--------------------------------------------------------------------------------------------------
           END IF
           !==================================================================================================
           IF (BC_1L .GT. 0) THEN
              i = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
                    ss = MOD(i+k+S22R+1,2)
!pgi$ unroll = n:8
                    DO j = ss+S22R, N22R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
                    ss = MOD(i+k+S22R+1,2)
!pgi$ unroll = n:8
                    DO j = ss+S22R, N22R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_1U .GT. 0) THEN
              i = N1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
                    ss = MOD(i+k+S22R+1,2)
!pgi$ unroll = n:8
                    DO j = ss+S22R, N22R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
                    ss = MOD(i+k+S22R+1,2)
!pgi$ unroll = n:8
                    DO j = ss+S22R, N22R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !==================================================================================================
           IF (BC_2L .GT. 0) THEN
              j = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_2U .GT. 0) THEN
              j = N2
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !==================================================================================================
           IF (RBGS_mode == 1) THEN
           !--------------------------------------------------------------------------------------------------
           IF (BC_3L .GT. 0) THEN
              k = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g)
                    END DO
                 END DO
              ELSE
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_3U .GT. 0) THEN
              k = N3
              IF (r == 1 .AND. init_yes) THEN
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g)
                    END DO
                 END DO
              ELSE
                 DO j = S22R, N22R
                    ss = MOD(j+k+S11R+1,2)
!pgi$ unroll = n:8
                    DO i = ss+S11R, N11R, 2
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           ELSE
           !--------------------------------------------------------------------------------------------------
           IF (BC_3L .GT. 0) THEN
              k = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g)
                    END DO
                 END DO
              ELSE
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           END IF
           !==================================================================================================
           !==================================================================================================
           
        ELSE
           
           !==================================================================================================
           !==================================================================================================
           IF (BC_1U .GT. 0) THEN
              i = N1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
!pgi$ unroll = n:8
                    DO j = S22R, N22R
                       rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
!pgi$ unroll = n:8
                    DO j = S22R, N22R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_2U .GT. 0) THEN
              j = N2
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_3U .GT. 0) THEN
              k = N3
              IF (r == 1 .AND. init_yes) THEN
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                    END DO
                 END DO
              ELSE
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !==================================================================================================
           IF (r == 1 .AND. init_yes) THEN
              DO k = N33R, S33R, -1
                 DO j = N22R, S22R, -1
!pgi$ unroll = n:8
                    DO i = N11R, S11R, -1
                       IF (SOR_yes) THEN
                          rel(i,j,k) = bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k) - cdg2( 1,j,g)*rel(i,j+1,k) - cdg3( 1,k,g)*rel(i,j,k+1)
                          rel(i,j,k) = omega*rel(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       ELSE
                          rel(i,j,k) = bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k) - cdg2( 1,j,g)*rel(i,j+1,k) - cdg3( 1,k,g)*rel(i,j,k+1)
                          rel(i,j,k) = rel(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
              END DO
           ELSE
              DO k = N33R, S33R, -1
                 DO j = N22R, S22R, -1
!pgi$ unroll = n:8
                    DO i = N11R, S11R, -1
                       IF (SOR_yes) THEN
                          rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                       ELSE
                          rel(i,j,k) =       (bb(i,j,k)                                                 &
                                      &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                      &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                      &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                      &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                       END IF
                    END DO
                 END DO
              END DO
           END IF
           !==================================================================================================
           IF (BC_1L .GT. 0) THEN
              i = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
!pgi$ unroll = n:8
                    DO j = S22R, N22R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
!pgi$ unroll = n:8
                    DO j = S22R, N22R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_2L .GT. 0) THEN
              j = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO k = S33R, N33R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g)
                    END DO
                 END DO
              ELSE
                 DO k = S33R, N33R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !--------------------------------------------------------------------------------------------------
           IF (BC_3L .GT. 0) THEN
              k = 1
              IF (r == 1 .AND. init_yes) THEN
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g)
                    END DO
                 END DO
              ELSE
                 DO j = S22R, N22R
!pgi$ unroll = n:8
                    DO i = S11R, N11R
                       rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                    END DO
                 END DO
              END IF
           END IF
           !==================================================================================================
           !==================================================================================================
        END IF
        
     END IF
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     IF (impl_dir(3) == 1) THEN
        
        IF (.NOT. (r == 1 .AND. init_yes)) CALL exchange_relax(g,0,0,0,0,.TRUE.,rel)
        !=====================================================================================================
        IF (BC_1U .GT. 0) THEN
           i = N1
           IF (r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO j = S22R, N22R
                    rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO j = S22R, N22R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2U .GT. 0) THEN
           j = N2
           IF (r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
        DO j = N22R, S22R, -1
           DO i = N11R, S11R, -1
              
              !--- Diagonalelement / rechte Seite aufbauen ---------------------------------------------------
              IF (r == 1 .AND. init_yes) THEN
!pgi$ unroll = n:8
                 DO k = S3R, N3R
                    band3(1,k) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                    band3(2,k) = bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)
                 END DO
              ELSE
!pgi$ unroll = n:8
                 DO k = S3R, N3R
                    band3(1,k) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                    band3(2,k) = bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)   &
                                    &      - cdg1( 1,i,g)*rel(i+1,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)
                 END DO
              END IF
              
              !--- SOR (1) -----------------------------------------------------------------------------------
              IF (SOR_yes .AND. .NOT. (r == 1 .AND. init_yes)) THEN
!pgi$ unroll = n:8
                 DO k = S3R, N3R
                    SOR3(k) = rel(i,j,k)
                 END DO
              END IF
              
              
              !--- mit RB überschreiben ----------------------------------------------------------------------
              IF (BC_3L .GT. 0) THEN
                 band3(1,S3R) = cdg3(0,S3R,g)
                 band3(2,S3R) = bb  (i,j,S3R)
              ELSE
                 band3(2,S3R) = band3(2,S3R) - cdg3(-1,S3R,g)*rel(i,j,S3R-1)
              END IF
              
              IF (BC_3U .GT. 0) THEN
                 band3(1,N3R) = cdg3(0,N3R,g)
                 band3(2,N3R) = bb  (i,j,N3R)
              ELSE
                 band3(2,N3R) = band3(2,N3R) - cdg3( 1,N3R,g)*rel(i,j,N3R+1)
              END IF
              
              
              !--- Gauss-Elimination hoch --------------------------------------------------------------------
!pgi$ unroll = n:8
              DO k = S3R+1, N3R
                 mult       = cdg3(-1,k,g) / band3(1,k-1)
                 band3(1,k) = band3(1,k) - mult * cdg3(1,k-1,g)
                 band3(2,k) = band3(2,k) - mult * band3(2,k-1)
              END DO
              
              
              !--- Gauss-Elimination runter ------------------------------------------------------------------
              rel(i,j,N3R) = band3(2,N3R) / band3(1,N3R)
!pgi$ unroll = n:8
              DO k = N3R-1, S3R, -1
                 rel(i,j,k) = (band3(2,k) - cdg3(1,k,g)*rel(i,j,k+1)) / band3(1,k)
              END DO
              
              
              !--- SOR (2) -----------------------------------------------------------------------------------
              IF (SOR_yes) THEN
                 IF (r == 1 .AND. init_yes) THEN
!pgi$ unroll = n:8
                    DO k = S3R, N3R
                       rel(i,j,k) = omega*rel(i,j,k)
                    END DO
                 ELSE
!pgi$ unroll = n:8
                    DO k = S3R, N3R
                       rel(i,j,k) = omega*rel(i,j,k) + (1.-omega)*SOR3(k)
                    END DO
                 END IF
              END IF
              
           END DO
        END DO
        !=====================================================================================================
        IF (BC_1L .GT. 0) THEN
           i = 1
           IF (r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO j = S22R, N22R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g)
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO j = S22R, N22R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L .GT. 0) THEN
           j = 1
           IF (r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g)
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
        
     END IF
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     IF (impl_dir(2) == 1) THEN
        
        IF (.NOT. (impl_dir(3) == 0 .AND. r == 1 .AND. init_yes)) CALL exchange_relax(g,0,0,0,0,.TRUE.,rel)
        !=====================================================================================================
        IF (BC_1U .GT. 0) THEN
           i = N1
           IF (impl_dir(3) == 0 .AND. r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO j = S22R, N22R
                    rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO j = S22R, N22R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U .GT. 0) THEN
           k = N3
           IF (impl_dir(3) == 0 .AND. r == 1 .AND. init_yes) THEN
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                 END DO
              END DO
           ELSE
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
        DO k = N33R, S33R, -1
           DO i = N11R, S11R, -1
              
              !--- Diagonalelement / rechte Seite aufbauen ---------------------------------------------------
              IF (impl_dir(3) == 0 .AND. r == 1 .AND. init_yes) THEN
!pgi$ unroll = n:8
                 DO j = S2R, N2R
                    band2(1,j) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                    band2(2,j) = bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)
                 END DO
              ELSE
!pgi$ unroll = n:8
                 DO j = S2R, N2R
                    band2(1,j) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                    band2(2,j) = bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)   &
                                    &      - cdg1( 1,i,g)*rel(i+1,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)
                 END DO
              END IF
              
              !--- SOR (1) -----------------------------------------------------------------------------------
              IF (SOR_yes .AND. .NOT. (impl_dir(3) == 0 .AND. r == 1 .AND. init_yes)) THEN
!pgi$ unroll = n:8
                 DO j = S2R, N2R
                    SOR2(j) = rel(i,j,k)
                 END DO
              END IF
              
              
              !--- mit RB überschreiben ----------------------------------------------------------------------
              IF (BC_2L .GT. 0) THEN
                 band2(1,S2R) = cdg2(0,S2R,g)
                 band2(2,S2R) = bb  (i,S2R,k)
              ELSE
                 band2(2,S2R) = band2(2,S2R) - cdg2(-1,S2R,g)*rel(i,S2R-1,k)
              END IF
              
              IF (BC_2U .GT. 0) THEN
                 band2(1,N2R) = cdg2(0,N2R,g)
                 band2(2,N2R) = bb  (i,N2R,k)
              ELSE
                 band2(2,N2R) = band2(2,N2R) - cdg2( 1,N2R,g)*rel(i,N2R+1,k)
              END IF
              
              
              !--- Gauss-Elimination hoch --------------------------------------------------------------------
!pgi$ unroll = n:8
              DO j = S2R+1, N2R
                 mult       = cdg2(-1,j,g) / band2(1,j-1)
                 band2(1,j) = band2(1,j) - mult * cdg2(1,j-1,g)
                 band2(2,j) = band2(2,j) - mult * band2(2,j-1)
              END DO
              
              
              !--- Gauss-Elimination runter ------------------------------------------------------------------
              rel(i,N2R,k) = band2(2,N2R) / band2(1,N2R)
!pgi$ unroll = n:8
              DO j = N2R-1, S2R, -1
                 rel(i,j,k) = (band2(2,j) - cdg2(1,j,g)*rel(i,j+1,k)) / band2(1,j)
              END DO
              
              
              !--- SOR (2) -----------------------------------------------------------------------------------
              IF (SOR_yes) THEN
                 IF (impl_dir(3) == 0 .AND. r == 1 .AND. init_yes) THEN
!pgi$ unroll = n:8
                    DO j = S2R, N2R
                       rel(i,j,k) = omega*rel(i,j,k)
                    END DO
                 ELSE
!pgi$ unroll = n:8
                    DO j = S2R, N2R
                       rel(i,j,k) = omega*rel(i,j,k) + (1.-omega)*SOR2(j)
                    END DO
                 END IF
              END IF
              
           END DO
        END DO
        !=====================================================================================================
        IF (BC_1L .GT. 0) THEN
           i = 1
           IF (impl_dir(3) == 0 .AND. r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO j = S22R, N22R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g)
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO j = S22R, N22R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L .GT. 0) THEN
           k = 1
           IF (impl_dir(3) == 0 .AND. r == 1 .AND. init_yes) THEN
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g)
                 END DO
              END DO
           ELSE
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
        
     END IF
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     IF (impl_dir(1) == 1) THEN
        
        !=====================================================================================================
        IF (BC_2U .GT. 0) THEN
           j = N2
           IF (impl_dir(2) == 0 .AND. impl_dir(3) == 0 .AND. r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3U .GT. 0) THEN
           k = N3
           IF (impl_dir(2) == 0 .AND. impl_dir(3) == 0 .AND. r == 1 .AND. init_yes) THEN
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                 END DO
              END DO
           ELSE
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
        DO k = N33R, S33R, -1
           DO j = N22R, S22R, -1
              
              !--- Diagonalelement / rechte Seite aufbauen ---------------------------------------------------
              IF (impl_dir(2) == 0 .AND. impl_dir(3) == 0 .AND. r == 1 .AND. init_yes) THEN
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    band1(1,i) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                    band1(2,i) = bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k) - cdg3( 1,k,g)*rel(i,j,k+1)
                 END DO
              ELSE
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    band1(1,i) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                    band1(2,i) = bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k) - cdg3(-1,k,g)*rel(i,j,k-1)   &
                                    &      - cdg2( 1,j,g)*rel(i,j+1,k) - cdg3( 1,k,g)*rel(i,j,k+1)
                 END DO
              END IF
              
              !--- SOR (1) -----------------------------------------------------------------------------------
              IF (SOR_yes .AND. .NOT. (impl_dir(2) == 0 .AND. impl_dir(3) == 0 .AND. r == 1 .AND. init_yes)) THEN
!pgi$ unroll = n:8
                 DO i = S1R, N1R
                    SOR1(i) = rel(i,j,k)
                 END DO
              END IF
              
              
              !--- mit RB überschreiben ----------------------------------------------------------------------
              IF (BC_1L .GT. 0) THEN
                 band1(1,S1R) = cdg1(0,S1R,g)
                 band1(2,S1R) = bb  (S1R,j,k)
              ELSE
                 band1(2,S1R) = band1(2,S1R) - cdg1(-1,S1R,g)*rel(S1R-1,j,k)
              END IF
              
              IF (BC_1U .GT. 0) THEN
                 band1(1,N1R) = cdg1(0,N1R,g)
                 band1(2,N1R) = bb  (N1R,j,k)
              ELSE
                 band1(2,N1R) = band1(2,N1R) - cdg1( 1,N1R,g)*rel(N1R+1,j,k)
              END IF
              
              
              !--- Gauss-Elimination hoch --------------------------------------------------------------------
!pgi$ unroll = n:8
              DO i = S1R+1, N1R
                 mult       = cdg1(-1,i,g) / band1(1,i-1)
                 band1(1,i) = band1(1,i) - mult * cdg1(1,i-1,g)
                 band1(2,i) = band1(2,i) - mult * band1(2,i-1)
              END DO
              
              
              !--- Gauss-Elimination runter ------------------------------------------------------------------
              rel(N1R,j,k) = band1(2,N1R) / band1(1,N1R)
!pgi$ unroll = n:8
              DO i = N1R-1, S1R, -1
                 rel(i,j,k) = (band1(2,i) - cdg1(1,i,g)*rel(i+1,j,k)) / band1(1,i)
              END DO
              
              
              !--- SOR (2) -----------------------------------------------------------------------------------
              IF (SOR_yes) THEN
                 IF (impl_dir(2) == 0 .AND. impl_dir(3) == 0 .AND. r == 1 .AND. init_yes) THEN
!pgi$ unroll = n:8
                    DO i = S1R, N1R
                       rel(i,j,k) = omega*rel(i,j,k)
                    END DO
                 ELSE
!pgi$ unroll = n:8
                    DO i = S1R, N1R
                       rel(i,j,k) = omega*rel(i,j,k) + (1.-omega)*SOR1(i)
                    END DO
                 END IF
              END IF
              
           END DO
        END DO
        !=====================================================================================================
        IF (BC_2L .GT. 0) THEN
           j = 1
           IF (impl_dir(2) == 0 .AND. impl_dir(3) == 0 .AND. r == 1 .AND. init_yes) THEN
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g)
                 END DO
              END DO
           ELSE
              DO k = S33R, N33R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BC_3L .GT. 0) THEN
           k = 1
           IF (impl_dir(2) == 0 .AND. impl_dir(3) == 0 .AND. r == 1 .AND. init_yes) THEN
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g)
                 END DO
              END DO
           ELSE
              DO j = S22R, N22R
!pgi$ unroll = n:8
                 DO i = S11R, N11R
                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                 END DO
              END DO
           END IF
        END IF
        !=====================================================================================================
     END IF
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     
  END DO
  !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************
  
  END IF
  
  
  IF (corner_yes) CALL handle_corner_Lap(g,rel)
  
  
  END SUBROUTINE relaxation_div_grad_inv
  
  
  
  
  
  
  
  
  
  
  !> sets the pressure to zero at corners and edges to not create a residuum in the iterative solver
  SUBROUTINE handle_corner_Lap(g,phi)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  g
  REAL   , INTENT(inout) ::  phi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Diese Routine dient dazu, den unbestimmten Druck in den Ecken und Kanten explizit zu      !
  !                behandeln und gleichzeitig das Konvergenzverhalten der Löser (BiCGstab oder Richardson)   !
  !                möglichst nicht zu beeinträchtigen.                                                       !
  !              - Der Druck wird hier direkt zu Null gesetzt, um innerhalb des iterativen Lösers kein Re-   !
  !                siduum zu erzeugen (anstelle über die RHS).                                               !
  !              - Der Druck wird erst bei Bedarf (z.B. vor einem Ausschrieb) auf einen sinnvollen Wert ge-  !
  !                setzt.                                                                                    !
  !              - Siehe dazu auch die korrespondierende Subroutine "handle_corner_rhs"!                     !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (BC(1,1,g) .GT. 0 .AND. BC(1,2,g) .GT. 0) phi(1      ,1      ,1:NN(3,g)) = 0. ! TEST!!! verifizieren ...
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
  
  
  END SUBROUTINE handle_corner_Lap
  
  
  
  
END MODULE mod_laplace
