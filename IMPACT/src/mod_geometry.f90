!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!*************************************************************************************************************

!> @file mod_geometry.f90
!! File holding the module mod_geometry.

!> module containing subroutine coordinates. It uses modules mod_dims, mod_vars.
MODULE mod_geometry
  
  
  USE mod_dims
  USE mod_vars
  
  
  PRIVATE
  
  PUBLIC coordinates
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> Subroutine that gets the coordinates by calling usr_geometry::get_coords. For every spacial dimension
  !! it sets up the coordinates with the respective boundary conditions globally and then distributes them
  !! onto the blocks.
  SUBROUTINE coordinates
  
  IMPLICIT NONE
  
  INTEGER               ::  g
  INTEGER               ::  Nmax, Mmax
  
  INTEGER               ::  i, ii, di, iiShift
  INTEGER               ::  j, jj, dj, jjShift
  INTEGER               ::  k, kk, dk, kkShift
  
  REAL                  ::  y1pR(b1L:(M1+b1U)), y1uR(b1L:(M1+b1U))
  REAL                  ::  y2pR(b2L:(M2+b2U)), y2vR(b2L:(M2+b2U))
  REAL                  ::  y3pR(b3L:(M3+b3U)), y3wR(b3L:(M3+b3U))
  
  CHARACTER(LEN=1)      ::  part, grid
  
  
  x1p   = 0.
  x2p   = 0.
  x3p   = 0.
  
  x1u   = 0.
  x2v   = 0.
  x3w   = 0.
  
  dx1p  = 0.
  dx2p  = 0.
  dx3p  = 0.
  
  dx1u  = 0.
  dx2v  = 0.
  dx3w  = 0.
  
  !=== Koordinaten bestimmen (global) ========================================================================
  CALL get_coords
  
  
  !=== grobe Gitter (global) =================================================================================
  DO g = 1, n_grids
     
     Nmax = NN(1,g)
     Mmax = NB(1,g)*(Nmax-1)+1
     
     di   = (M1-1) / (Mmax-1)
     
     y1pR = 0.
     y1uR = 0.
     
     IF (g == 1) THEN
        y1pR(1:M1) = y1p(1:M1)
        y1uR(0:M1) = y1u(0:M1)
     ELSE
        DO i = 1, Mmax
           y1pR(i) = y1p(1 + di*i - di  )
        END DO
        DO i = 1, Mmax-1
           y1uR(i) = y1p(1 + di*i - di/2)
        END DO
     END IF
     
     
     !--- periodische RB -------------------------------------------------------------------------------------
     IF (BC_1L_global == -1) THEN
        y1pR(b1L:0) = y1pR((Mmax-1+b1L):(Mmax-1)) - L1
        y1uR(b1L:0) = y1uR((Mmax-1+b1L):(Mmax-1)) - L1
        
        y1pR(Mmax:(Mmax+b1U)) = y1pR(1:(1+b1U)) + L1
        y1uR(Mmax:(Mmax+b1U)) = y1uR(1:(1+b1U)) + L1
     END IF
     
     
     !--- Symmetrie-RB oder feste Wände ----------------------------------------------------------------------
     IF (BC_1L_global == -2 .OR. BC_1L_global .GT. 0) THEN
        IF (BC_1L_global .GT. 0 .AND. g == 1) THEN
           y1pR(b1L: 0) = 2.*y1pR(1) - y1pR((2-b1L):2:-1)
           y1uR(b1L:-1) = 2.*y1pR(1) - y1uR((1-b1L):2:-1)
        ELSE
           y1pR(b1L: 0) = 2.*y1pR(1) - y1pR((2-b1L):2:-1)
           y1uR(b1L: 0) = 2.*y1pR(1) - y1uR((1-b1L):1:-1)
        END IF
     END IF
     
     IF (BC_1U_global == -2 .OR. BC_1U_global .GT. 0) THEN
        IF (BC_1L_global .GT. 0 .AND. g == 1) THEN
           y1pR((Mmax+1):(Mmax+b1U)) = 2.*y1pR(Mmax) - y1pR((Mmax-1):(Mmax  -b1U):-1)
           y1uR((Mmax+1):(Mmax+b1U)) = 2.*y1pR(Mmax) - y1uR((Mmax-2):(Mmax-1-b1U):-1)
        ELSE
           y1pR((Mmax+1):(Mmax+b1U)) = 2.*y1pR(Mmax) - y1pR((Mmax-1):(Mmax  -b1U):-1)
           y1uR( Mmax   :(Mmax+b1U)) = 2.*y1pR(Mmax) - y1uR((Mmax-1):(Mmax-1-b1U):-1)
        END IF
     END IF
     
     
     !--- Verteilen auf Blöcke -------------------------------------------------------------------------------
     iiShift = (iB(1,g)-1)*(Nmax-1)
     
     x1pR(b1L:(Nmax+b1U),g) = y1pR((iiShift+b1L):(iiShift+Nmax+b1U))
     x1uR(b1L:(Nmax+b1U),g) = y1uR((iiShift+b1L):(iiShift+Nmax+b1U))
     
     
     !--- Schreiben ------------------------------------------------------------------------------------------
     IF (rank == 0 .AND. write_test_yes) THEN
        
        WRITE(grid,'(i1.1)') g
        
        OPEN(10, FILE='test_coord_p1_grid'//grid//'_restart'//restart_char//'.txt', STATUS='UNKNOWN')
        OPEN(11, FILE='test_coord_u1_grid'//grid//'_restart'//restart_char//'.txt', STATUS='UNKNOWN')
        
        DO ii = b1L, Nmax+b1U
           i = 1 + di*(ii-1)
           WRITE(10,'(2E25.17)') REAL(i+iShift)             , x1pR(ii,g)
           WRITE(11,'(2E25.17)') REAL(i+iShift)+0.5*REAL(di), x1uR(ii,g)
        END DO
        
        CLOSE(10)
        CLOSE(11)
        
     END IF
     
  END DO
  
  x1p(:) = x1pR(:,1)
  x1u(:) = x1uR(:,1)
  
  ! Ist besser geeignet fuer Auswertung (Integrationsgewichte):
  DO i = 1, N1
     dx1p(i) = x1u(i)-x1u(i-1)
  END DO
  IF (BC_1L .GT. 0 .OR. BC_1L == -2) dx1p(1 ) = x1u(1 )-x1p(1   )
  IF (BC_1U .GT. 0 .OR. BC_1U == -2) dx1p(N1) = x1p(N1)-x1u(N1-1)
  
  DO i = 1, N1-1
     dx1u(i) = x1p(i+1)-x1p(i)
  END DO
  dx1u(0 ) = dx1u(1   ) ! Notwendig fuer Partikel-Rückkopplung, bbecsek 021214: TODO: remove if possible
  dx1u(N1) = dx1u(N1-1)
  IF (BC_1L == 0 .OR. BC_1L == -1) dx1u(0 ) = x1p(1   )-x1p(0 )
  IF (BC_1U == 0 .OR. BC_1U == -1) dx1u(N1) = x1p(N1+1)-x1p(N1)
 
  !===========================================================================================================
  
  
  
  !=== grobe Gitter (global) =================================================================================
  DO g = 1, n_grids
     
     Nmax = NN(2,g)
     Mmax = NB(2,g)*(Nmax-1)+1
     
     dj   = (M2-1) / (Mmax-1)
     
     y2pR = 0.
     y2vR = 0.
     
     IF (g == 1) THEN
        y2pR(1:M2) = y2p(1:M2)
        y2vR(0:M2) = y2v(0:M2)
     ELSE
        DO j = 1, Mmax
           y2pR(j) = y2p(1 + dj*j - dj  )
        END DO
        DO j = 1, Mmax-1
           y2vR(j) = y2p(1 + dj*j - dj/2)
        END DO
     END IF
     
     
     !--- periodische RB -------------------------------------------------------------------------------------
     IF (BC_2L_global == -1) THEN
        y2pR(b2L:0) = y2pR((Mmax-1+b2L):(Mmax-1)) - L2
        y2vR(b2L:0) = y2vR((Mmax-1+b2L):(Mmax-1)) - L2
        
        y2pR(Mmax:(Mmax+b2U)) = y2pR(1:(1+b2U)) + L2
        y2vR(Mmax:(Mmax+b2U)) = y2vR(1:(1+b2U)) + L2
     END IF
     
     
     !--- Symmetrie-RB oder feste Wände ----------------------------------------------------------------------
     IF (BC_2L_global == -2 .OR. BC_2L_global .GT. 0) THEN
        IF (BC_2L_global .GT. 0 .AND. g == 1) THEN
           y2pR(b2L: 0) = 2.*y2pR(1) - y2pR((2-b2L):2:-1)
           y2vR(b2L:-1) = 2.*y2pR(1) - y2vR((1-b2L):2:-1)
        ELSE
           y2pR(b2L: 0) = 2.*y2pR(1) - y2pR((2-b2L):2:-1)
           y2vR(b2L: 0) = 2.*y2pR(1) - y2vR((1-b2L):1:-1)
        END IF
     END IF
     
     IF (BC_2U_global == -2 .OR. BC_2U_global .GT. 0) THEN
        IF (BC_2U_global .GT. 0 .AND. g == 1) THEN
           y2pR((Mmax+1):(Mmax+b2U)) = 2.*y2pR(Mmax) - y2pR((Mmax-1):(Mmax  -b2U):-1)
           y2vR((Mmax+1):(Mmax+b2U)) = 2.*y2pR(Mmax) - y2vR((Mmax-2):(Mmax-1-b2U):-1)
        ELSE
           y2pR((Mmax+1):(Mmax+b2U)) = 2.*y2pR(Mmax) - y2pR((Mmax-1):(Mmax  -b2U):-1)
           y2vR( Mmax   :(Mmax+b2U)) = 2.*y2pR(Mmax) - y2vR((Mmax-1):(Mmax-1-b2U):-1)
        END IF
     END IF
     
     
     !--- Verteilen auf Blöcke -------------------------------------------------------------------------------
     jjShift = (iB(2,g)-1)*(Nmax-1)
     
     x2pR(b2L:(Nmax+b2U),g) = y2pR((jjShift+b2L):(jjShift+Nmax+b2U))
     x2vR(b2L:(Nmax+b2U),g) = y2vR((jjShift+b2L):(jjShift+Nmax+b2U))
     
     
     !--- Schreiben ------------------------------------------------------------------------------------------
     IF (rank == 0 .AND. write_test_yes) THEN
        
        WRITE(grid,'(i1.1)') g
        
        OPEN(10, FILE='test_coord_p2_grid'//grid//'_restart'//restart_char//'.txt', STATUS='UNKNOWN')
        OPEN(11, FILE='test_coord_v2_grid'//grid//'_restart'//restart_char//'.txt', STATUS='UNKNOWN')
        
        DO jj = b2L, Nmax+b2U
           j = 1 + dj*(jj-1)
           WRITE(10,'(2E25.17)') REAL(j+jShift)             , x2pR(jj,g)
           WRITE(11,'(2E25.17)') REAL(j+jShift)+0.5*REAL(dj), x2vR(jj,g)
        END DO
        
        CLOSE(10)
        CLOSE(11)
        
     END IF
     
  END DO
  
  x2p(:) = x2pR(:,1)
  x2v(:) = x2vR(:,1)
  
  DO j = 1, N2
     dx2p(j) = x2v(j)-x2v(j-1)
  END DO
  IF (BC_2L .GT. 0 .OR. BC_2L == -2) dx2p(1 ) = x2v(1 )-x2p(1   )
  IF (BC_2U .GT. 0 .OR. BC_2U == -2) dx2p(N2) = x2p(N2)-x2v(N2-1)
  
  DO j = 1, N2-1
     dx2v(j) = x2p(j+1)-x2p(j)
  END DO
  dx2v(0 ) = dx2v(1   ) ! Notwendig fuer Partikel-Rückkopplung, bbecsek 021214: TODO: remove if possible
  dx2v(N2) = dx2v(N2-1)
  IF (BC_2L == 0 .OR. BC_2L == -1) dx2v(0 ) = x2p(1   )-x2p(0 )
  IF (BC_2U == 0 .OR. BC_2U == -1) dx2v(N2) = x2p(N2+1)-x2p(N2)
  
  !===========================================================================================================
  
  
  IF (dimens == 3) THEN
  !=== grobe Gitter (global) =================================================================================
  DO g = 1, n_grids
     
     Nmax = NN(3,g)
     Mmax = NB(3,g)*(Nmax-1)+1
     
     dk   = (M3-1) / (Mmax-1)
     
     y3pR = 0.
     y3wR = 0.
     
     IF (g == 1) THEN
        y3pR(1:M3) = y3p(1:M3)
        y3wR(0:M3) = y3w(0:M3)
     ELSE
        DO k = 1, Mmax
           y3pR(k) = y3p(1 + dk*k - dk  )
        END DO
        DO k = 1, Mmax-1
           y3wR(k) = y3p(1 + dk*k - dk/2)
        END DO
     END IF
     
     
     !--- periodische RB -------------------------------------------------------------------------------------
     IF (BC_3L_global == -1) THEN
        y3pR(b3L:0) = y3pR((Mmax-1+b3L):(Mmax-1)) - L3
        y3wR(b3L:0) = y3wR((Mmax-1+b3L):(Mmax-1)) - L3
        
        y3pR(Mmax:(Mmax+b3U)) = y3pR(1:(1+b3U)) + L3
        y3wR(Mmax:(Mmax+b3U)) = y3wR(1:(1+b3U)) + L3
     END IF
     
     
     !--- Symmetrie-RB oder feste Wände ----------------------------------------------------------------------
     IF (BC_3L_global == -2 .OR. BC_3L_global .GT. 0) THEN
        IF (BC_3L_global .GT. 0 .AND. g == 1) THEN
           y3pR(b3L: 0) = 2.*y3pR(1) - y3pR((2-b3L):2:-1)
           y3wR(b3L:-1) = 2.*y3pR(1) - y3wR((1-b3L):2:-1)
        ELSE
           y3pR(b3L: 0) = 2.*y3pR(1) - y3pR((2-b3L):2:-1)
           y3wR(b3L: 0) = 2.*y3pR(1) - y3wR((1-b3L):1:-1)
        END IF
     END IF
     
     IF (BC_3U_global == -2 .OR. BC_3U_global .GT. 0) THEN
        IF (BC_3U_global .GT. 0 .AND. g == 1) THEN
           y3pR((Mmax+1):(Mmax+b3U)) = 2.*y3pR(Mmax) - y3pR((Mmax-1):(Mmax  -b3U):-1)
           y3wR((Mmax+1):(Mmax+b3U)) = 2.*y3pR(Mmax) - y3wR((Mmax-2):(Mmax-1-b3U):-1)
        ELSE
           y3pR((Mmax+1):(Mmax+b3U)) = 2.*y3pR(Mmax) - y3pR((Mmax-1):(Mmax  -b3U):-1)
           y3wR( Mmax   :(Mmax+b3U)) = 2.*y3pR(Mmax) - y3wR((Mmax-1):(Mmax-1-b3U):-1)
        END IF
     END IF
     
     
     !--- Verteilen auf Blöcke -------------------------------------------------------------------------------
     kkShift = (iB(3,g)-1)*(Nmax-1)
     
     x3pR(b3L:(Nmax+b3U),g) = y3pR((kkShift+b3L):(kkShift+Nmax+b3U))
     x3wR(b3L:(Nmax+b3U),g) = y3wR((kkShift+b3L):(kkShift+Nmax+b3U))
     
     
     !--- Schreiben ------------------------------------------------------------------------------------------
     IF (rank == 0 .AND. write_test_yes) THEN
        
        WRITE(grid,'(i1.1)') g
        
        OPEN(10, FILE='test_coord_p3_grid'//grid//'_restart'//restart_char//'.txt', STATUS='UNKNOWN')
        OPEN(11, FILE='test_coord_w3_grid'//grid//'_restart'//restart_char//'.txt', STATUS='UNKNOWN')
        
        DO kk = b3L, Nmax+b3U
           k = 1 + dk*(kk-1)
           WRITE(10,'(2E25.17)') REAL(k+kShift)             , x3pR(kk,g)
           WRITE(11,'(2E25.17)') REAL(k+kShift)+0.5*REAL(dk), x3wR(kk,g)
        END DO
        
        CLOSE(10)
        CLOSE(11)
        
     END IF
     
  END DO
  
  x3p(:) = x3pR(:,1)
  x3w(:) = x3wR(:,1)
  
  DO k = 1, N3
     dx3p(k) = x3w(k)-x3w(k-1)
  END DO
  IF (BC_3L .GT. 0 .OR. BC_3L == -2) dx3p(1 ) = x3w(1 )-x3p(1   )
  IF (BC_3U .GT. 0 .OR. BC_3U == -2) dx3p(N3) = x3p(N3)-x3w(N3-1)
  
  DO k = 1, N3-1
     dx3w(k) = x3p(k+1)-x3p(k)
  END DO
  dx3w(0 ) = dx3w(1   ) ! Notwendig fuer Partikel-Rückkopplung, bbecsek 021214: TODO: remove if possible
  dx3w(N3) = dx3w(N3-1)
  IF (BC_3L == 0 .OR. BC_3L == -1) dx3w(0 ) = x3p(1   )-x3p(0 )
  IF (BC_3U == 0 .OR. BC_3U == -1) dx3w(N3) = x3p(N3+1)-x3p(N3)
  
  !===========================================================================================================
  ELSE
  ! fuer Statistiken: ! TEST!!! ok?
  dy3p  = 1.
  dy3w  = 1.
  dx3p  = 1.
  dx3w  = 1.

  END IF
  
  
  END SUBROUTINE coordinates
  
  
  
  
END MODULE mod_geometry
