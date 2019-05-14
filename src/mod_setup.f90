!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!*************************************************************************************************************

!> @file mod_setup.f90
!! File containing mod_setup with initialization routines

!> Module containing public initialization routines for discretization parameters, parallel computation setup,
!! boundary setup, and index setup. It uses mod_dims, mod_vars, mod_lib.
MODULE mod_setup
  
  
  USE mod_dims
  USE mod_vars
  USE mod_lib ! (num_to_string)
  !USE mpi

  PRIVATE
  
  PUBLIC init_general, init_parallel, init_boundaries, init_limits
  
  INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> Subroutine initializing discretization parameters (temporal and spatial), allocating fields and defining
  !! output.
  SUBROUTINE init_general
  
  IMPLICIT NONE
  
  INTEGER                ::  b, g, m
  INTEGER                ::  stride(1:3)
  INTEGER                ::  MM(1:3,1:n_grids_max)
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Da bei den Helmholtz-Problemen für die Geschwindigkeiten auf dem feinsten Gitter          !
  !                mindestens soviele Punkte in Wand-normaler Richtung vorhanden sind wie bei dem Poisson-   !
  !                Problem bzw. die selben groben Gittern verwendet werden, kann es hier zu keinem Konflikt  !
  !                kommen.                                                                                   !
  !----------------------------------------------------------------------------------------------------------!
  
  ! bbecsek 021214: TODO: timeint_mode is changed here. Add message that informs of the change. Perhaps add
  ! 			  a logical switch that only executes when timein_mode is in the wrong state.
  !===========================================================================================================
  !=== Zeitintegration =======================================================================================
  !===========================================================================================================
  ! - Euler-explizit ist nicht erlaubt.
  ! - Euler_yes (== Re --> infty) ==> ist immer explizit ==> timeint_mode = 1.
  ! - thetaL ist nur bei timeint_mode == 0 von Bedeutung.
  
  IF (Euler_yes) timeint_mode = 1
  !===========================================================================================================
  
  ! bbecsek 021214: TODO: this command block is previously already executed in mod_test::test_parameter
  !===========================================================================================================
  !=== Dimension =============================================================================================
  !===========================================================================================================
  IF (M3 == 2) THEN
     dimens = 2
  ELSE
     dimens = 3
  END IF
  !===========================================================================================================
  
  
#ifdef ALLOC
  !===========================================================================================================
  !=== Domain- und Blockspezifikationen ======================================================================
  !===========================================================================================================
  N1 = 1+(M1-1)/NB1
  N2 = 1+(M2-1)/NB2
  N3 = 1+(M3-1)/NB3
  !===========================================================================================================
#endif

  !===========================================================================================================
  !=== Multigrid-Gitterspezifikationen =======================================================================
  !===========================================================================================================
  ! TEST!!! Test schreiben, dass Preconditioning von Helmholtz-Gleichungen (vel+conc) momentan nicht voll unterstützt ist!
  
  MM(1:3,1) = (/M1 ,M2 ,M3 /)
  NN(1:3,1) = (/N1 ,N2 ,N3 /)
  NB(1:3,1) = (/NB1,NB2,NB3/)
  
  n_grids  = 1
  n_gather = 1
  stride   = 1
  
  ! (alte Version geloescht am 24.01.2010)
  DO g = 2, n_grids_max
     DO m = 1, 3
        !--- global number of grid points ---
        IF (MOD(MM(m,g-1)-1,2) == 0 .AND. MM(m,g-1) .GE. 5) THEN ! 5=2*(Nmin-1)+1, Nmin=3 for Laplacian
           MM(m,g) = (MM(m,g-1)-1) / 2 + 1
           n_grids = g
        ELSE
           MM(m,g) = MM(m,g-1)
        END IF
        
        !--- find largest possible number of blocks ---
        NB(m,g) = 1
        DO b = 2, NB(m,g-1)
           IF (MOD(MM(m,g)-1,b) == 0 .AND. MOD(NB(m,g-1),b) == 0) THEN
              ! prolongation to next coarser grid must be possible on each block, independent of the other blocks!
              ! NN(m,g+1)-1 .GE. 1 (next coarser grid)
              IF (MOD((MM(m,g)-1)/b,2) == 0 .AND. (MM(m,g)-1)/b .GE. 2) NB(m,g) = b
           END IF
        END DO
        
        !--- enforce gathering of coarsest grid to one processor ---
        IF ((MOD(MM(m,g)-1,2) /= 0 .OR. MM(m,g) .LT. 5) .AND. NB(m,g) .GT. 1) NB(m,g) = 1
        
        !--- number of gathered blocks ---
        n_gather(m,g) = NB(m,g-1)/NB(m,g)
        
        !--- local grid size ---
        NN(m,g) = (MM(m,g)-1)/NB(m,g)+1
     END DO
     
  END DO
  
  
  IF (n_grids .GT. n_grids_limit) n_grids = n_grids_limit
  
  
  IF (rank == 0) THEN
     OPEN(10,FILE='test_multigrid_properties.txt', STATUS='UNKNOWN')
     WRITE(10,'(a     )') 'grid dimensions (total):'
     WRITE(10,'(a,20i6)') '         M1: ', (NN(1,1:n_grids)-1)*NB(1,1:n_grids)+1
     WRITE(10,'(a,20i6)') '         M2: ', (NN(2,1:n_grids)-1)*NB(2,1:n_grids)+1
     WRITE(10,'(a,20i6)') '         M3: ', (NN(3,1:n_grids)-1)*NB(3,1:n_grids)+1
     WRITE(10,*)
     WRITE(10,'(a     )') 'grid dimensions (per processor block):'
     WRITE(10,'(a,20i6)') '         N1: ', NN(1,1:n_grids)
     WRITE(10,'(a,20i6)') '         N2: ', NN(2,1:n_grids)
     WRITE(10,'(a,20i6)') '         N3: ', NN(3,1:n_grids)
     WRITE(10,*)
     WRITE(10,'(a     )') 'Processor blocks:'
     WRITE(10,'(a,20i6)') '        NB1: ', NB(1,1:n_grids)
     WRITE(10,'(a,20i6)') '        NB2: ', NB(2,1:n_grids)
     WRITE(10,'(a,20i6)') '        NB3: ', NB(3,1:n_grids)
     WRITE(10,*)
     WRITE(10,'(a     )') 'Processor block gathering:'
     WRITE(10,'(a,20i6)') 'n_gather(1): ', n_gather(1,1:n_grids)
     WRITE(10,'(a,20i6)') 'n_gather(2): ', n_gather(2,1:n_grids)
     WRITE(10,'(a,20i6)') 'n_gather(3): ', n_gather(3,1:n_grids)
     CLOSE(10)
  END IF
  
  
  !===========================================================================================================
  !=== Felder allocieren =====================================================================================
  !===========================================================================================================
  ! TEST!!! alloc.f90 ganz hineinziehen?
#ifdef ALLOC
   INCLUDE 'alloc.f90'

#ifndef FTOPY
  !--- C-FORTRAN interface ---
  dy1pptr = C_LOC(dy1p(1))
  dy2pptr = C_LOC(dy2p(1))
  dy3pptr = C_LOC(dy3p(1))

  dy1uptr = C_LOC(dy1u(0))
  dy2vptr = C_LOC(dy2v(0))
  dy3wptr = C_LOC(dy3w(0))

  dx1pptr = C_LOC(dx1p(1))
  dx2pptr = C_LOC(dx2p(1))
  dx3pptr = C_LOC(dx3p(1))

  dx1uptr = C_LOC(dx1u(0))
  dx2vptr = C_LOC(dx2v(0))
  dx3wptr = C_LOC(dx3w(0))

  velptr = C_LOC(vel(b1L,b2L,b3L,1))
  nlptr = C_LOC(nl(b1L,b2L,b3L,1))
  rhsptr = C_LOC(rhs(b1L,b2L,b3L,1))
  preptr = C_LOC(pre(b1L,b2L,b3L))

  work1ptr = C_LOC(work1(b1L,b2L,b3L))
  work2ptr = C_LOC(work2(b1L,b2L,b3L))
  work3ptr = C_LOC(work3(b1L,b2L,b3L))

  y1pptr = C_LOC(y1p(1))
  y2pptr = C_LOC(y2p(1))
  y3pptr = C_LOC(y3p(1))

  y1uptr = C_LOC(y1u(0))
  y2vptr = C_LOC(y2v(0))
  y3wptr = C_LOC(y3w(0))

  x1pptr = C_LOC(x1p(b1L))
  x2pptr = C_LOC(x2p(b2L))
  x3pptr = C_LOC(x3p(b3L))

  x1uptr = C_LOC(x1u(b1L))
  x2vptr = C_LOC(x2v(b2L))
  x3wptr = C_LOC(x3w(b3L))

  !=== IBM (bbecsek) ==========================================================================================
  !--- Force density ------------------------------------------------------------------------------------------
  fdptr = C_LOC(fd(b1L,b2L,b3L,1))
  vel_oldptr = C_LOC(vel_old(b1L,b2L,b3L,1))
#endif

#endif

  IF (n_grids .GE.  2) ALLOCATE(vec2A (b1L:(NN(1,2 )+b1U),b2L:(NN(2,2 )+b2U),b3L:(NN(3,2 )+b3U)))
  IF (n_grids .GE.  2) ALLOCATE(vec2B (b1L:(NN(1,2 )+b1U),b2L:(NN(2,2 )+b2U),b3L:(NN(3,2 )+b3U)))
  IF (n_grids .GE.  2) ALLOCATE(vec2C (b1L:(NN(1,2 )+b1U),b2L:(NN(2,2 )+b2U),b3L:(NN(3,2 )+b3U)))
  
  IF (n_grids .GE.  3) ALLOCATE(vec3A (b1L:(NN(1,3 )+b1U),b2L:(NN(2,3 )+b2U),b3L:(NN(3,3 )+b3U)))
  IF (n_grids .GE.  3) ALLOCATE(vec3B (b1L:(NN(1,3 )+b1U),b2L:(NN(2,3 )+b2U),b3L:(NN(3,3 )+b3U)))
  IF (n_grids .GE.  3) ALLOCATE(vec3C (b1L:(NN(1,3 )+b1U),b2L:(NN(2,3 )+b2U),b3L:(NN(3,3 )+b3U)))
  
  IF (n_grids .GE.  4) ALLOCATE(vec4A (b1L:(NN(1,4 )+b1U),b2L:(NN(2,4 )+b2U),b3L:(NN(3,4 )+b3U)))
  IF (n_grids .GE.  4) ALLOCATE(vec4B (b1L:(NN(1,4 )+b1U),b2L:(NN(2,4 )+b2U),b3L:(NN(3,4 )+b3U)))
  IF (n_grids .GE.  4) ALLOCATE(vec4C (b1L:(NN(1,4 )+b1U),b2L:(NN(2,4 )+b2U),b3L:(NN(3,4 )+b3U)))
  
  IF (n_grids .GE.  5) ALLOCATE(vec5A (b1L:(NN(1,5 )+b1U),b2L:(NN(2,5 )+b2U),b3L:(NN(3,5 )+b3U)))
  IF (n_grids .GE.  5) ALLOCATE(vec5B (b1L:(NN(1,5 )+b1U),b2L:(NN(2,5 )+b2U),b3L:(NN(3,5 )+b3U)))
  IF (n_grids .GE.  5) ALLOCATE(vec5C (b1L:(NN(1,5 )+b1U),b2L:(NN(2,5 )+b2U),b3L:(NN(3,5 )+b3U)))
  
  IF (n_grids .GE.  6) ALLOCATE(vec6A (b1L:(NN(1,6 )+b1U),b2L:(NN(2,6 )+b2U),b3L:(NN(3,6 )+b3U)))
  IF (n_grids .GE.  6) ALLOCATE(vec6B (b1L:(NN(1,6 )+b1U),b2L:(NN(2,6 )+b2U),b3L:(NN(3,6 )+b3U)))
  IF (n_grids .GE.  6) ALLOCATE(vec6C (b1L:(NN(1,6 )+b1U),b2L:(NN(2,6 )+b2U),b3L:(NN(3,6 )+b3U)))
  
  IF (n_grids .GE.  7) ALLOCATE(vec7A (b1L:(NN(1,7 )+b1U),b2L:(NN(2,7 )+b2U),b3L:(NN(3,7 )+b3U)))
  IF (n_grids .GE.  7) ALLOCATE(vec7B (b1L:(NN(1,7 )+b1U),b2L:(NN(2,7 )+b2U),b3L:(NN(3,7 )+b3U)))
  IF (n_grids .GE.  7) ALLOCATE(vec7C (b1L:(NN(1,7 )+b1U),b2L:(NN(2,7 )+b2U),b3L:(NN(3,7 )+b3U)))
  
  IF (n_grids .GE.  8) ALLOCATE(vec8A (b1L:(NN(1,8 )+b1U),b2L:(NN(2,8 )+b2U),b3L:(NN(3,8 )+b3U)))
  IF (n_grids .GE.  8) ALLOCATE(vec8B (b1L:(NN(1,8 )+b1U),b2L:(NN(2,8 )+b2U),b3L:(NN(3,8 )+b3U)))
  IF (n_grids .GE.  8) ALLOCATE(vec8C (b1L:(NN(1,8 )+b1U),b2L:(NN(2,8 )+b2U),b3L:(NN(3,8 )+b3U)))
  
  IF (n_grids .GE.  9) ALLOCATE(vec9A (b1L:(NN(1,9 )+b1U),b2L:(NN(2,9 )+b2U),b3L:(NN(3,9 )+b3U)))
  IF (n_grids .GE.  9) ALLOCATE(vec9B (b1L:(NN(1,9 )+b1U),b2L:(NN(2,9 )+b2U),b3L:(NN(3,9 )+b3U)))
  IF (n_grids .GE.  9) ALLOCATE(vec9C (b1L:(NN(1,9 )+b1U),b2L:(NN(2,9 )+b2U),b3L:(NN(3,9 )+b3U)))
  
  IF (n_grids .GE. 10) ALLOCATE(vec10A(b1L:(NN(1,10)+b1U),b2L:(NN(2,10)+b2U),b3L:(NN(3,10)+b3U)))
  IF (n_grids .GE. 10) ALLOCATE(vec10B(b1L:(NN(1,10)+b1U),b2L:(NN(2,10)+b2U),b3L:(NN(3,10)+b3U)))
  IF (n_grids .GE. 10) ALLOCATE(vec10C(b1L:(NN(1,10)+b1U),b2L:(NN(2,10)+b2U),b3L:(NN(3,10)+b3U)))
  
  IF (n_grids .GE. 11) ALLOCATE(vec11A(b1L:(NN(1,11)+b1U),b2L:(NN(2,11)+b2U),b3L:(NN(3,11)+b3U)))
  IF (n_grids .GE. 11) ALLOCATE(vec11B(b1L:(NN(1,11)+b1U),b2L:(NN(2,11)+b2U),b3L:(NN(3,11)+b3U)))
  IF (n_grids .GE. 11) ALLOCATE(vec11C(b1L:(NN(1,11)+b1U),b2L:(NN(2,11)+b2U),b3L:(NN(3,11)+b3U)))
  
  IF (n_grids .GE. 12) ALLOCATE(vec12A(b1L:(NN(1,12)+b1U),b2L:(NN(2,12)+b2U),b3L:(NN(3,12)+b3U)))
  IF (n_grids .GE. 12) ALLOCATE(vec12B(b1L:(NN(1,12)+b1U),b2L:(NN(2,12)+b2U),b3L:(NN(3,12)+b3U)))
  IF (n_grids .GE. 12) ALLOCATE(vec12C(b1L:(NN(1,12)+b1U),b2L:(NN(2,12)+b2U),b3L:(NN(3,12)+b3U)))
  
  IF (n_grids .GE. 13) ALLOCATE(vec13A(b1L:(NN(1,13)+b1U),b2L:(NN(2,13)+b2U),b3L:(NN(3,13)+b3U)))
  IF (n_grids .GE. 13) ALLOCATE(vec13B(b1L:(NN(1,13)+b1U),b2L:(NN(2,13)+b2U),b3L:(NN(3,13)+b3U)))
  IF (n_grids .GE. 13) ALLOCATE(vec13C(b1L:(NN(1,13)+b1U),b2L:(NN(2,13)+b2U),b3L:(NN(3,13)+b3U)))
  
  IF (n_grids .GE. 14) ALLOCATE(vec14A(b1L:(NN(1,14)+b1U),b2L:(NN(2,14)+b2U),b3L:(NN(3,14)+b3U)))
  IF (n_grids .GE. 14) ALLOCATE(vec14B(b1L:(NN(1,14)+b1U),b2L:(NN(2,14)+b2U),b3L:(NN(3,14)+b3U)))
  IF (n_grids .GE. 14) ALLOCATE(vec14C(b1L:(NN(1,14)+b1U),b2L:(NN(2,14)+b2U),b3L:(NN(3,14)+b3U)))
  
  IF (n_grids .GE. 15) ALLOCATE(vec15A(b1L:(NN(1,15)+b1U),b2L:(NN(2,15)+b2U),b3L:(NN(3,15)+b3U)))
  IF (n_grids .GE. 15) ALLOCATE(vec15B(b1L:(NN(1,15)+b1U),b2L:(NN(2,15)+b2U),b3L:(NN(3,15)+b3U)))
  
  
  IF (n_grids .GE.  2) ALLOCATE(psi_rel2 (b1L:(NN(1,2 )+b1U),b2L:(NN(2,2 )+b2U),b3L:(NN(3,2 )+b3U)))
  IF (n_grids .GE.  3) ALLOCATE(psi_rel3 (b1L:(NN(1,3 )+b1U),b2L:(NN(2,3 )+b2U),b3L:(NN(3,3 )+b3U)))
  IF (n_grids .GE.  4) ALLOCATE(psi_rel4 (b1L:(NN(1,4 )+b1U),b2L:(NN(2,4 )+b2U),b3L:(NN(3,4 )+b3U)))
  IF (n_grids .GE.  5) ALLOCATE(psi_rel5 (b1L:(NN(1,5 )+b1U),b2L:(NN(2,5 )+b2U),b3L:(NN(3,5 )+b3U)))
  IF (n_grids .GE.  6) ALLOCATE(psi_rel6 (b1L:(NN(1,6 )+b1U),b2L:(NN(2,6 )+b2U),b3L:(NN(3,6 )+b3U)))
  IF (n_grids .GE.  7) ALLOCATE(psi_rel7 (b1L:(NN(1,7 )+b1U),b2L:(NN(2,7 )+b2U),b3L:(NN(3,7 )+b3U)))
  IF (n_grids .GE.  8) ALLOCATE(psi_rel8 (b1L:(NN(1,8 )+b1U),b2L:(NN(2,8 )+b2U),b3L:(NN(3,8 )+b3U)))
  IF (n_grids .GE.  9) ALLOCATE(psi_rel9 (b1L:(NN(1,9 )+b1U),b2L:(NN(2,9 )+b2U),b3L:(NN(3,9 )+b3U)))
  IF (n_grids .GE. 10) ALLOCATE(psi_rel10(b1L:(NN(1,10)+b1U),b2L:(NN(2,10)+b2U),b3L:(NN(3,10)+b3U)))
  IF (n_grids .GE. 11) ALLOCATE(psi_rel11(b1L:(NN(1,11)+b1U),b2L:(NN(2,11)+b2U),b3L:(NN(3,11)+b3U)))
  IF (n_grids .GE. 12) ALLOCATE(psi_rel12(b1L:(NN(1,12)+b1U),b2L:(NN(2,12)+b2U),b3L:(NN(3,12)+b3U)))
  IF (n_grids .GE. 13) ALLOCATE(psi_rel13(b1L:(NN(1,13)+b1U),b2L:(NN(2,13)+b2U),b3L:(NN(3,13)+b3U)))
  IF (n_grids .GE. 14) ALLOCATE(psi_rel14(b1L:(NN(1,14)+b1U),b2L:(NN(2,14)+b2U),b3L:(NN(3,14)+b3U)))
  IF (n_grids .GE. 15) ALLOCATE(psi_rel15(b1L:(NN(1,15)+b1U),b2L:(NN(2,15)+b2U),b3L:(NN(3,15)+b3U)))
  !===========================================================================================================
  
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Konvergenzstatistik ===================================================================================
  !===========================================================================================================
  countO = 0
  countP = 0
  countH = 0
  
  ratioO = 0.
  ratioH = 0.
  ratioP = 0.
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Abbruchkriterium der äusseren Iteration ===============================================================
  !===========================================================================================================
  ALLOCATE(precOffset(1:RK_steps,1:n_it_outer))
  ALLOCATE(precRatio (1:RK_steps,1:n_it_outer))
  
  DO substep = 1, RK_steps
    precOffset(substep,:) = precOffset0(substep)
    precRatio (substep,:) = precRatio0 (substep)
  END DO
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Restart-Nr. als String für File-Namen =================================================================
  !===========================================================================================================
  CALL num_to_string(3,restart,restart_char)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Feld-I/O ==============================================================================================
  !===========================================================================================================
  write_large = .FALSE.
  write_med   = .FALSE.
  write_small = .FALSE.
  
  IF (stride_large(1) /= 0 .AND. stride_large(2) /= 0 .AND. stride_large(3) /= 0) write_large = .TRUE.
  IF (stride_med  (1) /= 0 .AND. stride_med  (2) /= 0 .AND. stride_med  (3) /= 0) write_med   = .TRUE.
  IF (stride_small(1) /= 0 .AND. stride_small(2) /= 0 .AND. stride_small(3) /= 0) write_small = .TRUE.
  !===========================================================================================================
  
  
  epsU0 = epsU ! TEST!!! Woanders hin?
  
  
  !IF (rank == 0) CALL SYSTEM('echo "IMPACT is run by $USER on $HOST which is a $HOSTTYPE" | mail -s "IMPACT usage notification" henniger@ifd.mavt.ethz.ch')
  
  
  END SUBROUTINE init_general
  
  
  
  
  
  
  
  
  
  
  !> Subroutine to initialize the parallel computation environment. It creates a new communicator COMM_CART that has
  !! cartesian topology of dimension 3 with (NB1, NB2, NB3) processes for each respective direction, including
  !! information on whether the domain edges are periodic or not. Then it creates
  !! new communicators to slice and bar subgroups. Finally it initializes error handlers for all cartesian 
  !! subgroups.
  SUBROUTINE init_parallel
  
  IMPLICIT NONE
  
  INTEGER                ::  ijkB(1:3)
  INTEGER                ::  periodic(1:3)
  
  
  IF (BC_1L_global == -1) THEN
     periodic(1) = 1
  ELSE
     periodic(1) = 0
  END IF
  
  IF (BC_2L_global == -1) THEN
     periodic(2) = 1
  ELSE
     periodic(2) = 0
  END IF
  
  IF (BC_3L_global == -1) THEN
     periodic(3) = 1
  ELSE
     periodic(3) = 0
  END IF
  
  
  ! Macht keinen messbaren Unterschied (falls doch irgendwann, dann sollte der CALL auch auf die Grobgitter-Kommunikatoren auch angewandt werden!):
  !CALL MPI_CART_CREATE(MPI_COMM_WORLD,3,(/NB1,NB2,NB3/),periodic,.FALSE.,COMM_CART,merror)
  !CALL MPI_CART_MAP(COMM_CART,3,(/NB1,NB2,NB3/),periodic,rank,merror)
  CALL MPI_CART_CREATE(MPI_COMM_WORLD,3,(/NB1,NB2,NB3/),periodic,.TRUE.,COMM_CART,merror)
  CALL MPI_COMM_RANK  (COMM_CART,rank,merror)
  CALL MPI_CART_COORDS(COMM_CART,rank,3,ijkB,merror)
  
  iB(1,1) = ijkB(1)+1
  iB(2,1) = ijkB(2)+1
  iB(3,1) = ijkB(3)+1
  
  iShift = (iB(1,1)-1)*(N1-1)
  jShift = (iB(2,1)-1)*(N2-1)
  kShift = (iB(3,1)-1)*(N3-1)
  
  CALL MPI_CART_SHIFT(COMM_CART,0,1,rank1L,rank1U,merror)
  CALL MPI_CART_SHIFT(COMM_CART,1,1,rank2L,rank2U,merror)
  CALL MPI_CART_SHIFT(COMM_CART,2,1,rank3L,rank3U,merror)
  
  
  CALL MPI_CART_SUB(COMM_CART,(/0,1,1/),COMM_SLICE1,merror)
  CALL MPI_CART_SUB(COMM_CART,(/1,0,1/),COMM_SLICE2,merror)
  CALL MPI_CART_SUB(COMM_CART,(/1,1,0/),COMM_SLICE3,merror)
  
  CALL MPI_CART_SUB(COMM_CART,(/1,0,0/),COMM_BAR1  ,merror)
  CALL MPI_CART_SUB(COMM_CART,(/0,1,0/),COMM_BAR2  ,merror)
  CALL MPI_CART_SUB(COMM_CART,(/0,0,1/),COMM_BAR3  ,merror)
  
  
  CALL MPI_COMM_RANK(COMM_SLICE1,rank_slice1,merror)
  CALL MPI_COMM_RANK(COMM_SLICE2,rank_slice2,merror)
  CALL MPI_COMM_RANK(COMM_SLICE3,rank_slice3,merror)
  
  CALL MPI_COMM_RANK(COMM_BAR1  ,rank_bar1  ,merror)
  CALL MPI_COMM_RANK(COMM_BAR2  ,rank_bar2  ,merror)
  CALL MPI_COMM_RANK(COMM_BAR3  ,rank_bar3  ,merror)
  
  
  ! Evtl. überflüssig ...
  CALL MPI_ERRHANDLER_SET(COMM_CART  ,MPI_ERRORS_ARE_FATAL,merror)
  CALL MPI_ERRHANDLER_SET(COMM_SLICE1,MPI_ERRORS_ARE_FATAL,merror)
  CALL MPI_ERRHANDLER_SET(COMM_SLICE2,MPI_ERRORS_ARE_FATAL,merror)
  CALL MPI_ERRHANDLER_SET(COMM_SLICE3,MPI_ERRORS_ARE_FATAL,merror)
  CALL MPI_ERRHANDLER_SET(COMM_BAR1  ,MPI_ERRORS_ARE_FATAL,merror)
  CALL MPI_ERRHANDLER_SET(COMM_BAR2  ,MPI_ERRORS_ARE_FATAL,merror)
  CALL MPI_ERRHANDLER_SET(COMM_BAR3  ,MPI_ERRORS_ARE_FATAL,merror)
  
  
  END SUBROUTINE init_parallel
  
  
  
  
  
  
  
  
  
  
  !> Subroutine for setting the BCs locally in each processor block. It catches 
  !! scenarios where there is only one process in a certain direction with periodic BCs.
  SUBROUTINE init_boundaries
  
  IMPLICIT NONE
  
  
  ! Default: Nachbarblöcke
  BC_1L = 0
  BC_1U = 0
  
  BC_2L = 0
  BC_2U = 0
  
  BC_3L = 0
  BC_3U = 0
  
  ! Spezialfall periodische RB bei einem Block:
  IF (BC_1L_global == -1 .AND. NB1 == 1) THEN
     BC_1L = -1
     BC_1U = -1
  END IF
  IF (BC_2L_global == -1 .AND. NB2 == 1) THEN
     BC_2L = -1
     BC_2U = -1
  END IF
  IF (BC_3L_global == -1 .AND. NB3 == 1) THEN
     BC_3L = -1
     BC_3U = -1
  END IF
  
  ! Wände oder Symmetrie-RB:
  IF (rank1L .LT. 0) BC_1L = BC_1L_global
  IF (rank1U .LT. 0) BC_1U = BC_1U_global
  
  IF (rank2L .LT. 0) BC_2L = BC_2L_global
  IF (rank2U .LT. 0) BC_2U = BC_2U_global
  
  IF (rank3L .LT. 0) BC_3L = BC_3L_global
  IF (rank3U .LT. 0) BC_3U = BC_3U_global
  
  
  END SUBROUTINE init_boundaries
  
  
  
  
  
  
  
  
  
  
  !> subroutine that initializes limit indices globally and locally inside each block. It also sets up MPI
  !! communicators for multigrid.
  SUBROUTINE init_limits
  
  IMPLICIT NONE
  
  INTEGER                ::  bar1_out(1:NB1)
  INTEGER                ::  bar2_out(1:NB2)
  INTEGER                ::  bar3_out(1:NB3)
  
  INTEGER                ::  g
  
  INTEGER                ::  rank_bar1
  INTEGER                ::  rank_bar2
  INTEGER                ::  rank_bar3
  
  !*************************************************************************
  INTEGER                ::  base_group, new_group, comm_temp
  INTEGER                ::  stride(1:3), counter
  INTEGER                ::  rank_comm1(1:n_grids_max), rank_comm2(1:n_grids_max)
  INTEGER, ALLOCATABLE   ::  new_ranks(:,:,:)
  INTEGER                ::  i, ii, iii, iimax, iiShift
  INTEGER                ::  j, jj, jjj, jjmax, jjShift
  INTEGER                ::  k, kk, kkk, kkmax, kkShift
  INTEGER                ::  offs_global(1:3,NB1*NB2*NB3)
  INTEGER                ::  sizs_global(1:3,NB1*NB2*NB3)
  INTEGER                ::  count_comm
  !*************************************************************************
  
  INTEGER                ::  periodic(1:3)
  LOGICAL                ::  member_yes
  
  IF (BC_1L_global == -1) THEN
     periodic(1) = 1
  ELSE
     periodic(1) = 0
  END IF
  
  IF (BC_2L_global == -1) THEN
     periodic(2) = 1
  ELSE
     periodic(2) = 0
  END IF
  
  IF (BC_3L_global == -1) THEN
     periodic(3) = 1
  ELSE
     periodic(3) = 0
  END IF
  
  
  
  !-----------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - S1p = S1R, N1p = N1R (in Relaxationsroutinen)                                                          !
  !              - Default sind periodische / Nachbarblock-RB                                                             !
  !-----------------------------------------------------------------------------------------------------------------------!
  
  !-----------------------------------------------------------------------------------------------------------------------!
  ! ls =  0: erster  Block ab Wand hat 2 Schichten wand-normale Geschwindigkeiten mehr als die folgenden Blöcke           !
  !                         -||-       1 Schicht   wand-parallele Geschw. / Druck mehr als die folgenden Blöcke           !
  !                                                                                                                       !
  ! ls = -1: erster / letzter Block ab Wand haben 1 Schicht wand-normale Geschwindigkeiten mehr als die Blöcke dazwischen !
  !                   letzter Block ab Wand hat   1 Schicht wand-parallele Geschw. / Druck mehr als die Blöcke davor      !
  ! bbecsek 021214: Anmerkung: 												  !
  ! ls wird gesetzt in mod_vars (ist also hardcoded!)									  !
  !-----------------------------------------------------------------------------------------------------------------------!
  
  
  IF (ls1 ==  0) ex1 =  1
  IF (ls1 == -1) ex1 = -1
  
  IF (ls2 ==  0) ex2 =  1
  IF (ls2 == -1) ex2 = -1
  
  IF (ls3 ==  0) ex3 =  1
  IF (ls3 == -1) ex3 = -1
  
  
  S1R  = 2 + ls1
  S2R  = 2 + ls2
  S3R  = 2 + ls3
  
  S11R = 2 + ls1
  S22R = 2 + ls2
  S33R = 2 + ls3
  
  d1R  = ls1
  d2R  = ls2
  d3R  = ls3
  
  d11R = ls1
  d22R = ls2
  d33R = ls3
  
  
  ! Druck / Konzentrationen (MIT Rand):
  S1p  = 2 + ls1
  S2p  = 2 + ls2
  S3p  = 2 + ls3
  
  N1p  = N1 + ls1
  N2p  = N2 + ls2
  N3p  = N3 + ls3
  
  
  ! Geschwindigkeiten (Feld MIT Rand)
  S11B = 2 + ls1
  S21B = 2 + ls2
  S31B = 2 + ls3
  
  S12B = 2 + ls1
  S22B = 2 + ls2
  S32B = 2 + ls3
  
  S13B = 2 + ls1
  S23B = 2 + ls2
  S33B = 2 + ls3
  
  N11B = N1 + ls1
  N21B = N2 + ls2
  N31B = N3 + ls3
  
  N12B = N1 + ls1
  N22B = N2 + ls2
  N32B = N3 + ls3
  
  N13B = N1 + ls1
  N23B = N2 + ls2
  N33B = N3 + ls3
  
  ! Geschwindigkeiten (Feld OHNE Rand)
  S11 = 2 + ls1
  S21 = 2 + ls2
  S31 = 2 + ls3
  
  S12 = 2 + ls1
  S22 = 2 + ls2
  S32 = 2 + ls3
  
  S13 = 2 + ls1
  S23 = 2 + ls2
  S33 = 2 + ls3
  
  N11 = N1 + ls1
  N21 = N2 + ls2
  N31 = N3 + ls3
  
  N12 = N1 + ls1
  N22 = N2 + ls2
  N32 = N3 + ls3
  
  N13 = N1 + ls1
  N23 = N2 + ls2
  N33 = N3 + ls3
  
  !===========================================================================================================
  IF (BC_1L .GT. 0) THEN
     S1p  = 1
     S1R  = 1
     S11R = 2
     
     S11B = 0
     S12B = 1
     S13B = 1
     
     S11  = 1
     S12  = 2
     S13  = 2
  END IF
  IF (BC_1L == -2) THEN
     S1p  = 1
     S1R  = 1
     S11R = 1
     
     S11B = 1
     S12B = 1
     S13B = 1
     
     S11  = 1
     S12  = 1
     S13  = 1
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_1U .GT. 0) THEN
     N1p  = N1
     d1R  =  0
     d11R = -1
     
     N11B = N1
     N12B = N1
     N13B = N1
     
     N11  = N1-1
     N12  = N1-1
     N13  = N1-1
  END IF
  IF (BC_1U == -2) THEN
     N1p  = N1
     d1R  =  0
     d11R =  0
     
     N11B = N1-1
     N12B = N1
     N13B = N1
     
     N11  = N1-1
     N12  = N1
     N13  = N1
  END IF
  !===========================================================================================================
  IF (BC_2L .GT. 0) THEN
     S2p  = 1
     S2R  = 1
     S22R = 2
     
     S21B = 1
     S22B = 0
     S23B = 1
     
     S21  = 2
     S22  = 1
     S23  = 2
  END IF
  IF (BC_2L == -2) THEN
     S2p  = 1
     S2R  = 1
     S22R = 1
     
     S21B = 1
     S22B = 1
     S23B = 1
     
     S21  = 1
     S22  = 1
     S23  = 1
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_2U .GT. 0) THEN
     N2p  = N2
     d2R  =  0
     d22R = -1
     
     N21B = N2
     N22B = N2
     N23B = N2
     
     N21  = N2-1
     N22  = N2-1
     N23  = N2-1
  END IF
  IF (BC_2U == -2) THEN
     N2p  = N2
     d2R  =  0
     d22R =  0
     
     N21B = N2
     N22B = N2-1
     N23B = N2
     
     N21  = N2
     N22  = N2-1
     N23  = N2
  END IF
  !===========================================================================================================
  IF (BC_3L .GT. 0) THEN
     S3p  = 1
     S3R  = 1
     S33R = 2
     
     S31B = 1
     S32B = 1
     S33B = 0
     
     S31  = 2
     S32  = 2
     S33  = 1
  END IF
  IF (BC_3L == -2) THEN
     S3p  = 1
     S3R  = 1
     S33R = 1
     
     S31B = 1
     S32B = 1
     S33B = 1
     
     S31  = 1
     S32  = 1
     S33  = 1
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_3U .GT. 0) THEN
     N3p  = N3
     d3R  =  0
     d33R = -1
     
     N31B = N3
     N32B = N3
     N33B = N3
     
     N31  = N3-1
     N32  = N3-1
     N33  = N3-1
  END IF
  IF (BC_3U == -2) THEN
     N3p  = N3
     d3R  =  0
     d33R =  0
     
     N31B = N3
     N32B = N3
     N33B = N3-1
     
     N31  = N3
     N32  = N3
     N33  = N3-1
  END IF
  
  
  
  !*********************************************************************
  participate_yes = .FALSE.
  
  stride     = 1
  count_comm = 0
  
  DO g = 1, n_grids
     
     !--------------------------------------------------------------
     stride(1) = stride(1)*n_gather(1,g)
     stride(2) = stride(2)*n_gather(2,g)
     stride(3) = stride(3)*n_gather(3,g)
     !--------------------------------------------------------------
     
     
     !--------------------------------------------------------------
     iB(1,g) = (iB(1,1)-1)*NB(1,g)/NB1 + 1
     iB(2,g) = (iB(2,1)-1)*NB(2,g)/NB2 + 1
     iB(3,g) = (iB(3,1)-1)*NB(3,g)/NB3 + 1
     !--------------------------------------------------------------
     
     
     !--------------------------------------------------------------
     ! Analog "init_boundaries"
     BC(1:2,1,g) = (/BC_1L,BC_1U/)
     BC(1:2,2,g) = (/BC_2L,BC_2U/)
     BC(1:2,3,g) = (/BC_3L,BC_3U/)
     
     IF (iB(1,g) == 1       .AND. .NOT. BC_1L_global == -1) BC(1,1,g) = BC_1L_global
     IF (iB(1,g) == NB(1,g) .AND. .NOT. BC_1U_global == -1) BC(2,1,g) = BC_1U_global
     
     IF (iB(2,g) == 1       .AND. .NOT. BC_2L_global == -1) BC(1,2,g) = BC_2L_global
     IF (iB(2,g) == NB(2,g) .AND. .NOT. BC_2U_global == -1) BC(2,2,g) = BC_2U_global
     
     IF (iB(3,g) == 1       .AND. .NOT. BC_3L_global == -1) BC(1,3,g) = BC_3L_global
     IF (iB(3,g) == NB(3,g) .AND. .NOT. BC_3U_global == -1) BC(2,3,g) = BC_3U_global
     
     
     IF (NB(1,g) == 1 .AND. BC_1L_global == -1) BC(1:2,1,g) = BC_1L_global
     IF (NB(2,g) == 1 .AND. BC_2L_global == -1) BC(1:2,2,g) = BC_2L_global
     IF (NB(3,g) == 1 .AND. BC_3L_global == -1) BC(1:2,3,g) = BC_3L_global
     !--------------------------------------------------------------
     
     
     !--------------------------------------------------------------
     SNB(1:2,1,g) = (/2+ls1,NN(1,g)+ls1/)
     SNB(1:2,2,g) = (/2+ls2,NN(2,g)+ls2/)
     SNB(1:2,3,g) = (/2+ls3,NN(3,g)+ls3/)
     
     IF (BC(1,1,g) .GT. 0) SNB(1,1,g) = 1
     IF (BC(1,1,g)  == -2) SNB(1,1,g) = 1
     IF (BC(2,1,g) .GT. 0) SNB(2,1,g) = NN(1,g)
     IF (BC(2,1,g)  == -2) SNB(2,1,g) = NN(1,g)
     
     IF (BC(1,2,g) .GT. 0) SNB(1,2,g) = 1
     IF (BC(1,2,g)  == -2) SNB(1,2,g) = 1
     IF (BC(2,2,g) .GT. 0) SNB(2,2,g) = NN(2,g)
     IF (BC(2,2,g)  == -2) SNB(2,2,g) = NN(2,g)
     
     IF (BC(1,3,g) .GT. 0) SNB(1,3,g) = 1
     IF (BC(1,3,g)  == -2) SNB(1,3,g) = 1
     IF (BC(2,3,g) .GT. 0) SNB(2,3,g) = NN(3,g)
     IF (BC(2,3,g)  == -2) SNB(2,3,g) = NN(3,g)
     !--------------------------------------------------------------
     
     
     !--------------------------------------------------------------
     SNF(1:2,1,g) = (/2+ls1,NN(1,g)+ls1/)
     SNF(1:2,2,g) = (/2+ls2,NN(2,g)+ls2/)
     SNF(1:2,3,g) = (/2+ls3,NN(3,g)+ls3/)
     
     IF (BC(1,1,g) .GT. 0) SNF(1,1,g) = 2
     IF (BC(1,1,g)  == -2) SNF(1,1,g) = 1
     IF (BC(2,1,g) .GT. 0) SNF(2,1,g) = NN(1,g)-1
     IF (BC(2,1,g)  == -2) SNF(2,1,g) = NN(1,g)
     
     IF (BC(1,2,g) .GT. 0) SNF(1,2,g) = 2
     IF (BC(1,2,g)  == -2) SNF(1,2,g) = 1
     IF (BC(2,2,g) .GT. 0) SNF(2,2,g) = NN(2,g)-1
     IF (BC(2,2,g)  == -2) SNF(2,2,g) = NN(2,g)
     
     IF (BC(1,3,g) .GT. 0) SNF(1,3,g) = 2
     IF (BC(1,3,g)  == -2) SNF(1,3,g) = 1
     IF (BC(2,3,g) .GT. 0) SNF(2,3,g) = NN(3,g)-1
     IF (BC(2,3,g)  == -2) SNF(2,3,g) = NN(3,g)
     !--------------------------------------------------------------
     
     
     !--------------------------------------------------------------
     CALL MPI_CART_SHIFT(COMM_CART,0,stride(1),ngb(1,1,g),ngb(2,1,g),merror)
     CALL MPI_CART_SHIFT(COMM_CART,1,stride(2),ngb(1,2,g),ngb(2,2,g),merror)
     CALL MPI_CART_SHIFT(COMM_CART,2,stride(3),ngb(1,3,g),ngb(2,3,g),merror)
     !--------------------------------------------------------------
     
     
     !--------------------------------------------------------------
     !--- coarse grid communicators (comm1) ---
     ! siehe auch Kommentare in alten Versionen bis 19.01.2010!
     IF (g .GE. 2 .AND. n_gather(1,g)*n_gather(2,g)*n_gather(3,g) .GT. 1) THEN
     
        ALLOCATE(new_ranks(1:NB(1,g),1:NB(2,g),1:NB(3,g)))
        
        !-----------------------------------------------------------------------------------
        CALL MPI_COMM_GROUP(COMM_CART,base_group,merror)
        participate_yes(g) = .FALSE.
        DO k = 1, NB(3,g)
           DO j = 1, NB(2,g)
              DO i = 1, NB(1,g)
                 CALL MPI_CART_RANK(COMM_CART,(/MOD((i-1)*stride(1),NB1),  &
                                     &          MOD((j-1)*stride(2),NB2),  &
                                     &          MOD((k-1)*stride(3),NB3)/),new_ranks(i,j,k),merror)
                 IF (rank == new_ranks(i,j,k)) participate_yes(g) = .TRUE.
              END DO
           END DO
        END DO
        
        CALL MPI_GROUP_INCL(base_group,NB(1,g)*NB(2,g)*NB(3,g),new_ranks,new_group,merror)
        CALL MPI_COMM_CREATE(COMM_CART,new_group,comm_temp,merror)
        CALL MPI_GROUP_FREE(base_group,merror)
        CALL MPI_GROUP_FREE(new_group ,merror)
        
        IF (participate_yes(g)) THEN
           CALL MPI_CART_CREATE(comm_temp,3,NB(1:3,g),periodic,.TRUE.,comm1(g),merror)
           CALL MPI_COMM_FREE(comm_temp,merror)
        END IF
        !-----------------------------------------------------------------------------------
        
        DEALLOCATE(new_ranks)
     ELSE
        IF (g == 1) THEN
           comm1          (g) = COMM_CART
           rank_comm1     (g) = rank
           participate_yes(g) = .TRUE.
        ELSE
           comm1          (g) = comm1          (g-1)
           rank_comm1     (g) = rank_comm1     (g-1)
           participate_yes(g) = participate_yes(g-1)
        END IF
     END IF
     !--------------------------------------------------------------
     
     
     !--------------------------------------------------------------
     !--- gathering communicators (comm2) ---
     IF (g .GE. 2 .AND. n_gather(1,g)*n_gather(2,g)*n_gather(3,g) .GT. 1) THEN
        
        ALLOCATE(new_ranks(1:n_gather(1,g),1:n_gather(2,g),1:n_gather(3,g)))
        
        !-----------------------------------------------------------------------------------------------
        CALL MPI_COMM_GROUP(COMM_CART,base_group,merror)
        
        DO kk = 1, NB(3,g)
           DO jj = 1, NB(2,g)
              DO ii = 1, NB(1,g)
                 
                 member_yes = .FALSE.
                 DO k = 1, n_gather(3,g)
                    DO j = 1, n_gather(2,g)
                       DO i = 1, n_gather(1,g)
                          ! TEST!!! Das sieht etwas komisch aus ...
                          CALL MPI_CART_RANK(COMM_CART,(/MOD(((ii-1)*n_gather(1,g)+i-1)*NB1/NB(1,g-1),NB1),                             &
                                            &            MOD(((jj-1)*n_gather(2,g)+j-1)*NB2/NB(2,g-1),NB2),                             &
                                            &            MOD(((kk-1)*n_gather(3,g)+k-1)*NB3/NB(3,g-1),NB3)/),new_ranks(i,j,k),merror)
                          IF (rank == new_ranks(i,j,k)) member_yes = .TRUE.
                       END DO
                    END DO
                 END DO
                 
                 CALL MPI_GROUP_INCL(base_group,n_gather(1,g)*n_gather(2,g)*n_gather(3,g),new_ranks,new_group,merror)
                 CALL MPI_COMM_CREATE(COMM_CART,new_group,comm_temp,merror)
                 CALL MPI_GROUP_FREE(new_group,merror)
                 
                 IF (member_yes) THEN
                    CALL MPI_CART_CREATE(comm_temp,3,n_gather(1:3,g),periodic,.TRUE.,comm2(g),merror)
                    CALL MPI_COMM_FREE(comm_temp,merror)
                    
                    !--- comm1 und comm2 synchronisieren ---
                    rank_comm2(g) = 0
                    IF (participate_yes(g)) CALL MPI_COMM_RANK(comm2(g),rank_comm2(g),merror) ! Nur fuer "rankc2(g)" gilt "participate_yes(g)"
                    CALL MPI_ALLREDUCE(rank_comm2(g),rankc2(g),1,MPI_INTEGER,MPI_SUM,comm2(g),merror)
                    CALL MPI_COMM_RANK(comm2(g),rank_comm2(g),merror)
                 END IF
                 
              END DO
           END DO
        END DO
        
        CALL MPI_GROUP_FREE(base_group,merror)
        !-----------------------------------------------------------------------------------------------
        
        DEALLOCATE(new_ranks)
        
     END IF
     !--------------------------------------------------------------
     
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     IF (g .GE. 2 .AND. n_gather(1,g)*n_gather(2,g)*n_gather(3,g) .GT. 1 ) THEN
     IF (participate_yes(g-1)) THEN
        
        iimax   = (NN(1,g)-1)/n_gather(1,g)+1
        jjmax   = (NN(2,g)-1)/n_gather(2,g)+1
        kkmax   = (NN(3,g)-1)/n_gather(3,g)+1
        
        iiShift = (iimax-1)*MOD(iB(1,g-1)-1,n_gather(1,g))
        jjShift = (jjmax-1)*MOD(iB(2,g-1)-1,n_gather(2,g))
        kkShift = (kkmax-1)*MOD(iB(3,g-1)-1,n_gather(3,g))
        
        CALL MPI_COMM_RANK(comm2(g),rank_comm2(g),merror)
        
        !--- Interpolation ---------------------------------------------------------------------
        offs_global = 0
        sizs_global = 0
        
        offs_global(1:3,rank_comm2(g)+1) = (/iiShift,jjShift,kkShift/)
        sizs_global(1:3,rank_comm2(g)+1) = (/iimax  ,jjmax  ,kkmax  /)
        
        CALL MPI_ALLREDUCE(offs_global,offsI(1,1,g),3*n_gather(1,g)*n_gather(2,g)*n_gather(3,g),MPI_INTEGER,MPI_SUM,comm2(g),merror)
        CALL MPI_ALLREDUCE(sizs_global,sizsI(1,1,g),3*n_gather(1,g)*n_gather(2,g)*n_gather(3,g),MPI_INTEGER,MPI_SUM,comm2(g),merror)
        
        counter = 0
        DO k = 1, n_gather(3,g)
           DO j = 1, n_gather(2,g)
              DO i = 1, n_gather(1,g)
                 recvI                (  i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)      &
                    &          = sizsI(1,i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)      &
                    &           *sizsI(2,i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)      &
                    &           *sizsI(3,i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)
                 
                 dispI                    (  i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g) = counter
                 counter = counter + recvI(  i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)
              END DO
           END DO
        END DO
        
        
        !--- Restriction -----------------------------------------------------------------------
        IF (ls1 ==  0 .AND. (BC(1,1,g-1) == 0 .OR. BC(1,1,g-1) == -1)) iimax = iimax-1
        IF (ls1 == -1 .AND. (BC(2,1,g-1) == 0 .OR. BC(2,1,g-1) == -1)) iimax = iimax-1
        
        IF (ls2 ==  0 .AND. (BC(1,2,g-1) == 0 .OR. BC(1,2,g-1) == -1)) jjmax = jjmax-1
        IF (ls2 == -1 .AND. (BC(2,2,g-1) == 0 .OR. BC(2,2,g-1) == -1)) jjmax = jjmax-1
        
        IF (ls3 ==  0 .AND. (BC(1,3,g-1) == 0 .OR. BC(1,3,g-1) == -1)) kkmax = kkmax-1
        IF (ls3 == -1 .AND. (BC(2,3,g-1) == 0 .OR. BC(2,3,g-1) == -1)) kkmax = kkmax-1
        
        IF (ls1 ==  0 .AND. (BC(1,1,g-1) == 0 .OR. BC(1,1,g-1) == -1)) iiShift = iiShift+1
        IF (ls2 ==  0 .AND. (BC(1,2,g-1) == 0 .OR. BC(1,2,g-1) == -1)) jjShift = jjShift+1
        IF (ls3 ==  0 .AND. (BC(1,3,g-1) == 0 .OR. BC(1,3,g-1) == -1)) kkShift = kkShift+1
        
        offs_global = 0
        sizs_global = 0
        
        offs_global(1:3,rank_comm2(g)+1) = (/iiShift,jjShift,kkShift/)
        sizs_global(1:3,rank_comm2(g)+1) = (/iimax  ,jjmax  ,kkmax  /)
        
        CALL MPI_ALLREDUCE(offs_global,offsR(1,1,g),3*n_gather(1,g)*n_gather(2,g)*n_gather(3,g),MPI_INTEGER,MPI_SUM,comm2(g),merror)
        CALL MPI_ALLREDUCE(sizs_global,sizsR(1,1,g),3*n_gather(1,g)*n_gather(2,g)*n_gather(3,g),MPI_INTEGER,MPI_SUM,comm2(g),merror)
        
        counter = 0
        DO k = 1, n_gather(3,g)
           DO j = 1, n_gather(2,g)
              DO i = 1, n_gather(1,g)
                 recvR                (  i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)      &
                    &          = sizsR(1,i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)      &
                    &           *sizsR(2,i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)      &
                    &           *sizsR(3,i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)
                 
                 dispR                    (  i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g) = counter
                 counter = counter + recvR(  i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)
              END DO
           END DO
        END DO
        !---------------------------------------------------------------------------------------
     END IF
     END IF
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     
  END DO
  !*********************************************************************
  
  !===========================================================================================================
  
  dim1 = M1
  dim2 = M2
  dim3 = M3
  
  IF (BC_1L_global == -1) dim1 = M1-1
  IF (BC_2L_global == -1) dim2 = M2-1
  IF (BC_3L_global == -1) dim3 = M3-1
  
  !===========================================================================================================
  
  bar1_size = 0
  bar2_size = 0
  bar3_size = 0

  bar1_offset = 0
  bar2_offset = 0
  bar3_offset = 0
  
  CALL MPI_COMM_RANK(COMM_BAR1,rank_bar1,merror)
  CALL MPI_COMM_RANK(COMM_BAR2,rank_bar2,merror)
  CALL MPI_COMM_RANK(COMM_BAR3,rank_bar3,merror)
  
  bar1_size(rank_bar1+1) = N1p-S1p+1
  bar2_size(rank_bar2+1) = N2p-S2p+1
  bar3_size(rank_bar3+1) = N3p-S3p+1
  
  CALL MPI_ALLREDUCE(bar1_size,bar1_out,NB1,MPI_INTEGER,MPI_SUM,COMM_BAR1,merror)
  CALL MPI_ALLREDUCE(bar2_size,bar2_out,NB2,MPI_INTEGER,MPI_SUM,COMM_BAR2,merror)
  CALL MPI_ALLREDUCE(bar3_size,bar3_out,NB3,MPI_INTEGER,MPI_SUM,COMM_BAR3,merror)
  
  bar1_size = bar1_out
  bar2_size = bar2_out
  bar3_size = bar3_out
  
  
  bar1_offset(rank_bar1+1) = iShift + 1 + ls1
  bar2_offset(rank_bar2+1) = jShift + 1 + ls2
  bar3_offset(rank_bar3+1) = kShift + 1 + ls3
  
  CALL MPI_ALLREDUCE(bar1_offset,bar1_out,NB1,MPI_INTEGER,MPI_SUM,COMM_BAR1,merror)
  CALL MPI_ALLREDUCE(bar2_offset,bar2_out,NB2,MPI_INTEGER,MPI_SUM,COMM_BAR2,merror)
  CALL MPI_ALLREDUCE(bar3_offset,bar3_out,NB3,MPI_INTEGER,MPI_SUM,COMM_BAR3,merror)
  
  bar1_offset = bar1_out
  bar2_offset = bar2_out
  bar3_offset = bar3_out
  
  
  END SUBROUTINE init_limits
  
  
  
END MODULE mod_setup
