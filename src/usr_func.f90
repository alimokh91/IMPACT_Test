!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!*************************************************************************************************************

!> module containing functions specified by the user
MODULE usr_func
  
  
  USE mod_dims
  USE mod_vars
  USE mod_lib
  USE mod_diff
  USE mod_exchange
  USE usr_vars
  USE ISO_C_BINDING !bbecsek
  USE mod_solvers !kschlegel
  USE mod_inout !kschlegel
  USE HDF5 !kschlegel
  USE mr_io_protocol
  USE mr_io_parallel_spacetime
  USE MPI  
  
  PRIVATE
  
  PUBLIC smooth_step !bbecsek
  PUBLIC linear_step !bbecsek
  PUBLIC fringe_coeff !bbecsek
  PUBLIC poiseuille_parabola !bbecsek
  PUBLIC interface, erf, atanh
  PUBLIC coord_tanh
  PUBLIC coord_gausslobatto !bbecsek
  PUBLIC init_particles
  PUBLIC init_hdf5, init_mpi, finl_hdf5, finl_mpi, pass_mpi_comm, print_fcomm_size, mpi_bcast_fort !bbecsek
  PUBLIC windkessel_integration_step, pdot3EWK, pdot4EWK
  PUBLIC flow_and_rate, flow_and_rate2D, flow_and_rate3D, fitted_flow 
  PUBLIC apply_fringe_forcing, apply_windkessel_loading
  PUBLIC load_vel_initcond_hdf5
  PUBLIC check_node_ids
  PUBLIC block_id, block_cart
  PUBLIC total_n_local_tet_elements, total_n_global_tet_elements
  PUBLIC total_n_local_hex_elements, total_n_global_hex_elements
  PUBLIC total_n_local_tri_elements, total_n_global_tri_elements
  PUBLIC total_n_local_quad_elements, total_n_global_quad_elements
  PUBLIC local_to_global_node_id_loc_con_allperiodic2
  PUBLIC global_to_local_node_id_loc_con_allperiodic2
  PUBLIC local_to_global_node_id_loc_con ! new general routines
  PUBLIC global_to_local_node_id_loc_con ! new general routines
  PUBLIC global_id2cart_loc_con, global_cart2id_loc_con
  PUBLIC n_local_nodes, n_global_nodes
  PUBLIC local_id2cart_loc_con, local_cart2id_loc_con
  PUBLIC local_tet_element_nodes
  PUBLIC local_hex_element_nodes
  PUBLIC local_quad_element_nodes
  PUBLIC local_tri_element_nodes
  PUBLIC interpolate_force_pre_vel
  PUBLIC residual2volume_force
  PUBLIC local_to_global_tet_elem_id
  PUBLIC local_to_global_hex_elem_id
  PUBLIC local_to_global_quad_elem_id
  PUBLIC local_to_global_tri_elem_id
  !PUBLIC set_pointers_to_non_allocatable_arrays
  PUBLIC rhs_L2_norm, vel_l2_norm, force_l2_norm
  PUBLIC compute_and_store_global_block_sizes
  PUBLIC reset_force_density_component
  PUBLIC open_log_iterations_fortran, close_log_iterations_fortran
  PUBLIC save_old_velocity, restore_old_velocity

  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  !=========================================================================================================!
  !--- user-specific subroutines and functions ---                                                          !
  !=========================================================================================================!
  
  
  FUNCTION smooth_step(step_yes , xx) RESULT(fn_val)

  IMPLICIT NONE

  LOGICAL                :: step_yes
  REAL   , INTENT(IN)    :: xx
  REAL                   :: fn_val
  
  IF (step_yes) THEN
     IF       (xx .LE. 0.) THEN
       fn_val = 0.
     ELSE IF  ( (xx .GT. 0.) .AND. (xx .LT. 1.) ) THEN
       fn_val = 1./( 1. + EXP( (1./(xx-1.)) + (1./xx)) )
     ELSE IF  (xx .GE. 1.) THEN
       fn_val = 1.
     ELSE
       write(*,*)'     WARNING: invalid input to function smooth_step.'
     END IF
  ELSE
    fn_val = 1.
  END IF

  RETURN

  END FUNCTION smooth_step


  
  FUNCTION linear_step(xx) RESULT(fn_val)

  IMPLICIT NONE

  REAL   , INTENT(IN)     :: xx
  REAL                    :: fn_val

  IF ( xx .LE. 0 ) THEN
    fn_val = 0.
  ELSEIF ( (xx .GT. 0. ) .AND. (xx .LT. 1.) ) THEN
    fn_val = xx
  ELSEIF (xx .GE. 1.) THEN
    fn_val = 1.
  ELSE
    write(*,*)'     WARNING: invalid input to function linear_step.'
  END IF

  END FUNCTION linear_step



  !> function that computes the fringe_coeff lambda(x)
  !! @note: 
  SUBROUTINE fringe_coeff(lam_max, x_start, x_end, d_rise, d_fall, xx, lam_fringe)

  IMPLICIT NONE

  REAL   , INTENT(IN)    :: lam_max
  REAL   , INTENT(IN)    :: x_start, x_end
  REAL   , INTENT(IN)    :: d_rise, d_fall
  REAL   , INTENT(IN)    :: xx

  REAL   , INTENT(OUT)   :: lam_fringe

  IF ((d_rise .NE. d_fall) .AND. (x_start .LT. L1) .AND. (x_start .LT. L2))  THEN
    !write(*,*)'      WARNING: make sure rise and fall intervals equal.'
  END IF

  lam_fringe = lam_max*( smooth_step(.TRUE.,(xx-x_start)/d_rise) - smooth_step(.TRUE.,1. + ( (xx-x_end)/d_fall )) )

  IF ( (x_end .LE. (x_start + d_rise)) .OR. (x_start .GE. (x_end - d_fall)) ) THEN
    lam_fringe = lam_fringe + lam_max
  END IF

  END SUBROUTINE fringe_coeff





  SUBROUTINE poiseuille_parabola(x_start , x_end , xx, parab_val)
 
  IMPLICIT NONE

  REAL   , INTENT(IN)    :: x_start
  REAL   , INTENT(IN)    :: x_end
  REAL   , INTENT(IN)    :: xx
  REAL                   :: a, b, c

  REAL   , INTENT(OUT)   :: parab_val


  IF ((xx .GE. x_start) .AND. (xx .LE. x_end)) THEN
    a = -4.                 /((x_end - x_start)**2.)
    b = (4.*(x_end+x_start))/((x_end - x_start)**2.)
    c = -(4.*x_end*x_start) /((x_end - x_start)**2.)

    parab_val = a*(xx**2.) + (b*xx) + c
  ELSE
    parab_val = 0.
    
  END IF


  END SUBROUTINE poiseuille_parabola



  SUBROUTINE poiseuille_paraboloid(xx, yy, x_center, y_center, radius, parab_val)

  IMPLICIT NONE

  REAL   , INTENT(IN)    :: xx, yy
  REAL   , INTENT(IN)    :: x_center, y_center
  REAL   , INTENT(IN)    :: radius

  REAL   , INTENT(OUT)   :: parab_val


  parab_val = 0.

  IF ( ((xx - x_center)**2. + (yy - y_center)**2.) .LE. radius**2. ) THEN

    parab_val = -1. * ( ( (xx - x_center)**2. + (yy - y_center)**2. ) &
                         /(radius**2.) - 1. )

  END IF


  END SUBROUTINE poiseuille_paraboloid



  SUBROUTINE windkessel_integration_step

  IMPLICIT NONE

  REAL                  :: WK_predot_old(1:3)

  WK_predot_old = WK_predot

  IF (WK_type .EQ. 3) THEN
    CALL pdot3EWK( WK_pre(2) , WK_predot(2) )
  ELSE IF (WK_type .EQ. 4) THEN
    CALL pdot4EWK( WK_pre , WK_predot )
  ELSE
    WRITE(*,*)'WARNING: Invalid Windkessel model specified. Aborting ...'
    CALL MPI_FINALIZE(merror)
    STOP
  END IF

  WK_pre = WK_pre + dtime*aRK(substep)*WK_predot + dtime*bRK(substep)*WK_predot_old
  
  !--- Algebraic equation for 3EWK, Q_AV is global and already computed ---
  IF (WK_type .EQ. 3) WK_pre(1) = WK_pre(2) + R_c*Q_AV

!*** debugging
!print*, subtime, WK_pre(1), WK_pre(2), WK_pre(3), Q_AV, dQ_AV_dt
  END SUBROUTINE windkessel_integration_step






  !> subroutine that returns the RHS of the 3EWK model
  SUBROUTINE pdot3EWK( pb , dpb_dt )

  IMPLICIT NONE

  REAL   , INTENT(IN)    :: pb
  REAL   , INTENT(OUT)   :: dpb_dt

  Q_AV = flow_and_rate( WK_flow_dir , WK_flow_pos , WK_flow_center , WK_flow_radius )
  Q_AV = Q_AV*U_ref*L_ref**2. ! make dimensional
  !Q_AV = fitted_flow(MOD(subtime,1.))

  dpb_dt = -1./(R_p * C_art) * pb + Q_AV/C_art

  END SUBROUTINE pdot3EWK




  !> subroutine that returns the RHS of the 4EWK model
  SUBROUTINE pdot4EWK( p , dp_dt )

  IMPLICIT NONE

  REAL   , INTENT(IN)    :: p(1:3)
  REAL   , INTENT(OUT)   :: dp_dt(1:3)
  REAL                   :: Q_AV_old

  !Q_AV_old = Q_AV
  !Q_AV     = fitted_flow(MOD(subtime,1.))

!=== This computes the flow rate by integrating the RHS ===
  dQ_AV_dt = flow_and_rate( WK_flow_dir , WK_flow_pos , WK_flow_center , WK_flow_radius )
  dQ_AV_dt = dQ_AV_dt*L_ref*U_ref**2. ! make dimensional  
!=== This computes the flow rate from the flow by FD differentiation ===
!--- forward Euler ---
!  dQ_AV_dt = (Q_AV - Q_AV_old) / (dtime*(aRK(substep)+bRK(substep)))

!--- inverted RK3 scheme ---
!  IF (substep == 1) THEN
!    Q_AV_0   = Q_AV
!    !Q_AV    = flow( WK_flow_dir , WK_flow_pos , WK_flow_center , WK_flow_radius )
!    Q_AV_1   = fitted_flow(MOD(subtime,1.))
!    dQ_AV_dt = (Q_AV_1 - Q_AV_0)/(dtime*aRK(1))
!  ELSE IF (substep == 2) THEN
!    Q_AV_2   = fitted_flow(MOD(subtime,1.))
!    dQ_AV_dt = (-(aRK(1)+bRK(2))*(Q_AV_1-Q_AV_0)/(aRK(1)*aRK(2)) + (Q_AV_2 - Q_AV_0)/aRK(2) )/dtime
!  ELSE IF (substep == 3) THEN
!    Q_AV_3   = fitted_flow(MOD(subtime,1.))
!    dQ_AV_dt = ( bRK(3)*(aRK(1) + bRK(2))*(Q_AV_1 - Q_AV_0)/(aRK(1)*aRK(2)*aRK(3)) &
!                       -(aRK(2) + bRK(3))*(Q_AV_2 - Q_AV_0)/(aRK(2)*aRK(3)) &
!                       +                  (Q_AV_3 - Q_AV_0)/aRK(3))/dtime
!    Q_AV   = Q_AV_3
!  END IF
!========================================================================

  dp_dt(1) = -R_c/L_art*p(1) + R_c/L_art*p(2) +               p(3) + R_c * dQ_AV_dt 
  dp_dt(2) =                                                      p(3)
  dp_dt(3) =                                      -1./(R_p*C_art)*p(3) + dQ_AV_dt / C_art

  END SUBROUTINE pdot4EWK



  !> function returns the correct flow or flow rate depending on the WK type at a desired location
  !! WK4 -> flow rate; WK3 -> flow
  FUNCTION flow_and_rate(dir , pos , center , radius) RESULT(fn_val)

  IMPLICIT NONE
 
  INTEGER, INTENT(IN)    :: dir
  REAL   , INTENT(IN)    :: pos
  REAL   , INTENT(IN)    :: center(1:3) !< center(3)==1 in 2D
  REAL   , INTENT(IN)    :: radius
  REAL                   :: fn_val
  REAL                   :: fn_val_global

  IF (dimens .EQ. 2) THEN
    fn_val = flow_and_rate2D(dir, pos, center, radius)
  ELSE IF (dimens .EQ. 3) THEN
    fn_val = flow_and_rate3D(dir, pos, center, radius)
  END IF

  !--- Sum up the flow contribution of all MPI processes
  fn_val_global = 0.
  CALL MPI_ALLREDUCE(fn_val, fn_val_global, 1, MPI_REAL8, MPI_SUM, COMM_CART, merror)
  fn_val = fn_val_global

  END FUNCTION flow_and_rate


  !> function returns the 2D flow along an axis
  FUNCTION flow_and_rate2D(dir, pos , center, radius) RESULT(fn_val)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: dir
  REAL   , INTENT(IN)    :: pos
  REAL   , INTENT(IN)    :: center(1:3) !< center(3)==1 in 2D
  REAL   , INTENT(IN)    :: radius
  REAL                   :: x1start, x1end, x2start, x2end
  INTEGER                :: x1loc, x2loc
  REAL                   :: fn_val
  INTEGER                :: i,j,k


  fn_val = 0.
 
  !=== in x-direction =======================================================================================
  IF (dir .EQ. 1) THEN
    x2start = center(2) - radius
    x2end   = center(2) + radius
    !--- 3EWK -> return flow --------------------------------------------------------------------------------
    IF (WK_type .EQ. 3) THEN
      IF ((pos .GE. x1p(S1p)) .AND. (pos .LE. x1p(N1p))) THEN ! check if desired position is in process
        x1loc   = MINLOC(ABS(x1p - pos), 1) !< finds the index closest to the specified location

        DO k = S3p, N3p ! this should be 1 anyway in 2D
          DO j = S2p, N2p
            DO i = x1loc, x1loc
              IF ( (x2p(j) .GE. x2start) .AND. (x2p(j) .LE. x2end) ) THEN 
                fn_val     = fn_val + work1(i,j,k)*dx1p(i)*dx2p(j)*dx3p(k)
              END IF
            END DO 
          END DO
        END DO
      ELSE
       fn_val = 0.
      END IF
    !--- 4EWK -> return flow rate ---------------------------------------------------------------------------
    ELSE IF (WK_type .EQ. 4) THEN
      IF ((pos .GE. x1u(S11B)) .AND. (pos .LE. x1u(N11B))) THEN ! check if desired position is in process
        x1loc   = MINLOC(ABS(x1u - pos), 1) !< finds the index closest to the specified location
     
        DO k = S31B, N31B
          DO j = S21B, N21B
            DO i = x1loc, x1loc
              IF ( (x2p(j) .GE. x2start) .AND. (x2p(j) .LE. x2end) ) THEN
                fn_val     = fn_val + rhs(i,j,k,1)*dx1u(i)*dx1p(j)*dx3p(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    END IF
  !=== in y-direction ======================================================================================= 
  ELSE IF (dir .EQ. 2) THEN
    x1start = center(1) - radius
    x1end   = center(1) + radius
    !--- 3EWK -> return flow --------------------------------------------------------------------------------
    IF (WK_type .EQ. 3) THEN
      IF ((pos .GE. x2p(S2p)) .AND. (pos .LE. x2p(N2p))) THEN ! check if desired position is in process
        x2loc   = MINLOC(ABS(x2p - pos), 1) !< finds the index closest to the specified location 

        DO k = S3p, N3p ! this should be 1 anyway in 2D
          DO j = x2loc, x2loc
            DO i = S1p, N1p
              IF ( (x1p(i) .GE. x1start) .AND. (x1p(i) .LE. x1end) ) THEN
                fn_val     = fn_val + work2(i,j,k)*dx1p(i)*dx2p(j)*dx3p(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    !--- 4EWK -> return flow rate ---------------------------------------------------------------------------
    ELSE IF (WK_type .EQ. 4) THEN
      IF ((pos .GE. x2v(S22B)) .AND. (pos .LE. x2v(N22B))) THEN ! check if desired position is in process
        x2loc   = MINLOC(ABS(x2v - pos), 1) !< finds the index closest to the specified location 

        DO k = S32B, N32B ! this should be 1 anyway in 2D
          DO j = x2loc, x2loc
            DO i = S12B, N12B
              IF ( (x1p(i) .GE. x1start) .AND. (x1p(i) .LE. x1end) ) THEN
                fn_val     = fn_val + rhs(i,j,k,2)*dx1p(i)*dx2v(j)*dx3p(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    END IF
  !=== in z-direction: NOT POSSIBLE =========================================================================
  ELSE IF (dir .EQ. 3 .AND. dimens .NE. 3) THEN
     WRITE(*,*)'WARNING: cannot compute the flux in 3-direction for a 2D flow. Aborting...'
     CALL MPI_FINALIZE(merror)
     STOP
  !==========================================================================================================
  ELSE
     WRITE(*,*)'WARNING: something went wrong with calculating the flow. Aborting...'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF

  END FUNCTION flow_and_rate2D



  !> returns the flow in 3D through a circular opening in the direction of an axis
  FUNCTION flow_and_rate3D(dir , pos , center , radius) RESULT(fn_val)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: dir
  REAL   , INTENT(IN)    :: pos
  REAL   , INTENT(IN)    :: center(1:3) !< center(3)==1 in 2D
  REAL   , INTENT(IN)    :: radius
  REAL                   :: x1start, x1end, x2start, x2end, x3start, x3end
  INTEGER                :: x1loc, x2loc, x3loc
  REAL                   :: fn_val
  INTEGER                :: i,j,k

  !==========================================================================================================
  ! Depending on wether we use 4EWK or 3EWK we need either the flow or the flow rate to be returned.
  ! We make use of the fact that the integral over du/dt = rhs equals the flow rate.
  ! Because the components u,v,w are stored in containers work1,2,3 on the pressure grid, we can use the 
  ! pressure indices to loop. However, the RHS is only given on the velocity grids, i.e. we need to loop 
  ! over different indices. It may at a certain point be nicer to reduce the code and include the WK type 
  ! checking into the triple loops and use vel(:,:,:,i) insted of the worki(:,:,:)
  !==========================================================================================================


  fn_val = 0.

  !=== in x-direction =======================================================================================
  IF (dir .EQ. 1) THEN
    !--- 3EWK -> return flow --------------------------------------------------------------------------------
    IF (WK_type .EQ. 3) THEN
      IF ((pos .GE. x1p(S1p)) .AND. (pos .LE. x1p(N1p))) THEN ! check if desired position is in process
        x1loc = MINLOC(ABS(x1p - pos), 1) !< finds the index closest to the specified location

        DO k = S3p, N3p
          DO j = S2p, N2p
            DO i = x1loc, x1loc
              IF ( (x2p(j) - center(2))**2 + (x3p(k) - center(3))**2 .LE. radius**2) THEN
                fn_val     = fn_val + work1(i,j,k)*dx1p(i)*dx2p(j)*dx3p(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    !--- 4EWK -> return flow rate ---------------------------------------------------------------------------
    ELSE IF (WK_type .EQ. 4) THEN
      IF ((pos .GE. x1u(S11B)) .AND. (pos .LE. x1u(N11B))) THEN ! check if desired position is in process
        x1loc = MINLOC(ABS(x1u - pos), 1) !< finds the index closest to the specified location

        DO k = S31B, N31B
          DO j = S21B, N21B
            DO i = x1loc, x1loc
              IF ( (x2p(j) - center(2))**2 + (x3p(k) - center(3))**2 .LE. radius**2) THEN
                fn_val     = fn_val + rhs(i,j,k,1)*dx1u(i)*dx2p(j)*dx3p(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    END IF
  
  !=== in y-direction =======================================================================================
  ELSE IF (dir .EQ. 2) THEN
    !--- 3EWK -> return flow --------------------------------------------------------------------------------
    IF (WK_type .EQ. 3) THEN
      IF ((pos .GE. x2p(S2p)) .AND. (pos .LE. x2p(N2p))) THEN ! check if desired position is in process
        x2loc = MINLOC(ABS(x2p - pos), 1) !< finds the index closest to the specified location

        DO k = S3p, N3p
          DO j = x2loc, x2loc
            DO i = S1p, N1p
              IF ( (x1p(i) - center(1))**2 + (x3p(k) - center(3))**2 .LE. radius**2) THEN
                fn_val     = fn_val + work2(i,j,k)*dx1p(i)*dx2p(j)*dx3p(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    !--- 4EWK -> return flow rate ---------------------------------------------------------------------------
    ELSE IF (WK_type .EQ. 4) THEN
      IF ((pos .GE. x2v(S22B)) .AND. (pos .LE. x2v(N22B))) THEN ! check if desired position is in process
        x2loc = MINLOC(ABS(x2v - pos), 1) !< finds the index closest to the specified location

        DO k = S32B, N32B
          DO j = x2loc, x2loc
            DO i = S12B, N12B
              IF ( (x1p(i) - center(1))**2 + (x3p(k) - center(3))**2 .LE. radius**2) THEN
                fn_val     = fn_val + rhs(i,j,k,2)*dx1p(i)*dx2v(j)*dx3p(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    END IF

  !=== in z-direction =======================================================================================
  ELSE IF (dir .EQ. 3) THEN
    !--- 3EWK -> return flow --------------------------------------------------------------------------------
    IF (WK_type .EQ. 3) THEN
      IF ((pos .GE. x3p(S3p)) .AND. (pos .LE. x3p(N3p))) THEN ! check if desired position is in process
        x3loc = MINLOC(ABS(x3p - pos), 1) !< finds the index closest to the specified location

        DO k = x3loc, x3loc
          DO j = S2p, N2p
            DO i = S1p, N1p
              IF ( (x1p(i) - center(1))**2 + (x2p(j) - center(2))**2 .LE. radius**2) THEN
                fn_val     = fn_val + work3(i,j,k)*dx1p(i)*dx2p(j)*dx3p(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    !--- 4EWK -> return flow rate ---------------------------------------------------------------------------
    ELSE IF (WK_type .EQ. 4) THEN
      IF ((pos .GE. x3w(S33B)) .AND. (pos .LE. x3w(N33B))) THEN ! check if desired position is in process
        x3loc = MINLOC(ABS(x3w - pos), 1) !< finds the index closest to the specified location

        DO k = x3loc, x3loc
          DO j = S23B, N23B
            DO i = S13B, N13B
              IF ( (x1p(i) - center(1))**2 + (x2p(j) - center(2))**2 .LE. radius**2) THEN
                fn_val     = fn_val + rhs(i,j,k,3)*dx1p(i)*dx2p(j)*dx3w(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    END IF

  !===========================================================================================================
  ELSE
     WRITE(*,*)'WARNING: something went wrong with calculating the flow. Aborting...'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF

  END FUNCTION flow_and_rate3D



  FUNCTION fitted_flow(xx) RESULT(fn_val)

  IMPLICIT NONE

  REAL                   ::  alpha1, beta1, gamma1
  REAL                   ::  alpha2, beta2, gamma2
  REAL                   ::  fn_val
  REAL, INTENT(IN)       ::  xx

  alpha1 = 0.0003117
  beta1  = 0.3399
  gamma1 = 0.04485
  alpha2 = 0.0004466
  beta2  = 0.2415
  gamma2 = 0.07956
  
  fn_val = alpha1*EXP(-((xx-beta1)/gamma1)**2) + alpha2*EXP(-((xx-beta2)/gamma2)**2)

  END FUNCTION fitted_flow

  
  
  
  FUNCTION fitted_pressure(xx) RESULT(fn_val)
 
  IMPLICIT NONE
 
  REAL                   :: alpha1, beta1, gamma1
  REAL                   :: alpha2, beta2, gamma2
  REAL                   :: fn_val
  REAL, INTENT(IN)       :: xx

  alpha1 = 1.647e04
  beta1  = 0.2898
  gamma1 = 0.2067
  alpha2 = 2.637e16
  beta2  = 8.607
  gamma2 = 1.387

  fn_val = alpha1*EXP(-((xx-beta1)/gamma1)**2) + alpha2*EXP(-((xx-beta2)/gamma2)**2)
  
  END FUNCTION fitted_pressure



  SUBROUTINE apply_fringe_forcing

  IMPLICIT NONE

  INTEGER            :: i,j,k,ii,jj,kk,m,n
  INTEGER, DIMENSION(2,3) :: bounds
  REAL :: vel_voxel(1:3), pulse,time_inst,time_2,time_1
  INTEGER :: ph_1,ph_2

  bounds(1,1) = lbound(mri_inst%mri%velocity_mean%array,3)
  bounds(2,1) = ubound(mri_inst%mri%velocity_mean%array,3)
  bounds(1,2) = lbound(mri_inst%mri%velocity_mean%array,4)
  bounds(2,2) = ubound(mri_inst%mri%velocity_mean%array,4)
  bounds(1,3) = lbound(mri_inst%mri%velocity_mean%array,5)
  bounds(2,3) = ubound(mri_inst%mri%velocity_mean%array,5)

  phase = mod(write_kalm_count,intervals) + 1
  pulse = write_kalm_count/intervals
  time_inst = (time+dtime) - pulse*kalman_mri_input_attr_t_heart_cycle_period 
  ph_2 = phase
  time_2 = mri_inst%mri%t_coordinates(ph_2)
  if (phase.eq.1) then
     ph_1 = intervals
     time_1 = mri_inst%mri%t_coordinates(ph_1) - kalman_mri_input_attr_t_heart_cycle_period
  else
     ph_1 = phase-1
     time_1 = mri_inst%mri%t_coordinates(ph_1)
  end if

  DO k = bounds(1,3), bounds(2,3)
    DO j = bounds(1,2), bounds(2,2)
      DO i = bounds(1,1), bounds(2,1)
        if (mri_inst%mri%segmentation_prob%array(phase,i,j,k) .ge. 0.0 ) then

          m =      i-bounds(1,1)
          m = m + (j-bounds(1,2))*size(mri_inst%mri%velocity_mean%array,3)
          m = m + (k-bounds(1,3))*size(mri_inst%mri%velocity_mean%array,3)*size(mri_inst%mri%velocity_mean%array,4)
 
          DO kk = (k-1)*kalman_num_spatial_refinements(3)+1,k*kalman_num_spatial_refinements(3)+1
            DO jj = (j-1)*kalman_num_spatial_refinements(2)+1,j*kalman_num_spatial_refinements(2)+1
              DO ii = (i-1)*kalman_num_spatial_refinements(1)+1,i*kalman_num_spatial_refinements(1)+1
                n =      ii-((bounds(1,1)-1)*kalman_num_spatial_refinements(1)+1)
                n = n + (jj-((bounds(1,2)-1)*kalman_num_spatial_refinements(2)+1))* &
                        (size(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)+1)
                n = n + (kk-((bounds(1,3)-1)*kalman_num_spatial_refinements(3)+1))* &
                        (size(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)+1) * &
                        (size(mri_inst%mri%velocity_mean%array,4)*kalman_num_spatial_refinements(2)+1)

                vel_voxel(1:3) =      (time_inst-time_1)/(time_2-time_1) *mri_inst%mri%velocity_mean%array(1:3,ph_2,i,j,k) + &
                                 (1.0-(time_inst-time_1)/(time_2-time_1))*mri_inst%mri%velocity_mean%array(1:3,ph_1,i,j,k)

                nl(ii,jj,kk,1) = nl(ii,jj,kk,1)-(vel_voxel(1)-vel(ii,jj,kk,1))/(dtime*aRK(substep))/wgt_interp(ii,jj,kk)
                nl(ii,jj,kk,2) = nl(ii,jj,kk,2)-(vel_voxel(2)-vel(ii,jj,kk,2))/(dtime*aRK(substep))/wgt_interp(ii,jj,kk)
                nl(ii,jj,kk,3) = nl(ii,jj,kk,3)-(vel_voxel(3)-vel(ii,jj,kk,3))/(dtime*aRK(substep))/wgt_interp(ii,jj,kk)

              END DO
            END DO
          END DO

        end if
      END DO
    END DO
  END DO

  END SUBROUTINE apply_fringe_forcing


  SUBROUTINE periodic_pressure(time, p)

  IMPLICIT NONE

  REAL, INTENT(IN)    ::  time
  REAL, INTENT(OUT)   ::  p
  REAL, PARAMETER     ::  pi = 2.*ABS(ACOS(0.))


  p = p_amp * COS(p_freq * 2. * pi * time)

  END SUBROUTINE periodic_pressure



  SUBROUTINE apply_windkessel_loading

  IMPLICIT NONE

  INTEGER            :: i,j,k
  REAL               :: lamb_fringe
  INTEGER            :: x1loc, x2loc, x3loc

  !--- x-direction ------------------------------------------------------------------------------------------
  IF (WK_flow_dir .EQ. 1) THEN
    DO k = S31, N31
      DO j = S21, N21
        DO i = S11, N11
          CALL fringe_coeff(WK_frge_amp , WK_frge_start , WK_frge_end , WK_frge_rise , &
                WK_frge_fall , x1u(i) , lamb_fringe)

          IF ( dimens .NE. 3) WK_flow_center(3) = x3p(k) ! transforms circle equation into linear interval
          IF ( ((x2p(j) - WK_flow_center(2))**2 + (x3p(k) - WK_flow_center(3))**2 .LE. WK_flow_radius**2 ) ) THEN
                 nl(i,j,k,1) = nl(i,j,k,1) + lamb_fringe*( WK_pre(1) ) &
                             / ( WK_frge_end - WK_frge_start ) &
                             / ( rho_fluid * U_ref**2. )
!                             / ((Re**2.)*(mu_fluid**2.))*(rho_fluid*(L_ref**2.))
          END IF
        END DO
      END DO
    END DO
  END IF
  !--- y-direction ------------------------------------------------------------------------------------------
  IF (WK_flow_dir .EQ. 2) THEN
    DO k = S32, N32
      DO j = S22, N22
        DO i = S12, N12
          CALL fringe_coeff(WK_frge_amp , WK_frge_start , WK_frge_end , WK_frge_rise , &
                WK_frge_fall , x2v(j) , lamb_fringe)

          IF ( dimens .NE. 3) WK_flow_center(3) = x3p(k) ! transforms circle equation into linear interval
          IF ( ((x1p(i) - WK_flow_center(1))**2 + (x3p(k) - WK_flow_center(3))**2 .LE. WK_flow_radius**2 ) ) THEN
                 nl(i,j,k,2) = nl(i,j,k,2) + lamb_fringe*( WK_pre(1) ) &
                             / ( WK_frge_end - WK_frge_start ) &
                             / ( rho_fluid * U_ref**2. )
!                             / ((Re**2.)*(mu_fluid**2.))*(rho_fluid*(L_ref**2.))
          END IF
        END DO
      END DO
    END DO
  END IF
  !--- z-direction ------------------------------------------------------------------------------------------
  IF (WK_flow_dir .EQ. 3 .AND. dimens .EQ. 3) THEN
    DO k = S33, N33
      DO j = S23, N23
        DO i = S13, N13
          CALL fringe_coeff(WK_frge_amp , WK_frge_start , WK_frge_end , WK_frge_rise , &
                WK_frge_fall , x3w(k) , lamb_fringe)

          IF ( ((x1p(i) - WK_flow_center(1))**2 + (x2p(j) - WK_flow_center(2))**2 .LE. WK_flow_radius**2 ) ) THEN
                 nl(i,j,k,3) = nl(i,j,k,3) + lamb_fringe*( WK_pre(1) ) & 
                             /( WK_frge_end - WK_frge_start ) &
                             / ( rho_fluid * U_ref**2. )
!                             / ((Re**2.)*(mu_fluid**2.))*(rho_fluid*(L_ref**2.))
          END IF
        END DO
      END DO
    END DO
  END IF
  !----------------------------------------------------------------------------------------------------------


  END SUBROUTINE apply_windkessel_loading





  !> function that provides a smoothed step function (erf) from 0 (if input<-cutoff=-3) to 1 (if input>cutoff=3)
  !! @param[in] xx coordinate
  FUNCTION interface(xx) RESULT(fn_val)
  ! (sample function)
  
  IMPLICIT NONE
  
  REAL   , INTENT(in   ) ::  xx
  REAL                   ::  fn_val
  REAL                   ::  cutoff
  
  
  ! underflow ...
  cutoff = 3.
  
  
  IF      (xx .GT.  cutoff) THEN
     fn_val = 1.
  ELSE IF (xx .LT. -cutoff) THEN
     fn_val = 0.
  ELSE
     fn_val = 0.5*(1.+erf(SQRT(pi)*xx))
  END IF
  
  
  END FUNCTION interface
  
  
  
  
  
  
  
  
  
  
  !> subroutine that transforms and maps coordinates to user-specified distributions
  !! @param[in] Lmax physical domain extent
  !! @param[in] iimax Number of gridpoints in the spacial dimension
  !! @param[in] ii0L lower fix bound of coordinate transform
  !! @param[in] ii0U upper fix bound of coordinate transform
  !! @param[in] ii current discrete coordinate point
  !! @param[out] xx transformed and mapped coordinate
  !! @param[out] dx immediate transformed and mapped coordinate step
  SUBROUTINE coord_tanh(Lmax,iimax,ii0L,ii0U,ii,xx,dx)
  ! (sample subroutine)
  
  IMPLICIT NONE
  
  REAL   , INTENT(in)    ::  Lmax
  REAL   , INTENT(in)    ::  iimax
  
  REAL   , INTENT(in)    ::  ii0L
  REAL   , INTENT(in)    ::  ii0U
  
  REAL   , INTENT(in)    ::  ii
  REAL   , INTENT(out)   ::  xx
  REAL   , INTENT(out)   ::  dx
  
  REAL                   ::  yy, cmin, cmax
  
  
  IF (ii0U == ii0L) THEN
     ! equidistant grid:
     xx = (ii-1.)*Lmax/(iimax-1.)
     dx =         Lmax/(iimax-1.)
  ELSE
     cmax =  TANH(ii0U)
     cmin = -TANH(ii0L)
     
     ! coordinate transformation (index i=1. is the origin):
     ! y = (i-1.)/(imax-1.)*(cmax-cmin) + cmin
     yy = (ii-1.)/(iimax-1.)*(cmax-cmin) + cmin
     
     ! mapping funktion f(yy)
     ! x = L * (f(y)-f(cmin)) / (f(cmax)-f(cmin))
     xx = Lmax * (atanh(yy)-atanh(cmin)) / (ii0U+ii0L)
     
     ! dx/di = L / (f(cmax)-f(cmin)) * dy/di                 * df(y(i))/dy
     !       = L / (f(cmax)-f(cmin)) * (cmax-cmin)/(imax-1.) * df(y(i))/dy
     dx = Lmax / (atanh(cmax)-atanh(cmin)) * (cmax-cmin)/(iimax-1.) * 1./(1.-yy**2)
  END IF
  
  
  END SUBROUTINE coord_tanh




  SUBROUTINE coord_gausslobatto(Lmax, iimax, ii, xi, xx)

  IMPLICIT NONE

  REAL, INTENT(in)    ::  Lmax      ! Physical Domain extent, i.e. L1, L2, or L3 
  REAL, INTENT(in)    ::  iimax     ! Number of gridpoints in spatial direction
  REAL, INTENT(in)    ::  ii        ! Coordinate index, starts at ii = 1, ..., Mi

  REAL, INTENT(out)   ::  xx        ! Coordinate x_i of the grid from 0 to Lmax

  REAL                ::  yy        ! Gauss-Lobatto Coordinate from 0 to 2

  REAL                ::  xi        ! Parameter for refinement level (0: Gauss lobatto points, inf: equidistant points)
  REAL                ::  kappa
  REAL, PARAMETER     ::  pi = 2.*ABS(ACOS(0.))

  kappa = pi/(iimax - 1 + 2*xi)

  ! Gauss-Lobatto Points from 0 to 2 (as in Henniger 2010)
  yy = 1. - cos(kappa*(xi + ii - 1.))/cos(kappa*xi)

  ! Mapping of the Gauss-Lobatto points from 0 to Lmax
  xx = yy/2.*Lmax


!  REAL   , INTENT(in)    ::  Lmax
!  REAL   , INTENT(in)    ::  iimax
!  REAL   , INTENT(in)    ::  Lhalf
!  REAL   , INTENT(in)    ::  Lamp
!
!  REAL   , INTENT(in)    ::  ii
!  REAL   , INTENT(out)   ::  xx
!  REAL   , INTENT(out)   ::  dx
!
!  REAL                   ::  A, B, eta, detadi
!  REAL   , PARAMETER     ::  pi = 2.*ABS(ACOS(0.))
!
!
!  IF (1==0) THEN
!  ! create Gauss-Lobatto points within [-1,1]
!  eta  = COS( (ii - 1.)*pi/(iimax - 1.) )
!  detadi = - SIN( (ii - 1.)*pi/(iimax - 1.) )*pi / (iimax - 1.)
!  END IF
!
!  IF (acos_yes) THEN
!    ! create something like "inverted" Gauss-Lobatto points
!    ! with refinement in the middle within [-1,1]
!    eta    = 1. - 2.*ACOS(-1. + 2.*(ii-1.)/(iimax-1.))/pi
!    detadi = 4./( pi * (iimax-1.) * SQRT(1. - (2.*(ii-1.)/(iimax-1.) - 1.)**2. ) + EPSILON(ii))
!
!    ! map them to [0, Lmax]
!    B = Lmax/(Lmax - 2.*Lhalf + EPSILON(Lhalf))
!    A = B * Lhalf
!  
!    xx = A*(1. + eta)/(B - eta)
!    dx = A*(B + 1.)*detadi/( (B-eta)**2. )
!  ELSE
!    ! create linear points within [0,1]
!    eta    = (ii - 1.)/(iimax - 1.)
!    detadi = 1./(iimax - 1.)
!
!    ! map them to [0, Lmax] and create 3rd order polynomial
!    xx = Lmax*eta + Lamp*(Lhalf - Lmax*eta)*(1.-eta)*eta
!    dx = detadi*(Lmax - 2.*Lamp*(Lhalf+Lmax)*eta + Lamp*Lhalf + 3.*Lamp*Lmax*eta**2.)
!  END IF
  

  END SUBROUTINE coord_gausslobatto
  
  
  SUBROUTINE load_vel_initcond_hdf5

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER       ::  filename_initcond_velx = 'velX_initcond'
  CHARACTER(LEN=*), PARAMETER       ::  filename_initcond_vely = 'velY_initcond'
  CHARACTER(LEN=*), PARAMETER       ::  filename_initcond_velz = 'velZ_initcond'

  REAL                              ::  vel_inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  LOGICAL                           ::  attr_yes
  INTEGER(HID_T)                    ::  file_id, plist_id, dset_id
  INTEGER(HID_T)                    ::  memtypeREAL, memtypeINT

  LOGICAL                           ::  regular_grid_yes_initcond
  REAL                              ::  grid_refinemnt_xi_initcond


  !--- Load initial conditions from hdf5-File ----------------------------------------------------------------

  IF(rank == 0 .AND. write_stout_yes) WRITE(*,'(a)') 'Start reading velocity initial conditions ...'

  !--- Initial conditions are in the pressure-coordinates
  !--- Need to be interpolated to velocity-coordinates
  !--- exchange() might not strictly be necessary

  !--- Load X-Velocity ---
  CALL read2_hdf(filename_initcond_velx,'velX_initcond',S1p,S2p,S3p,N1p,N2p,N3p,0,vel(:,:,:,1))
  vel_inter = 0.
  CALL exchange(1,0,vel(b1L,b2L,b3L,1))
  CALL exchange(2,0,vel(b1L,b2L,b3L,1))
  CALL exchange(3,0,vel(b1L,b2L,b3L,1))
  CALL interpolate2_pre_vel(.FALSE.,1,vel(:,:,:,1),vel_inter)
  vel(:,:,:,1) = vel_inter(:,:,:)
  CALL exchange(1,1,vel(b1L,b2L,b3L,1))
  CALL exchange(2,1,vel(b1L,b2L,b3L,1))
  CALL exchange(3,1,vel(b1L,b2L,b3L,1))
  !--- Load Y-Velocity ---
  CALL read2_hdf(filename_initcond_vely,'velY_initcond',S1p,S2p,S3p,N1p,N2p,N3p,0,vel(:,:,:,2))
  vel_inter = 0.
  CALL exchange(1,0,vel(b1L,b2L,b3L,2))
  CALL exchange(2,0,vel(b1L,b2L,b3L,2))
  CALL exchange(3,0,vel(b1L,b2L,b3L,2))
  CALL interpolate2_pre_vel(.FALSE.,2,vel(:,:,:,2),vel_inter)
  vel(:,:,:,2) = vel_inter(:,:,:)
  CALL exchange(1,2,vel(b1L,b2L,b3L,2))
  CALL exchange(2,2,vel(b1L,b2L,b3L,2))
  CALL exchange(3,2,vel(b1L,b2L,b3L,2))
  !--- Load Z-Velocity ---
  CALL read2_hdf(filename_initcond_velz,'velZ_initcond',S1p,S2p,S3p,N1p,N2p,N3p,0,vel(:,:,:,3))
  vel_inter = 0.
  CALL exchange(1,0,vel(b1L,b2L,b3L,3))
  CALL exchange(2,0,vel(b1L,b2L,b3L,3))
  CALL exchange(3,0,vel(b1L,b2L,b3L,3))
  CALL interpolate2_pre_vel(.FALSE.,3,vel(:,:,:,3),vel_inter)
  vel(:,:,:,3) = vel_inter(:,:,:)
  CALL exchange(1,3,vel(b1L,b2L,b3L,3))
  CALL exchange(2,3,vel(b1L,b2L,b3L,3))
  CALL exchange(3,3,vel(b1L,b2L,b3L,3))
  !--- Check if initial conditions use the same grid as IMPACT -----------------------------------------------
  !--- We check only in the x-Velocity...
  IF (rank == 0) WRITE(*,'(a)') 'Check if initconds have same grid as used in IMPACT...'
  !--- Copied from mod_inout.f90 SUBROUTINE read2_hdf and slightly adapted ... -------------------------------
  !===========================================================================================================
  !=== File oeffnen ==========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)

  ! Open the file collectively
  CALL h5fopen_f(filename_initcond_velx//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)

  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open file '//filename_initcond_velx//'.h5 !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF

  CALL h5pclose_f(plist_id,herror)
  !===========================================================================================================

  !===========================================================================================================
  !=== Attribute lesen =======================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,'velX_initcond',dset_id,herror)

  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open dataset velX_initcond !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  !-----------------------------------------------------------------------------------------------------------
  ! Read the attributes:
  attr_yes = .TRUE.
  CALL read_hdf_infoREAL(1,.TRUE. ,attr_yes,'grid_refinemnt_xi' ,scalar=grid_refinemnt_xi_initcond )
  CALL read_hdf_infoLOG (1,.TRUE. ,attr_yes,'regular_grid_yes',scalar=regular_grid_yes_initcond)
  !-----------------------------------------------------------------------------------------------------------
  CALL h5dclose_f(dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------

  !--- Check if parameters agree...
  IF(regular_grid_yes_initcond .NEQV. refine_2dir_yes .AND. grid_refinemnt_xi_initcond == refine_2dir_xi) THEN
    IF(rank == 0) WRITE(*,'(a)') 'OK...'
  ELSE
    IF(rank == 0) WRITE(*,'(a)') 'WARNING: not the same grid used for generating the inital-conditions!'
  END IF


  IF(rank == 0 .AND. write_stout_yes) THEN
    WRITE(*,'(a)') 'Finished reading velocity initial conditions ...'
    WRITE(*,'(a)') ''
  END IF




  END SUBROUTINE load_vel_initcond_hdf5

  
  
  
  
  
  
  
  
  !> function that calculates the inverse of the hyperbolic tangent. This function has not
  !! been made an intrinsic fortran function until the 2008 standard. The way it is calculated
  !! here, only the real parts are corresponding.
  !! @warning beware of the imaginary part if ever needed!
  FUNCTION atanh(x) RESULT(fn_val)
  ! (sample function)
  
  IMPLICIT NONE
  
  REAL   , INTENT(IN)  :: x
  REAL                 :: fn_val
  
  
  IF (x == 0.) THEN
     fn_val = 0.
     RETURN
  END IF
  
  fn_val = 0.5* LOG((1.+x)/(1.-x))
  
  
  END FUNCTION atanh
  
  
  
  
  
  
  
  
  
  
  !> Evaluation of the real error function. Based upon a Fortran 66 routine in the Naval Surface Warfare
  !! Center's Mathematics Library (1993 version). Adapten by Alan.Miller@vic.cmis.csiro.au.
  FUNCTION erf(x) RESULT(fn_val)
  ! (sample function)
  
  !-----------------------------------------------------------------------
  !             EVALUATION OF THE REAL ERROR FUNCTION
  ! Based upon a Fortran 66 routine in the Naval Surface Warfare Center's
  ! Mathematics Library (1993 version).
  ! Adapted by Alan.Miller @ vic.cmis.csiro.au
  !-----------------------------------------------------------------------
  
  IMPLICIT NONE
  
  REAL   , INTENT(IN)  ::  x
  REAL                 ::  fn_val
  
  REAL   , PARAMETER   ::  c = .564189583547756, one = 1.0, half = 0.5, zero = 0.0
  REAL   , PARAMETER   ::  &
           a(5) = (/  7.71058495001320D-05, -1.33733772997339D-03, 3.23076579225834D-02, 4.79137145607681D-02, 1.28379167095513D-01 /),  &
           b(3) = (/  3.01048631703895D-03,  5.38971687740286D-02, 3.75795757275549D-01 /),  &
           p(8) = (/ -1.36864857382717D-07,  5.64195517478974D-01, 7.21175825088309D+00, 4.31622272220567D+01, 1.52989285046940D+02, 3.39320816734344D+02, 4.51918953711873D+02, 3.00459261020162D+02 /), &
           q(8) = (/  1.00000000000000D+00,  1.27827273196294D+01, 7.70001529352295D+01, 2.77585444743988D+02, 6.38980264465631D+02, 9.31354094850610D+02, 7.90950925327898D+02, 3.00459260956983D+02 /), &
           r(5) = (/  2.10144126479064D+00,  2.62370141675169D+01, 2.13688200555087D+01, 4.65807828718470D+00, 2.82094791773523D-01 /),  &
           s(4) = (/  9.41537750555460D+01,  1.87114811799590D+02, 9.90191814623914D+01, 1.80124575948747D+01 /)
  REAL                 ::  ax, bot, t, top, x2
  
  
  ax = ABS(x)
  
  IF (ax <= half) THEN
     t = x*x
     top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + one
     bot = ((b(1)*t + b(2))*t + b(3))*t + one
     fn_val = x*(top/bot)
     RETURN
  END IF
  
  IF (ax <= 4.0) THEN
     top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax  &
           + p(6))*ax + p(7))*ax + p(8)
     bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax  &
           + q(6))*ax + q(7))*ax + q(8)
     fn_val = half + (half - EXP(-x*x)*top/bot)
     IF (x < zero) fn_val = -fn_val
     RETURN
  END IF
  
  IF (ax < 5.8) THEN
     x2 = x*x
     t = one / x2
     top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
     bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + one
     fn_val = (c - top/(x2*bot)) / ax
     fn_val = half + (half - EXP(-x2)*fn_val)
     IF (x < zero) fn_val = -fn_val
     RETURN
  END IF
  
  fn_val = SIGN(one, x)
  
  
  END FUNCTION erf
 
  
  SUBROUTINE init_hdf5()
  
  USE HDF5

  IMPLICIT NONE

  CALL h5open_f(herror)

  END SUBROUTINE init_hdf5


  SUBROUTINE finl_hdf5()
  
  USE HDF5

  IMPLICIT NONE

  CALL h5close_f(herror)

  END SUBROUTINE finl_hdf5



  SUBROUTINE init_mpi()
  
  IMPLICIT NONE

  CALL MPI_INIT(merror)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,merror)

  END SUBROUTINE init_mpi


  SUBROUTINE finl_mpi()
  
  IMPLICIT NONE

  CALL MPI_FINALIZE(merror)

  END SUBROUTINE finl_mpi

  SUBROUTINE print_fcomm_size()

  IMPLICIT NONE

  INTEGER :: fcomm_size
  INTEGER :: fcomm_rank
  INTEGER :: name_len
  CHARACTER(len=MPI_MAX_OBJECT_NAME) :: comm_name

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,fcomm_size,merror)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,fcomm_rank,merror)
  CALL MPI_COMM_GET_NAME(MPI_COMM_WORLD,comm_name,name_len,merror)

  print*, 'MPI_COMM_WORLD with name ',comm_name,' in Fortran has size: ', fcomm_size, ' in process ',fcomm_rank,'.'

  END SUBROUTINE print_fcomm_size

  SUBROUTINE mpi_bcast_fort()

  IMPLICIT NONE

  CALL MPI_BCAST(finish_yes,1,MPI_LOGICAL,0,COMM_CART,merror)


  END SUBROUTINE mpi_bcast_fort


  SUBROUTINE check_node_ids

  IMPLICIT NONE

  INTEGER                ::  msize
  INTEGER                ::  i,j,k
  CHARACTER(LEN=2)       ::  schar1,schar2,schar3
  INTEGER                ::  ii,jj,kk,nid,gnid
  INTEGER                ::  nn1l, nn1u, nn2l, nn2u, nn3l, nn3u
  INTEGER                ::  dummy1, dummy2
  INTEGER                ::  nn1, nn2, nn3
  INTEGER                ::  dir
  LOGICAL                ::  boundary_yes, vel_grid_yes
  INTEGER                ::  n_local, n_ghost
  INTEGER                ::  n_elem_local

  ii = n1p
  jj = 1
  kk = 1
  dir = 1
  boundary_yes = .TRUE.
  vel_grid_yes = .TRUE.

  CALL get_block_dims(nn1, nn2, nn3, .FALSE., .FALSE., -1)
  
  CALL compute_and_store_global_block_sizes(.FALSE.,.FALSE.,-1)

  write(*,*) 'process',rank,'nn1:',nn1,'N1p-S1p+1',N1p-S1p+1
  write(*,*) 'process',rank,'nn2:',nn2,'N2p-S2p+1',N2p-S2p+1
  write(*,*) 'process',rank,'nn3:',nn3,'N3p-S3p+1',N3p-S3p+1

  CALL get_block_dims(nn1, nn2, nn3, .TRUE., .FALSE., 1)

  write(*,*) 'process',rank,'nn1:',nn1,'N11-S11+1',N11-S11+1
  write(*,*) 'process',rank,'nn2:',nn2,'N21-S21+1',N21-S21+1
  write(*,*) 'process',rank,'nn3:',nn3,'N31-S31+1',N31-S31+1

  CALL get_block_dims(nn1, nn2, nn3, .TRUE., .TRUE., 1)

  write(*,*) 'process',rank,'nn1:',nn1,'N11B-S11B+1',N11B-S11B+1
  write(*,*) 'process',rank,'nn2:',nn2,'N21B-S21B+1',N21B-S21B+1
  write(*,*) 'process',rank,'nn3:',nn3,'N31B-S31B+1',N31B-S31B+1

  !--- Get the global number of nodes on which the values are stored
  CALL get_boundary_block_dims(nn1l, nn1u, nn2l, nn2u, nn3l, nn3u, vel_grid_yes, boundary_yes, dir)
  nn1 = nn1l + (NB1 - 2)*(N1 - 1) + nn1u
  nn2 = nn2l + (NB2 - 2)*(N2 - 1) + nn2u
  nn3 = nn3l + (NB3 - 2)*(N3 - 1) + nn3u

  !--- If in any direction only one block exists override dimension with local one
  IF (NB1 .EQ. 1) THEN
    CALL get_block_dims(nn1, dummy1, dummy2, vel_grid_yes, boundary_yes, dir)
  ENDIF
  IF (NB2 .EQ. 1) THEN
    CALL get_block_dims(dummy1, nn2, dummy2, vel_grid_yes, boundary_yes, dir)
  ENDIF
  IF (NB3 .EQ. 1) THEN
    CALL get_block_dims(dummy1, dummy2, nn3, vel_grid_yes, boundary_yes, dir)
  END IF

  

  CALL num_to_string(2,iB(1,1),schar1)
  CALL num_to_string(2,iB(2,1),schar2)
  CALL num_to_string(2,iB(3,1),schar3)
  OPEN(rank, FILE = "node_ids_proc_"//schar1//"_"//schar2//"_"//schar3//".txt",STATUS='UNKNOWN')
  OPEN(rank+100, FILE = "node_ids_local_cont_proc_"//schar1//"_"//schar2//"_"//schar3//".txt",STATUS='UNKNOWN')
  OPEN(rank+200, FILE = "local_node_ids_local_cont_proc_"//schar1//"_"//schar2//"_"//schar3//".txt",STATUS='UNKNOWN')


  IF (dimens .EQ. 3) THEN
    DO k = s3p-1, n3p
      DO j = s2p-1, n2p
        DO i = s1p-1, n1p

        !=== local coordinates with ghost cells ===
        CALL local_cart2id_loc_con(nid,i,j,k,n_local,n_ghost)
        write(rank+200,'(3i4,1a,1i8)') i,j,k, '-->', nid
        CALL local_id2cart_loc_con(nid,ii,jj,kk,n_local,n_ghost)
        write(rank+200,'(1i8,1a,3i4)') nid, '-->', ii,jj,kk
        write(rank+200,'(a)') '============================='
        FLUSH(rank+200)
        END DO
      END DO
    END DO
  ELSE
    DO k = s3p, n3p
      DO j = s2p-1, n2p
        DO i = s1p-1, n1p

        !=== local coordinates with ghost cells ===
        CALL local_cart2id_loc_con(nid,i,j,k,n_local,n_ghost)
        write(rank+200,'(3i4,1a,1i8)') i,j,k, '-->', nid
        CALL local_id2cart_loc_con(nid,ii,jj,kk,n_local,n_ghost)
        write(rank+200,'(1i8,1a,3i4)') nid, '-->', ii,jj,kk
        write(rank+200,'(a)') '============================='
        FLUSH(rank+200)
        END DO
      END DO
    END DO
  END IF



  DO k = s3p, n3p
    DO j = s2p, n2p
      DO i = s1p, n1p

      !=== global coordinates locally contiguous ===
      CALL global_cart2id_loc_con(nid,i,j,k)
      write(rank+100,'(3i4,1a,1i8)') i,j,k, '-->', nid
      CALL global_id2cart_loc_con(nid,ii,jj,kk)
      write(rank+100,'(1i8,1a,3i4)') nid, '-->', ii,jj,kk
      write(rank+100,'(a)') '============================='
      !=== output all nodes and their global ids ===
      CALL global_cart2id_glob_con(nid,i,j,k,.FALSE.,.FALSE.,dir)
      write(rank,'(3i4,1a,1i8)') i,j,k, '-->', nid
      !write(*,*) '=== Pressure Grid ==='
      !write(*,*) 'Process ',rank,': (',iB(1,1),',',iB(2,1),',',iB(3,1),'), Node (',i,',',j,',',k,') --> NodeID: ',nid
      CALL global_id2cart_glob_con(nid,ii,jj,kk,.FALSE.,.FALSE.,dir)
      write(rank,'(1i8,1a,3i4)') nid, '-->', ii,jj,kk
      write(rank,'(a)') '============================='
      !write(*,*) 'Process ',rank,': (',iB(1,1),',',iB(2,1),',',iB(3,1),'), NodeID: ',nid,' --> Node (',i,',',j,',',k,')' 
      !write(*,*) 'Node has ID: ', nid,' of a total ', nn1, 'x', nn2, 'x',nn3 , '=',nn1*nn2*nn3
      FLUSH(rank)
      FLUSH(rank+100)

      END DO
    END DO
  END DO

  CLOSE(rank)
  CLOSE(rank+100)
  CLOSE(rank+200)

  !=== Do Element node checking ============================
  IF (dimens .EQ. 3) THEN
    IF (elem_type .EQ. 1) n_elem_local = total_n_local_tet_elements()
    IF (elem_type .EQ. 2) n_elem_local = total_n_local_hex_elements()
  ELSE
    IF (elem_type .EQ. 1) n_elem_local = total_n_local_tri_elements()
    IF (elem_type .EQ. 2) n_elem_local = total_n_local_quad_elements()
  END IF
  OPEN(rank+300, FILE = "local_elem_node_ids_local_cont_proc_"//schar1//"_"//schar2//"_"//schar3//".txt",STATUS='UNKNOWN')
  DO i = 0, n_elem_local-1
    IF (elem_type .EQ. 1 .AND. dimens .EQ. 2) write(rank+300,'(1a,1i8)') "element", local_to_global_tri_elem_id(i)
    IF (elem_type .EQ. 1 .AND. dimens .EQ. 3) write(rank+300,'(1a,1i8)') "element", local_to_global_tet_elem_id(i)
    IF (elem_type .EQ. 2 .AND. dimens .EQ. 2) write(rank+300,'(1a,1i8)') "element", local_to_global_quad_elem_id(i)
    IF (elem_type .EQ. 2 .AND. dimens .EQ. 3) write(rank+300,'(1a,1i8)') "element", local_to_global_hex_elem_id(i)
    write(rank+300,'(a)') "=========="
    IF (dimens .EQ. 3) THEN
      IF(elem_type .EQ. 1) THEN
      DO j = 0, 3
        CALL local_tet_element_nodes(i,j,nid)
        CALL local_id2cart_loc_con(nid,ii,jj,kk,dummy1,dummy2)
        write(rank+300,'(1a,1i4,1a,1i8,1a,3i4)') 'node#',j,', node_id',nid,'  -->',ii,jj,kk
        CALL local_to_global_node_id_loc_con(nid, gnid)
        write(rank+300,'(1a,1i4,1a,1i8,1a,3i4)') 'node#',j,', node_id_loc',nid,' node_id_glob', gnid
      END DO
      ELSE IF (elem_type .EQ. 2) THEN
      DO j = 0, 7
        CALL local_hex_element_nodes(i,j,nid)
        CALL local_id2cart_loc_con(nid,ii,jj,kk,dummy1,dummy2)
        write(rank+300,'(1a,1i4,1a,1i8,1a,3i4)') 'node#',j,', node_id',nid,'  -->',ii,jj,kk
        CALL local_to_global_node_id_loc_con(nid, gnid)
        write(rank+300,'(1a,1i4,1a,1i8,1a,3i4)') 'node#',j,', node_id_loc',nid,' node_id_glob', gnid
      END DO
      END IF
    ELSE
      IF(elem_type .EQ. 1) THEN
      DO j = 0, 2
        CALL local_tri_element_nodes(i,j,nid)
        CALL local_id2cart_loc_con(nid,ii,jj,kk,dummy1,dummy2)
        write(rank+300,'(1a,1i4,1a,1i8,1a,3i4)') 'node#',j,', node_id',nid,'  -->',ii,jj,kk
        CALL local_to_global_node_id_loc_con(nid, gnid)
        write(rank+300,'(1a,1i4,1a,1i8,1a,3i4)') 'node#',j,', node_id_loc',nid,' node_id_glob', gnid
      END DO
      ELSE IF (elem_type .EQ. 2) THEN
      DO j = 0, 3
        CALL local_quad_element_nodes(i,j,nid)
        CALL local_id2cart_loc_con(nid,ii,jj,kk,dummy1,dummy2)
        write(rank+300,'(1a,1i4,1a,1i8,1a,3i4)') 'node#',j,', node_id',nid,'  -->',ii,jj,kk
        CALL local_to_global_node_id_loc_con(nid, gnid)
        write(rank+300,'(1a,1i4,1a,1i8,1a,3i4)') 'node#',j,', node_id_loc',nid,' node_id_glob', gnid
      END DO
      END IF
    END IF
    write(rank+300,'(a)') "=========="
    FLUSH(rank+300)
  END DO

  CLOSE(rank+300)

 !=== Do local -> global checking (nodes)=========================
 OPEN(rank+400, FILE = "local_to_global_node_ids_local_cont_proc_nodes_"//schar1//"_"//schar2//"_"//schar3//".txt",STATUS='UNKNOWN')
 IF (dimens .EQ. 3) THEN
   DO kk = s3p - 1, n3p
     DO jj = s2p - 1, n2p
       DO ii = s1p - 1, n1p
          CALL local_cart2id_loc_con(nid,ii,jj,kk,dummy1,dummy2)
          CALL local_to_global_node_id_loc_con(nid,gnid)
          write(rank+400,'(1a,1i8,1a,3i4,1a,1i8)') 'local:',nid,'; i,j,k:',ii,jj,kk,'; global conversion:',gnid
          CALL global_to_local_node_id_loc_con(nid,gnid)
          write(rank+400,'(1a,1i8,1a,3i4,1a,1i8)') 'global:',gnid,'; i,j,k:',ii,jj,kk,'; local conversion:',nid
          write(rank+400,'(a)') "==============================="
          FLUSH(rank+400)
        END DO
      END DO
    END DO
  ELSE
     DO kk = s3p, n3p
       DO jj = s2p - 1, n2p
         DO ii = s1p - 1, n1p
            CALL local_cart2id_loc_con(nid,ii,jj,kk,dummy1,dummy2)
            CALL local_to_global_node_id_loc_con(nid,gnid)
            write(rank+400,'(1a,1i8,1a,3i4,1a,1i8)') 'local:',nid,'; i,j,k:',ii,jj,kk,'; global conversion:',gnid
            CALL global_to_local_node_id_loc_con(nid,gnid)
            write(rank+400,'(1a,1i8,1a,3i4,1a,1i8)') 'global:',gnid,'; i,j,k:',ii,jj,kk,'; local conversion:',nid
            write(rank+400,'(a)') "==============================="
            FLUSH(rank+400)
          END DO
        END DO
      END DO
  END IF

  !=== Do local -> global checking (elems)=========================
  OPEN(rank+800, FILE = "local_to_global_node_ids_local_cont_proc_elems_"//schar1//"_"//schar2//"_"//schar3//".txt",STATUS='UNKNOWN')
  IF(dimens .NE. 3 .AND. elem_type .EQ. 1) THEN
	  write(rank+800,'(a)') "=== TRI elems ============================"
	  DO ii = 0, total_n_local_tri_elements()-1
	    write(rank+800,'(1a,1i8,1a,1i8,1a,1i8,1a,1i8)') 'local:', ii ,'; global conversion:', local_to_global_tri_elem_id(ii) ,'; total local:',total_n_local_tri_elements(),'; total global:', total_n_global_tri_elements()
	  END DO
  ELSE IF(dimens .NE. 3 .AND. elem_type .EQ. 2) THEN
	  write(rank+800,'(a)') "=== QUAD elems ============================"
	  DO ii = 0, total_n_local_quad_elements()-1
	    write(rank+800,'(1a,1i8,1a,1i8,1a,1i8,1a,1i8)') 'local:', ii ,'; global conversion:', local_to_global_quad_elem_id(ii) ,'; total local:',total_n_local_quad_elements(),'; total global:', total_n_global_quad_elements()
	  END DO
  END IF
  IF(dimens .EQ. 3 .AND. elem_type .EQ. 1) THEN
	  write(rank+800,'(a)') "=== TET elems ============================"
	  DO ii = 0, total_n_local_tet_elements()-1
	    write(rank+800,'(1a,1i8,1a,1i8,1a,1i8,1a,1i8)') 'local:', ii ,'; global conversion:', local_to_global_tet_elem_id(ii) ,'; total local:',total_n_local_tet_elements(),'; total global:', total_n_global_tet_elements()
	  END DO
  ELSE IF(dimens .EQ. 3 .AND. elem_type .EQ. 2) THEN
	  write(rank+800,'(a)') "=== HEX elems ============================"
	  DO ii = 0, total_n_local_hex_elements()-1
	    write(rank+800,'(1a,1i8,1a,1i8,1a,1i8,1a,1i8)') 'local:', ii ,'; global conversion:', local_to_global_hex_elem_id(ii) ,'; total local:',total_n_local_hex_elements(),'; total global:', total_n_global_hex_elements()
	  END DO
  END IF

  CLOSE(rank+400)
  CLOSE(rank+800)

  print*, 'x1p(N1p)', x1p(N1p), 'y1p(N1p)',y1p(2*N1p)

  print*, 'n_global_nodes: ', n_global_nodes()

  print* , mod(-43, 40), modulo(-43,40)
 
  print*, "process ", rank, "BC_1L ", BC_1L, "BC_2_L ", BC_2L, "BC3L ", BC_3L

  print*, "process ", rank, "BC_1L_global ", BC_1L_global, "BC_2_L_global ", BC_2L_global, "BC_3L_global ", BC_3L_global

  print*, "process ", rank, "s1p ", s1p, "s2p", s2p, "s3p ", s3p, "n1p ", n1p, "n2p", n2p, "n3p ", n3p
  END SUBROUTINE check_node_ids










  SUBROUTINE local_to_global_node_id_loc_con_allperiodic2(loc_node_id, glob_node_id)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    ::  loc_node_id
  INTEGER, INTENT(OUT)   ::  glob_node_id
  INTEGER                ::  nn1, nn2, nn3
  INTEGER                ::  nn_block
  INTEGER                ::  ii, jj, kk
  INTEGER                ::  min_current_global_id
  INTEGER                ::  n_local, n_ghost
  INTEGER                ::  to_block, from_block
  INTEGER                ::  i_to_block, j_to_block, k_to_block

  nn1 = N1 - 1
  nn2 = N2 - 1
  nn3 = N3 - 1

  nn_block = nn1 * nn2 * nn3

  i_to_block = iB(1,1)
  j_to_block = iB(2,1)
  k_to_block = iB(3,1)

  !--- convert to i, j, k
  CALL local_id2cart_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)

  !--- check if indices are greater or lesser than nip, sip and add/substract
  !--- nip.
  IF (ii .GT. N1p) THEN
    ii = ii - N1p
    i_to_block = i_to_block + 1
  ELSE IF (ii .LT. S1p) THEN
    ii = ii + N1p
    i_to_block = i_to_block - 1
  END IF
  IF (jj .GT. N2p) THEN
    jj = jj - N2p
    j_to_block = j_to_block + 1
  ELSE IF (jj .LT. S2p) THEN
    jj = jj + N2p
    j_to_block = j_to_block - 1
  END IF
  IF (kk .GT. N3p) THEN
    kk = kk - N3p
    k_to_block = k_to_block + 1
  ELSE IF (kk .LT. S3p) THEN
    kk = kk + N3p
    k_to_block = k_to_block - 1
  END IF

  !--- create local global variable ---
  ii = ii - 1
  jj = jj - 1
  kk = kk - 1
  glob_node_id = ii + jj * nn1 + kk * nn1 * nn2

  !--- augment global variable to be of global scope ---
  CALL global_cart2id_loc_con(min_current_global_id, S1p, S2p, S3p)
  from_block = block_id(iB(1,1), iB(2,1), iB(3,1))
  to_block   = block_id(i_to_block, j_to_block, k_to_block)
  glob_node_id = glob_node_id + min_current_global_id + (to_block - from_block) * nn_block

  END SUBROUTINE local_to_global_node_id_loc_con_allperiodic2



  !> This routine should work for arbitrary grids/BCs, does it work?
  SUBROUTINE local_to_global_node_id_loc_con(loc_node_id, glob_node_id)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    ::  loc_node_id
  INTEGER, INTENT(OUT)   ::  glob_node_id
  INTEGER                ::  nn1, nn2, nn3
  INTEGER                ::  nn_block
  INTEGER                ::  ii, jj, kk
  INTEGER                ::  min_current_global_id
  INTEGER                ::  n_local, n_ghost
  INTEGER                ::  to_block, from_block
  INTEGER                ::  i_to_block, j_to_block, k_to_block

  i_to_block = iB(1,1)
  j_to_block = iB(2,1)
  k_to_block = iB(3,1)

  !--- convert to i, j, k
  CALL local_id2cart_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)

  !--- check if indices are greater or lesser than nip, sip and check which
  !--- they have to be transferred to.
  !--- 1-direction
  IF (ii .GT. N1p) THEN
    i_to_block = i_to_block + 1
  ELSE IF (ii .LT. S1p) THEN
    i_to_block = i_to_block - 1
  END IF
  !--- 2-direction
  IF (jj .GT. N2p) THEN
    j_to_block = j_to_block + 1
  ELSE IF (jj .LT. S2p) THEN
    j_to_block = j_to_block - 1
  END IF
  !--- 3-direction
  IF (kk .GT. N3p) THEN
    k_to_block = k_to_block + 1
  ELSE IF (kk .LT. S3p) THEN
    k_to_block = k_to_block - 1
  END IF
  !--- get block ids
  from_block = block_id(iB(1,1), iB(2,1), iB(3,1))
  to_block   = block_id(i_to_block, j_to_block, k_to_block)

  !--- now transform indices to local global values on the destination block
  IF (ii .GT. N1p) THEN
    ii = ii - block_indices(to_block*6 + 3)
  ELSE IF (ii .LT. S1p) THEN
    ii = ii + block_indices(to_block*6 + 3)
  END IF
  IF (jj .GT. N2p) THEN
    jj = jj - block_indices(to_block*6 + 4)
  ELSE IF (jj .LT. S2p) THEN
    jj = jj + block_indices(to_block*6 + 4)
  END IF
  IF (kk .GT. N3p) THEN
    kk = kk - block_indices(to_block*6 + 5)
  ELSE IF (kk .LT. S3p) THEN
    kk = kk + block_indices(to_block*6 + 5)
  END IF
  !--- create local global variable on destination block ---
  ii = ii - 1
  jj = jj - 1
  kk = kk - 1
  nn1 = block_indices(to_block*6 + 3) - block_indices(to_block*6    ) + 1 
  nn2 = block_indices(to_block*6 + 4) - block_indices(to_block*6 + 1) + 1 
  nn3 = block_indices(to_block*6 + 5) - block_indices(to_block*6 + 2) + 1 
  glob_node_id = ii + jj * nn1 + kk * nn1 * nn2

  !--- augment global variable to be of global scope ---
  CALL global_cart2id_loc_con(min_current_global_id, S1p, S2p, S3p) !TODO
  glob_node_id = glob_node_id + min_current_global_id + global_node_difference(to_block, from_block)

  END SUBROUTINE local_to_global_node_id_loc_con







  SUBROUTINE global_to_local_node_id_loc_con_allperiodic2(loc_node_id, glob_node_id)

  IMPLICIT NONE

  INTEGER, INTENT(OUT)    :: loc_node_id
  INTEGER, INTENT(IN)     :: glob_node_id
  INTEGER                 :: nn1, nn2, nn3
  INTEGER                 :: ii, jj, kk
  INTEGER                 :: nn_block
  INTEGER                 :: min_current_global_id, max_current_global_id
  INTEGER                 :: n_local, n_ghost
  INTEGER                 :: temp_id, rel_id
  INTEGER                 :: nn_block_mult
  LOGICAL                 :: i_block_shift, j_block_shift, k_block_shift
  INTEGER                 :: current_block_id, dest_block_id
  INTEGER                 :: i_block_dest, j_block_dest, k_block_dest

  nn1 = N1 - 1
  nn2 = N2 - 1
  nn3 = N3 - 1

  nn_block = nn1 * nn2 * nn3

  !--- Get the current process' minimum and maximum id ---
  CALL global_cart2id_loc_con(min_current_global_id, S1p, S2p, S3p)
  CALL global_cart2id_loc_con(max_current_global_id, N1p, N2p, N3p)

  !--- Check whether our glob_node_id is within this process' global ids ---
  !--- If yes the conversion is straight forward ---
  IF (glob_node_id .GE. min_current_global_id .AND. glob_node_id .LE. max_current_global_id) THEN
    CALL global_id2cart_loc_con(glob_node_id, ii, jj, kk)
    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
  ELSE
  !--- relate the global id to the current block's minimum id ---
  rel_id = glob_node_id - min_current_global_id

  !--- calculate the remainder (local relative global index) ---
  temp_id = MOD(rel_id, nn_block) ! >/< 0

  !--- if temp_id is negative, account for that ---
  temp_id = MODULO(temp_id, nn_block)

  !--- now convert to i,j,k ---
  CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)

  !--- check which of the block-in-question's positions lie shifted
  !--- we have to do this because the mapping is not bijective by itself
  !--- floating point arithmetic
  i_block_shift = .FALSE.
  j_block_shift = .FALSE.
  k_block_shift = .FALSE.
  IF (rel_id .LT. 1) rel_id = rel_id - nn_block
  ! This tells us how many block we move forward(+)/backward(-)
  nn_block_mult = rel_id / nn_block ! this should be an integer!
  ! Get current block's id
  current_block_id = block_id(iB(1,1), iB(2,1), iB(3,1))
  ! Calculate the desination block_id
  dest_block_id = current_block_id + nn_block_mult
  ! Revert the destination block id to i,j,k notation
  CALL block_cart(dest_block_id, i_block_dest, j_block_dest, k_block_dest)

  if (rank == 7) then
  open(UNIT=787, FILE="periodic_test.txt", STATUS="unknown",position="append")
  write(787,*) "rel_id: ", rel_id, " temp_id: ", temp_id," nn_block_mult: ", nn_block_mult, " glob_id: ", glob_node_id
  close(787)
  end if

  ! set the shift flags depending on whether the desination block's 
  ! and the current block's indices differ
  IF (iB(1,1) .NE. i_block_dest) i_block_shift = .TRUE.
  IF (iB(2,1) .NE. j_block_dest) j_block_shift = .TRUE.
  IF (iB(3,1) .NE. k_block_dest) k_block_shift = .TRUE.


  !--- check if any of the indices match sip or nip, if they do, 
  !--- add or substract nip. In any case add 1 ---
  ii = ii + 1
  jj = jj + 1
  kk = kk + 1
  
  IF (ii .EQ. S1p) THEN
    IF (i_block_shift) ii = ii + N1p
  ELSE IF (ii .EQ. N1p) THEN
    IF (i_block_shift) ii = ii - N1p
  END IF
  IF (jj .EQ. S2p) THEN
    IF (j_block_shift) jj = jj + N2p
  ELSE IF (jj .EQ. N2p) THEN
    IF (j_block_shift) jj = jj - N2p
  END IF
  IF (kk .EQ. S3p) THEN
    IF (k_block_shift) kk = kk + N3p
  ELSE IF (kk .EQ. N3p) THEN
    IF (k_block_shift) kk = kk - N3p
  END IF

  !--- now convert back to local index ---
  CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)

  END IF

  END SUBROUTINE global_to_local_node_id_loc_con_allperiodic2


  !> This routine should work for arbitrary grids, does it work?
  SUBROUTINE global_to_local_node_id_loc_con(loc_node_id, glob_node_id)

  IMPLICIT NONE

  INTEGER, INTENT(OUT)    :: loc_node_id
  INTEGER, INTENT(IN)     :: glob_node_id
  INTEGER                 :: nn1, nn2, nn3
  INTEGER                 :: ii, jj, kk
  INTEGER                 :: nn_block
  INTEGER                 :: min_current_global_id, max_current_global_id
  INTEGER                 :: n_local, n_ghost
  INTEGER                 :: temp_id, rel_id
  INTEGER                 :: nn_block_mult
  LOGICAL                 :: i_block_shift, j_block_shift, k_block_shift
  INTEGER                 :: current_block_id, dest_block_id
  INTEGER                 :: i_block_dest, j_block_dest, k_block_dest

  !--- Get the current process' minimum and maximum id ---
  CALL global_cart2id_loc_con(min_current_global_id, S1p, S2p, S3p) !TODO
  CALL global_cart2id_loc_con(max_current_global_id, N1p, N2p, N3p) !TODO

  !--- Check whether our glob_node_id is within this process' global ids ---
  !--- If yes the conversion is straight forward ---
  IF (glob_node_id .GE. min_current_global_id .AND. glob_node_id .LE. max_current_global_id) THEN
    CALL global_id2cart_loc_con(glob_node_id, ii, jj, kk) !TODO
    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
  ELSE
  !--- relate the global id to the current block's minimum id ---
  rel_id = glob_node_id - min_current_global_id

  !--- check which of the block-in-question's positions lie shifted
  !--- we have to do this because the mapping is not bijective by itself
  !--- floating point arithmetic
  i_block_shift = .FALSE.
  j_block_shift = .FALSE.
  k_block_shift = .FALSE.
  ! IF (rel_id .LT. 1) rel_id = rel_id - nn_block
  ! Get current block's id
  current_block_id = block_id(iB(1,1), iB(2,1), iB(3,1))

  nn_block_mult = 0
  temp_id = rel_id
  DO
    IF (temp_id .LT. 0) THEN
      nn_block_mult = nn_block_mult - 1
      temp_id = temp_id + block_sizes(MODULO(current_block_id + nn_block_mult, NB1*NB2*NB3))
      IF (temp_id .GE. 0) THEN
        IF (temp_id .EQ. 0) THEN
          nn_block_mult = nn_block_mult - 1
        END IF
        EXIT
      END IF
    ELSE
      IF (temp_id .LT. block_sizes(current_block_id + nn_block_mult)) EXIT
      temp_id = temp_id - block_sizes(MODULO(current_block_id + nn_block_mult, NB1*NB2*NB3))
      nn_block_mult = nn_block_mult + 1
    END IF
  END DO

  ! This tells us how many block we move forward(+)/backward(-)
  !nn_block_mult = rel_id / nn_block ! this should be an integer!
  ! Calculate the desination block_id
  dest_block_id = MODULO(current_block_id + nn_block_mult, NB1*NB2*NB3)
  nn1 = block_indices(dest_block_id*6 + 3) - block_indices(dest_block_id*6    ) + 1 
  nn2 = block_indices(dest_block_id*6 + 4) - block_indices(dest_block_id*6 + 1) + 1 
  nn3 = block_indices(dest_block_id*6 + 5) - block_indices(dest_block_id*6 + 2) + 1 

  !--- now convert to i,j,k ---
  CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)

  ! Revert the destination block id to i,j,k notation
  CALL block_cart(dest_block_id, i_block_dest, j_block_dest, k_block_dest)

  ! set the shift flags depending on whether the desination block's 
  ! and the current block's indices differ
  IF (iB(1,1) .NE. i_block_dest) i_block_shift = .TRUE.
  IF (iB(2,1) .NE. j_block_dest) j_block_shift = .TRUE.
  IF (iB(3,1) .NE. k_block_dest) k_block_shift = .TRUE.


  !--- check if any of the indices match sip or nip, if they do, 
  !--- add or substract nip. In any case add 1 ---
  ii = ii + 1
  jj = jj + 1
  kk = kk + 1

  IF (ii .EQ. block_indices(dest_block_id*6    ) ) THEN
    IF (i_block_shift) ii = ii + block_indices(dest_block_id*6 + 3)
  ELSE IF (ii .EQ. block_indices(dest_block_id*6 + 3) ) THEN
    IF (i_block_shift) ii = ii - block_indices(dest_block_id*6 + 3)
  END IF
  IF (jj .EQ. block_indices(dest_block_id*6 + 1) ) THEN
    IF (j_block_shift) jj = jj + block_indices(dest_block_id*6 + 4) 
  ELSE IF (jj .EQ. block_indices(dest_block_id*6 + 4) ) THEN
    IF (j_block_shift) jj = jj - block_indices(dest_block_id*6 + 4)
  END IF
  IF (dimens .EQ. 3) THEN
    IF (kk .EQ. block_indices(dest_block_id*6 + 2) ) THEN
      IF (k_block_shift) kk = kk + block_indices(dest_block_id*6 + 5)
    ELSE IF (kk .EQ. block_indices(dest_block_id*6 + 5) ) THEN
      IF (k_block_shift) kk = kk - block_indices(dest_block_id*6 + 5) 
    END IF
  END IF

  !--- now convert back to local index ---
  CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)


  END IF

  END SUBROUTINE global_to_local_node_id_loc_con



  FUNCTION block_id(ib1,ib2,ib3) RESULT(fn_val)

  IMPLICIT NONE 

  INTEGER                ::  fn_val
  INTEGER, INTENT(IN)    ::  ib1, ib2, ib3
  INTEGER                ::  ii, jj, kk
 
  !--- respect periodicity ---
  ii = MODULO(ib1-1, NB1) + 1
  jj = MODULO(ib2-1, NB2) + 1
  kk = MODULO(ib3-1, NB3) + 1

  !--- column major ---
  !fn_val = (ii - 1) + (jj - 1) * NB1 + (kk - 1) * NB1 * NB2

  !--- row major ---
  fn_val = (kk - 1) + (jj - 1) * NB3 + (ii - 1) * NB2 * NB3

  END FUNCTION block_id


  SUBROUTINE block_cart(block_id, ib, jb, kb)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: block_id
  INTEGER, INTENT(OUT)   :: ib, jb, kb
  INTEGER                :: bid

  bid = block_id

  !--- column major ---
  !CALL id2cart(bid, ib, jb, kb, NB1, NB2, NB3)

  !--- row major ---
  CALL id2cart(bid, kb, jb, ib, NB3, NB2, NB1)

  ib = ib + 1
  jb = jb + 1
  kb = kb + 1

  END SUBROUTINE block_cart 


  SUBROUTINE global_cart2id_glob_con(node_id, i, j, k, vel_grid_yes, boundary_yes, dir)

  IMPLICIT NONE

  INTEGER, INTENT(IN)     ::  i,j,k
  INTEGER                 ::  ii,jj,kk
  INTEGER                 ::  nn1,nn2,nn3
  INTEGER                 ::  nn1l, nn1u, nn2l, nn2u, nn3l, nn3u
  INTEGER                 ::  dummy1, dummy2
  LOGICAL, INTENT(IN)     ::  vel_grid_yes, boundary_yes
  INTEGER, INTENT(IN)     ::  dir

  INTEGER, INTENT(OUT)    ::  node_id

  !--- Get the global number of nodes on which the values are stored
  CALL get_boundary_block_dims(nn1l, nn1u, nn2l, nn2u, nn3l, nn3u, vel_grid_yes, boundary_yes, dir)
  nn1 = nn1l + (NB1 - 2)*(N1 - 1) + nn1u
  nn2 = nn2l + (NB2 - 2)*(N2 - 1) + nn2u
  nn3 = nn3l + (NB3 - 2)*(N3 - 1) + nn3u

  !--- If in any direction only one block exists override dimension with local one
  IF (NB1 .EQ. 1) THEN
    CALL get_block_dims(nn1, dummy1, dummy2, vel_grid_yes, boundary_yes, dir)
  ENDIF
  IF (NB2 .EQ. 1) THEN
    CALL get_block_dims(dummy1, nn2, dummy2, vel_grid_yes, boundary_yes, dir)
  ENDIF
  IF (NB3 .EQ. 1) THEN
    CALL get_block_dims(dummy1, dummy2, nn3, vel_grid_yes, boundary_yes, dir)
  END IF

  !--- correct indices i,j,k locally
  IF (.NOT. vel_grid_yes) THEN
    !--- pressure grid indices ----------------------------------------------------------------------
    ! bbecsek: 
    ! if we want to convert the pressure grid indices into nodal IDs we have to
    ! consider that 
    ! if  BC_1U <= 0:  i = [1, N1-1]
    ! if  BC_1U >  0:  i = [1, N1  ]
    ! if  BC_2U <= 0:  j = [1, N2-1]
    ! if  BC_2U >  0:  j = [1, N2  ]
    ! if  BC_3U <= 0:  k = [1, N3-1]
    ! if  BC_3U >  0:  k = [1, N3  ]
    ! nodal IDs run from 0 to nn1xnn2xnn3-1
    ii  = i - 1
    jj  = j - 1
    kk  = k - 1

  ELSE
    !--- velocity grid indices ----------------------------------------------------------------------
    ! bbecsek: 
    ! if we want to convert the velocity grid indices into nodal IDs we have to
    ! consider that
    ! if we exclude boundary points, i.e. use sii and nii
    !   i = [1, N1-1]
    !   j = [1, N2-1]
    !   k = [1, N3-1]
    ! if we include boundary points, i.e. use siib and niib
    ! if  BC_1L <=0 and BC_1U <= 0:  i = [1, N1-1]
    ! if  BC_1L <=0 and BC_1U >  0:  i = [1, N1  ]
    ! if  BC_1L > 0 and BC_1U <= 0:  i = [0, N1-1]
    ! if  BC_1L > 0 and BC_1U >  0:  i = [0, N1  ]
    ! ... and analogously for the other spacial directions
    ! nodal IDs run from 0 to nn1xnn2xnn3-1
    IF ((BC_1L .LE. 0 .AND. dir .EQ. 1) .OR. (.NOT. boundary_yes .AND. dir .EQ. 1) .OR. (dir .NE. 1)) THEN
      ii = i - 1
    ELSE
      ii = i
    END IF

    IF ((BC_2L .LE. 0 .AND. dir .EQ. 2) .OR. (.NOT. boundary_yes .AND. dir .EQ. 2) .OR. (dir .NE. 2)) THEN
      jj = j - 1
    ELSE
      jj = j
    END IF

    IF ((BC_3L .LE. 0 .AND. dir .EQ. 3) .OR. (.NOT. boundary_yes .AND. dir .EQ. 3) .OR. (dir .NE. 3)) THEN
      kk = k - 1
    ELSE
      kk = k
    END IF

  END IF

  !--- This corrects the local indices to become global ones
  !--- 1-direction ---
  IF (iB(1,1) .GT. 1) THEN
    ii = ii + nn1l + (iB(1,1) - 2)*(N1 - 1)
  END IF
  !--- 2-direction ---
  IF (iB(2,1) .GT. 1) THEN
    jj = jj + nn2l + (iB(2,1) - 2)*(N2 - 1)
  END IF
  !--- 3-direction ---
  IF (iB(3,1) .GT. 1) THEN
    kk = kk + nn3l + (iB(3,1) - 2)*(N3 - 1)
  END IF

  node_id = ii + jj * nn1 + kk * nn1 * nn2

  END SUBROUTINE global_cart2id_glob_con


  SUBROUTINE global_id2cart_glob_con(node_id, i, j, k, vel_grid_yes, boundary_yes, dir)

  IMPLICIT NONE

  INTEGER, INTENT(IN)     ::  node_id
  INTEGER                 ::  ijn1
  INTEGER                 ::  nn1,nn2,nn3
  INTEGER                 ::  nn1l, nn1u, nn2l, nn2u, nn3l, nn3u
  INTEGER                 ::  dummy1, dummy2, dummy3
  LOGICAL, INTENT(IN)     ::  vel_grid_yes, boundary_yes
  INTEGER, INTENT(IN)     ::  dir

  INTEGER, INTENT(OUT)    ::  i,j,k

  !--- Get the global number of nodes on which the values are stored
  CALL get_boundary_block_dims(nn1l, nn1u, nn2l, nn2u, nn3l, nn3u, vel_grid_yes, boundary_yes, dir)
  nn1 = nn1l + (NB1 - 2)*(N1 - 1) + nn1u
  nn2 = nn2l + (NB2 - 2)*(N2 - 1) + nn2u
  nn3 = nn3l + (NB3 - 2)*(N3 - 1) + nn3u

  !--- If in any direction only one block exists override dimension with local one
  IF (NB1 .EQ. 1) THEN
    CALL get_block_dims(nn1, dummy1, dummy2, vel_grid_yes, boundary_yes, dir)
  ENDIF
  IF (NB2 .EQ. 1) THEN
    CALL get_block_dims(dummy1, nn2, dummy2, vel_grid_yes, boundary_yes, dir)
  ENDIF
  IF (NB3 .EQ. 1) THEN
    CALL get_block_dims(dummy1, dummy2, nn3, vel_grid_yes, boundary_yes, dir)
  END IF

  !IF (node_id .LT. nn1) THEN
  !  i = node_id
  !  j = 0
  !  k = 0
  !ELSE IF (node_id .GE. nn1 .AND. node_id .LT. nn1*nn2) THEN
  !  i = MOD(node_id, nn1)
  !  j = (node_id - i)/nn1
  !  k = 0
  !ELSE IF (node_id .GE. nn1*nn2) THEN
  !  ijn1 = MOD(node_id, nn1*nn2)
  !  i = MOD(ijn1, nn1)
  !  j = (ijn1 - i)/nn1
  !  k = (node_id - i - j*nn1)/(nn1*nn2)
  !END IF
  CALL id2cart(node_id, i, j, k, nn1, nn2, nn3)

  !--- This corrects the local indices to become global ones
  !--- 1-direction ---
  IF (iB(1,1) .GT. 1) THEN
    i = i - nn1l - (iB(1,1) - 2)*(N1 - 1)
  END IF
  !--- 2-direction ---
  IF (iB(2,1) .GT. 1) THEN
    j = j - nn2l - (iB(2,1) - 2)*(N2 - 1)
  END IF
  !--- 3-direction ---
  IF (iB(3,1) .GT. 1) THEN
    k = k - nn3l - (iB(3,1) - 2)*(N3 - 1)
  END IF

  !--- correct indices i,j,k locally
  IF (.NOT. vel_grid_yes) THEN
    i = i + 1
    j = j + 1
    k = k + 1
  ELSE
    IF ((BC_1L .LE. 0 .AND. dir .EQ. 1) .OR. (.NOT. boundary_yes .AND. dir .EQ. 1) .OR. (dir .NE. 1)) THEN
      i = i + 1
    END IF

    IF ((BC_2L .LE. 0 .AND. dir .EQ. 2) .OR. (.NOT. boundary_yes .AND. dir .EQ. 2) .OR. (dir .NE. 2)) THEN
      j = j + 1
    END IF

    IF ((BC_3L .LE. 0 .AND. dir .EQ. 3) .OR. (.NOT. boundary_yes .AND. dir .EQ. 3) .OR. (dir .NE. 3)) THEN
      k = k + 1
    END IF

  END IF

  END SUBROUTINE global_id2cart_glob_con


  ! Cartesian to index translation, only works for the pressure grid at the
  ! moment with periodic BCs
  SUBROUTINE global_cart2id_loc_con(node_id, i, j, k)

  IMPLICIT NONE

  INTEGER, INTENT(IN)     ::  i,j,k
  INTEGER                 ::  ii,jj,kk
  INTEGER                 ::  nn1,nn2,nn3
  INTEGER                 ::  nn_block, blockid
  INTEGER, INTENT(OUT)    ::  node_id

  CALL get_block_dims(nn1, nn2, nn3, .FALSE., .FALSE., -1)

  nn_block = nn1 * nn2 * nn3

  !--- Pressure grid indices
  ii = i - 1
  jj = j - 1
  kk = k - 1

  !--- global node id in this current block w/o consideration of other blocks
  node_id = ii + jj * nn1 + kk * nn1 * nn2

  !--- get the current block's id
  blockid = block_id(iB(1,1), iB(2,1), iB(3,1))

  !--- augment the node id to consider its location in the entire domain
  node_id = node_id + SUM( block_sizes(0:(blockid-1)) )

  END SUBROUTINE global_cart2id_loc_con




  SUBROUTINE global_id2cart_loc_con(node_id, i, j, k)

  IMPLICIT NONE

  INTEGER, INTENT(OUT)     ::  i,j,k
  INTEGER                  ::  ii,jj,kk
  INTEGER                  ::  nn1,nn2,nn3
  INTEGER                  ::  nn_block, blockid, nid
  INTEGER, INTENT(IN)      ::  node_id

  CALL get_block_dims(nn1, nn2, nn3, .FALSE., .FALSE., -1)

  nn_block = nn1 * nn2 * nn3

  !--- get the current block's id
  blockid = block_id(iB(1,1), iB(2,1), iB(3,1))

  !--- reduce the node id to be within this blocks ids
  nid = node_id - SUM( block_sizes(0:(blockid-1)) )

  !--- regain indices
  CALL id2cart(nid, ii, jj, kk, nn1, nn2, nn3)

  !--- pressure grid indices
  i = ii + 1
  j = jj + 1
  k = kk + 1



  END SUBROUTINE global_id2cart_loc_con




  FUNCTION n_local_nodes() RESULT (fn_val)

  IMPLICIT NONE
  INTEGER                :: nn1, nn2, nn3
  INTEGER                :: fn_val

  CALL get_block_dims(nn1, nn2, nn3, .FALSE., .FALSE., -1)

  fn_val =  nn1 * nn2 * nn3


  END FUNCTION n_local_nodes

  FUNCTION n_global_nodes() RESULT(fn_val)

  IMPLICIT NONE
  INTEGER                :: nn1, nn2, nn3
  INTEGER                :: fn_val

  fn_val = SUM(block_sizes)

  END FUNCTION n_global_nodes


  !> local node IDs with one row of ghost cells included
  SUBROUTINE local_cart2id_loc_con(node_id,i,j,k,n_local,n_ghost)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: i, j, k
  INTEGER, INTENT(OUT)   :: node_id
  INTEGER, INTENT(OUT)   :: n_local, n_ghost
  INTEGER                :: nn1, nn2, nn3
  INTEGER                :: ii, jj, kk

  CALL get_block_dims(nn1,nn2,nn3,.FALSE.,.FALSE.,-1)
  n_local = nn1*nn2*nn3

  !--- correct for rows of ghost vertices (one on each side) ---
  !--- 1-direction
  IF (iB(1,1) .EQ. 1 .AND. BC_1L_global .LT. 0) THEN
    nn1 = nn1
    ii  = i - 1
  ELSE IF (BC_1L .EQ. 0) THEN
    nn1 = nn1 + 1
    ii  = i
  ELSE
    ii  = i - 1
  END IF
  !--- 2-direction
  IF (iB(2,1) .EQ. 1 .AND. BC_2L_global .LT. 0) THEN
    nn2 = nn2
    jj  = j - 1
  ELSE IF (BC_2L .EQ. 0) THEN
    nn2 = nn2 + 1
    jj  = j
  ELSE
    jj  = j - 1
  END IF
  !--- 3-direction
  IF (dimens .EQ. 3) THEN
    IF (iB(3,1) .EQ. 1 .AND. BC_3L_global .LT. 0) THEN
      nn3 = nn3
      kk  = k - 1
    ELSE IF (BC_3L .EQ. 0) THEN
      nn3 = nn3 + 1
      kk  = k
    ELSE
      kk  = k - 1
    END IF
  ELSE
    nn3 = nn3
    kk  = 0
  END IF

  !--- calculate number of ghost nodes ---
  n_ghost = nn1*nn2*nn3 - n_local

  !--- compute node id ---
  node_id = ii + jj*nn1 + kk*nn1*nn2


  END SUBROUTINE local_cart2id_loc_con


  SUBROUTINE local_id2cart_loc_con(node_id,i,j,k,n_local,n_ghost)

  IMPLICIT NONE

  INTEGER, INTENT(OUT)   :: i, j, k
  INTEGER, INTENT(IN)    :: node_id
  INTEGER, INTENT(OUT)   :: n_local, n_ghost
  INTEGER                :: nn1, nn2, nn3
  INTEGER                :: ijn1

  CALL get_block_dims(nn1,nn2,nn3,.FALSE.,.FALSE.,-1)
  n_local = nn1*nn2*nn3

  !--- correct for rows of ghost vertices (one on each side) ---
  !--- 1-direction
  IF (iB(1,1) .EQ. 1 .AND. BC_1L_global .LT. 0) THEN
    nn1 = nn1
  ELSE IF (BC_1L .EQ. 0) THEN
    nn1 = nn1 + 1
  END IF
  !--- 2-direction
  IF (iB(2,1) .EQ. 1 .AND. BC_2L_global .LT. 0) THEN
    nn2 = nn2
  ELSE IF (BC_2L .LE. 0) THEN
    nn2 = nn2 + 1
  END IF
  !--- 3-direction
  IF (dimens .EQ. 3) THEN
    IF (iB(3,1) .EQ. 1 .AND. BC_3L_global .LT. 0) THEN
      nn3 = nn3
    ELSE IF (BC_3L .LE. 0) THEN
      nn3 = nn3 + 1
    END IF
  ELSE
    nn3 = nn3
  END IF

  !--- calculate number of ghost nodes ---
  n_ghost = nn1*nn2*nn3 - n_local

  !--- convert back ---
  !IF (node_id .LT. nn1) THEN
  !  i = node_id
  !  j = 0
  !  k = 0
  !ELSE IF (node_id .GE. nn1 .AND. node_id .LT. nn1*nn2) THEN
  !  i = MOD(node_id, nn1)
  !  j = (node_id - i)/nn1
  !  k = 0
  !ELSE IF (node_id .GE. nn1*nn2) THEN
  !  ijn1 = MOD(node_id, nn1*nn2)
  !  i = MOD(ijn1, nn1)
  !  j = (ijn1 - i)/nn1
  !  k = (node_id - i - j*nn1)/(nn1*nn2)
  !END IF
  CALL id2cart(node_id, i, j, k, nn1, nn2, nn3)

  !--- correct indices for rows of ghostcells ---
  IF (BC_1L .GT. 0 .OR. (iB(1,1) .EQ. 1 .AND. BC_1L_global .LT. 0)) THEN
    i = i + 1
  END IF
  IF (BC_2L .GT. 0 .OR. (iB(2,1) .EQ. 1 .AND. BC_2L_global .LT. 0)) THEN
    j = j + 1
  END IF
  IF (BC_3L .GT. 0 .AND. dimens .EQ.  3 .OR. (iB(3,1) .EQ. 1 .AND. BC_3L_global .LT. 0 .AND. dimens .EQ. 3)) THEN
    k = k + 1
  END IF

  IF (dimens .NE. 3) k = 1

  END SUBROUTINE local_id2cart_loc_con




  SUBROUTINE local_elem_base_loc_con(node_id,i,j,k)

  IMPLICIT NONE

  INTEGER, INTENT(OUT)   :: i, j, k
  INTEGER, INTENT(IN)    :: node_id
  INTEGER                :: nn1, nn2, nn3
  INTEGER                :: ijn1

  CALL get_block_dims(nn1,nn2,nn3,.FALSE.,.FALSE.,-1)

  !--- correct for 1 rows of ghost vertices (one on each side) ---
  ! per default we want one less node per direction as the base i
  ! (POTENTIALLY DEPRECATED)
  ! nn1 = nn1 - 1
  ! nn2 = nn2 - 1
  ! nn3 = nn3 - 1
  ! per default we want two less nodes per direction to avoid duplicate
  ! elements across process interfaces
  nn1 = nn1 - 1
  nn2 = nn2 - 1
  nn3 = nn3 - 1
  IF (iB(1,1) .EQ. 1 .AND. BC_1L_global .LT. 0) THEN
    nn1 = nn1
  ELSE IF (BC_1L .LE. 0) THEN
    nn1 = nn1 + 1
  END IF

  IF (iB(2,1) .EQ. 1 .AND. BC_2L_global .LT. 0) THEN
    nn2 = nn2
  ELSE IF (BC_2L .LE. 0) THEN
    nn2 = nn2 + 1
  END IF

  IF (dimens .EQ. 3) THEN
    IF (iB(3,1) .EQ. 1 .AND. BC_3L_global .LT. 0) THEN
      nn3 = nn3
    ELSE IF (BC_3L .LE. 0) THEN
      nn3 = nn3 + 1
    END IF
  END IF


  !--- convert back ---
  !IF (node_id .LT. nn1) THEN
  !  i = node_id
  !  j = 0
  !  k = 0
  !ELSE IF (node_id .GE. nn1 .AND. node_id .LT. nn1*nn2) THEN
  !  i = MOD(node_id, nn1)
  !  j = (node_id - i)/nn1
  !  k = 0
  !ELSE IF (node_id .GE. nn1*nn2) THEN
  !  ijn1 = MOD(node_id, nn1*nn2)
  !  i = MOD(ijn1, nn1)
  !  j = (ijn1 - i)/nn1
  !  k = (node_id - i - j*nn1)/(nn1*nn2)
  !END IF
  CALL id2cart(node_id, i, j, k, nn1, nn2, nn3)

  !--- correct indices for rows of ghostcells ---
  IF (BC_1L .GT. 0) THEN
    i = i + 1
  ELSE IF (iB(1,1) .EQ. 1 .AND. BC_1L_global .LT. 0) THEN
    i = i + 1
  END IF
  IF (BC_2L .GT. 0) THEN
    j = j + 1
  ELSE IF (iB(2,1) .EQ. 1 .AND. BC_2L_global .LT. 0) THEN
    j = j + 1
  END IF
  IF (BC_3L .GT. 0 .AND. dimens .EQ. 3) THEN
    k = k + 1
  ELSE IF (iB(3,1) .EQ. 1 .AND. BC_3L_global .LT. 0 .AND. dimens .EQ. 3) THEN
    k = k + 1
  END IF

  IF (dimens .NE. 3) k = 1

  END SUBROUTINE local_elem_base_loc_con



  SUBROUTINE id2cart(id, ii, jj, kk, nn1, nn2, nn3)

  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: id
  INTEGER, INTENT(OUT)  :: ii, jj, kk
  INTEGER, INTENT(IN)   :: nn1, nn2, nn3
  INTEGER               :: ijn1


  IF (id .LT. nn1) THEN
    ii = id
    jj = 0
    kk = 0
  ELSE IF (id .GE. nn1 .AND. id .LT. nn1*nn2) THEN
    ii = MOD(id, nn1)
    jj = (id - ii)/nn1
    kk = 0
  ELSE IF (id .GE. nn1*nn2) THEN
    ijn1 = MOD(id, nn1*nn2)
    ii = MOD(ijn1, nn1)
    jj = (ijn1 - ii)/nn1
    kk = (id - ii - jj*nn1)/(nn1*nn2)
  END IF

  END SUBROUTINE id2cart


  !> returns the local node_id to the node_no-th node of elem_no element
  SUBROUTINE local_tet_element_nodes(elem_no, node_no, node_id)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: elem_no, node_no
  INTEGER, INTENT(OUT)   :: node_id
  INTEGER                :: base_node_id, tet_type
  INTEGER                :: ii,jj,kk
  INTEGER                :: n_local, n_ghost
  INTEGER                :: n_elem

  !--- check whether the entered elem_no is actually within this process ---
  n_elem = total_n_local_tet_elements()
  
  IF (elem_no .GE. n_elem) THEN
      WRITE(*,*)"WARNING: requested tetrahedral element out of range on process ",rank,". Aborting ..."
      CALL MPI_FINALIZE(merror)
      STOP
  END IF

  !--- compute the type of tetrahedron inside the hexahedron ---
  tet_type = MOD(elem_no,6)
  
  !--- now that we know the tet-type in [0,5], compute the base node ID ---
  base_node_id = (elem_no - tet_type)/6
  
  !--- convert the local base_node_id to local i,j,k ---
  CALL local_elem_base_loc_con(base_node_id,ii,jj,kk)

  !--- based on the tet_type compute the node_no node's node_id ---
  SELECT CASE (tet_type)
    CASE (0)
      SELECT CASE (node_no)
        CASE (0)
          ii = ii
          jj = jj
          kk = kk
        CASE (1)
          ii = ii + 1
          jj = jj + 1
          kk = kk + 1
        CASE (2)
          ii = ii
          jj = jj + 1
          kk = kk + 1
        CASE (3)
          ii = ii
          jj = jj
          kk = kk + 1
        CASE DEFAULT
          !--- check whether requested element node is existing ---
          WRITE(*,*)"WARNING: requested tetrahedron node out of range on process",rank,". Aborting ..."
          CALL MPI_FINALIZE(merror)
          STOP
      END SELECT
    CASE (1)
      SELECT CASE (node_no)
        CASE (0)
          ii = ii
          jj = jj
          kk = kk
        CASE (1)
          ii = ii + 1
          jj = jj + 1
          kk = kk + 1
        CASE (2)
          ii = ii
          jj = jj + 1
          kk = kk
        CASE (3)
          ii = ii
          jj = jj + 1
          kk = kk + 1
        CASE DEFAULT
          !--- check whether requested element node is existing ---
          WRITE(*,*)"WARNING: requested tetrahedron node out of range on process",rank,". Aborting ..."
          CALL MPI_FINALIZE(merror)
          STOP
      END SELECT
    CASE (2)
      SELECT CASE (node_no)
        CASE (0)
          ii = ii
          jj = jj
          kk = kk
        CASE (1)
          ii = ii + 1
          jj = jj
          kk = kk + 1
        CASE (2)
          ii = ii + 1
          jj = jj + 1
          kk = kk + 1
        CASE (3)
          ii = ii
          jj = jj
          kk = kk + 1
        CASE DEFAULT
          !--- check whether requested element node is existing ---
          WRITE(*,*)"WARNING: requested tetrahedron node out of range on process",rank,". Aborting ..."
          CALL MPI_FINALIZE(merror)
          STOP
      END SELECT
    CASE (3)
      SELECT CASE (node_no)
        CASE (0)
          ii = ii
          jj = jj
          kk = kk
        CASE (1)
          ii = ii + 1
          jj = jj + 1
          kk = kk
        CASE (2)
          ii = ii
          jj = jj + 1
          kk = kk
        CASE (3)
          ii = ii + 1
          jj = jj + 1
          kk = kk + 1
        CASE DEFAULT
          !--- check whether requested element node is existing ---
          WRITE(*,*)"WARNING: requested tetrahedron node out of range on process",rank,". Aborting ..."
          CALL MPI_FINALIZE(merror)
          STOP
      END SELECT
    CASE (4)
      SELECT CASE (node_no)
        CASE (0)
          ii = ii
          jj = jj
          kk = kk
        CASE (1)
          ii = ii + 1
          jj = jj
          kk = kk
        CASE (2)
          ii = ii + 1
          jj = jj + 1
          kk = kk
        CASE (3)
          ii = ii + 1
          jj = jj + 1
          kk = kk + 1
        CASE DEFAULT
          !--- check whether requested element node is existing ---
          WRITE(*,*)"WARNING: requested tetrahedron node out of range on process",rank,". Aborting ..."
          CALL MPI_FINALIZE(merror)
          STOP
      END SELECT
    CASE (5)
      SELECT CASE (node_no)
        CASE (0)
          ii = ii
          jj = jj
          kk = kk
        CASE (1)
          ii = ii + 1
          jj = jj
          kk = kk
        CASE (2)
          ii = ii + 1
          jj = jj + 1
          kk = kk + 1
        CASE (3)
          ii = ii + 1
          jj = jj
          kk = kk + 1
        CASE DEFAULT
          !--- check whether requested element node is existing ---
          WRITE(*,*)"WARNING: requested tetrahedron node out of range on process",rank,". Aborting ..."
          CALL MPI_FINALIZE(merror)
          STOP
      END SELECT
    CASE DEFAULT
      WRITE(*,*)"WARNING: invalid tetrahedron type encountered. Aborting..."
      CALL MPI_FINALIZE(merror)
      STOP
  END SELECT

  !--- convert new i,j,k back to node_id ---
  CALL local_cart2id_loc_con(node_id, ii, jj, kk, n_local, n_ghost)

  END SUBROUTINE local_tet_element_nodes


  SUBROUTINE local_hex_element_nodes(elem_no, node_no, node_id)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: elem_no, node_no
  INTEGER, INTENT(OUT)   :: node_id
  INTEGER                :: base_node_id
  INTEGER                :: ii,jj,kk
  INTEGER                :: n_local, n_ghost
  INTEGER                :: n_elem

  !--- check whether the entered elem_no is actually within this process ---
  n_elem = total_n_local_hex_elements()
  
  IF (elem_no .GE. n_elem) THEN
      WRITE(*,*)"WARNING: requested hexahedral element out of range on process ",rank,". Aborting ..."
      CALL MPI_FINALIZE(merror)
      STOP
  END IF

  !--- the elements base node is equal to the elem_no since there is one element
  !--- per node
  base_node_id = elem_no

  !--- convert the local base_node_id to local i,j,k ---
  CALL local_elem_base_loc_con(base_node_id,ii,jj,kk)

  !--- compute the node_no node's node_id ---
  SELECT CASE (node_no)
    CASE (0)
     ii = ii
     jj = jj
     kk = kk
    CASE (1)
     ii = ii + 1
     jj = jj
     kk = kk
    CASE (2)
     ii = ii + 1
     jj = jj + 1
     kk = kk
    CASE (3)
     ii = ii
     jj = jj + 1
     kk = kk
    CASE (4)
     ii = ii
     jj = jj
     kk = kk + 1
    CASE (5)
     ii = ii + 1
     jj = jj
     kk = kk + 1
    CASE (6)
     ii = ii + 1
     jj = jj + 1
     kk = kk + 1
    CASE (7)
     ii = ii
     jj = jj + 1
     kk = kk + 1
    CASE DEFAULT
     !--- check whether requested element node is existing ---
     WRITE(*,*)"WARNING: requested hexahedron node out of range on process",rank,". Aborting ..."
     CALL MPI_FINALIZE(merror)
     STOP
  END SELECT

  !--- convert new i,j,k back to node_id ---
  CALL local_cart2id_loc_con(node_id, ii, jj, kk, n_local, n_ghost)

  END SUBROUTINE local_hex_element_nodes



  FUNCTION local_to_global_tet_elem_id(loc_elem_id) RESULT(fn_val)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: loc_elem_id
  INTEGER                :: fn_val

  fn_val = loc_elem_id &
         + SUM(block_elem_sizes( 0:(block_id(iB(1,1), iB(2,1), iB(3,1))-1) ))

  END FUNCTION local_to_global_tet_elem_id

  FUNCTION local_to_global_hex_elem_id(loc_elem_id) RESULT(fn_val)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: loc_elem_id
  INTEGER                :: fn_val

  fn_val = loc_elem_id &
         + SUM(block_elem_sizes( 0:(block_id(iB(1,1), iB(2,1), iB(3,1))-1) ))

  END FUNCTION local_to_global_hex_elem_id

  FUNCTION local_to_global_tri_elem_id(loc_elem_id) RESULT(fn_val)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: loc_elem_id
  INTEGER                :: fn_val

  fn_val = loc_elem_id &
         + SUM(block_elem_sizes( 0:(block_id(iB(1,1), iB(2,1), iB(3,1))-1) ))

  END FUNCTION local_to_global_tri_elem_id

  FUNCTION local_to_global_quad_elem_id(loc_elem_id) RESULT(fn_val)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: loc_elem_id
  INTEGER                :: fn_val

  fn_val = loc_elem_id &
         + SUM(block_elem_sizes( 0:(block_id(iB(1,1), iB(2,1), iB(3,1))-1) ))

  END FUNCTION local_to_global_quad_elem_id


  SUBROUTINE local_tri_element_nodes(elem_no, node_no, node_id)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: elem_no, node_no
  INTEGER, INTENT(OUT)   :: node_id
  INTEGER                :: base_node_id, tri_type
  INTEGER                :: ii,jj,kk
  INTEGER                :: n_elem
  INTEGER                :: n_local, n_ghost

  !--- check whether the entered elem_no is actually within this process ---
  n_elem = total_n_local_tri_elements()

  IF (elem_no .GE. n_elem) THEN
    WRITE(*,*)"WARNING: requested triangular element out of range on process ",rank,". Aborting ..."
    CALL MPI_FINALIZE(merror)
    STOP
  END IF

  !--- compute the type of triangle inside the quaadrilateral ---
  tri_type = MOD(elem_no,2)

  !--- now that we know the tri-type in [0,1], compute the base node ID ---
  base_node_id = (elem_no - tri_type)/2

  !--- convert the local base_node_id to local i,j,k ---
  CALL local_elem_base_loc_con(base_node_id,ii,jj,kk)

  !--- based on the tri_type compute the node_no node's node_id ---
  SELECT CASE (tri_type)
    CASE (0)
      SELECT CASE (node_no)
        CASE (0)
          ii = ii
          jj = jj
        CASE (1)
          ii = ii + 1
          jj = jj + 1
        CASE (2)
          ii = ii
          jj = jj + 1
        CASE DEFAULT
          !--- check whether requested element node is existing ---
          WRITE(*,*)"WARNING: requested triangle node out of range on process",rank,". Aborting ..."
          CALL MPI_FINALIZE(merror)
          STOP
      END SELECT
    CASE (1)
      SELECT CASE (node_no)
        CASE (0)
          ii = ii
          jj = jj
        CASE (1)
          ii = ii + 1
          jj = jj
        CASE (2)
          ii = ii + 1
          jj = jj + 1
        CASE DEFAULT
          !--- check whether requested element node is existing ---
          WRITE(*,*)"WARNING: requested triangle node out of range on process",rank,". Aborting ..."
          CALL MPI_FINALIZE(merror)
          STOP
      END SELECT
    CASE DEFAULT
      WRITE(*,*)"WARNING: invalid triangle type encountered. Aborting..."
      CALL MPI_FINALIZE(merror)
      STOP
  END SELECT

  !--- convert new i,j,k back to node_id ---
  CALL local_cart2id_loc_con(node_id, ii, jj, kk, n_local, n_ghost)


  END SUBROUTINE local_tri_element_nodes



  SUBROUTINE local_quad_element_nodes(elem_no, node_no, node_id)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: elem_no, node_no
  INTEGER, INTENT(OUT)   :: node_id
  INTEGER                :: base_node_id
  INTEGER                :: ii,jj,kk
  INTEGER                :: n_elem
  INTEGER                :: n_local, n_ghost

  !--- check whether the entered elem_no is actually within this process ---
  n_elem = total_n_local_quad_elements()

  IF (elem_no .GE. n_elem) THEN
    WRITE(*,*)"WARNING: requested quadrilateral element out of range on process ",rank,". Aborting ..."
    CALL MPI_FINALIZE(merror)
    STOP
  END IF
 
  !--- compute the base node ID, in this case it is equal to the element number
  !--- because there is exactly one element per node, just for completeness
  base_node_id = elem_no

  !--- convert the local base_node_id to local i,j,k ---
  CALL local_elem_base_loc_con(base_node_id,ii,jj,kk)

  !--- compute the node_no node's node_id ---
  SELECT CASE (node_no)
    CASE (0)
      ii = ii
      jj = jj 
    CASE (1)
      ii = ii + 1
      jj = jj
    CASE (2)
      ii = ii + 1
      jj = jj + 1
    CASE (3)
      ii = ii
      jj = jj + 1
    CASE DEFAULT
      !--- check whether requested element node is existing ---
      WRITE(*,*)"WARNING: requested quadrilateral node out of range on process",rank,". Aborting ..."
      CALL MPI_FINALIZE(merror)
      STOP
  END SELECT
  
  !--- convert new i,j,k back to node_id ---
  CALL local_cart2id_loc_con(node_id, ii, jj, kk, n_local, n_ghost)

  END SUBROUTINE local_quad_element_nodes



  FUNCTION total_n_local_tri_elements() RESULT(fn_val)

  IMPLICIT NONE

  INTEGER                :: fn_val
  INTEGER                :: nn1, nn2, nn3

  IF (dimens .EQ. 3) THEN
    WRITE(*,*)"WARNING: Tri elements only supported for 2D simulations. Aborting ..."
    CALL MPI_FINALIZE(merror)
    STOP
  END IF

  CALL get_block_dims(nn1,nn2,nn3,.FALSE.,.FALSE.,-1)

  !--- correct for rows of ghost vertices (one on each side) ---
  !--- 1-direction ---
  IF (iB(1,1) .EQ. 1 .AND. BC_1L_global .LT. 0) THEN
    nn1 = nn1
  ELSE IF (BC_1L .LE. 0) THEN
    nn1 = nn1 + 1
  END IF
  !--- 2-direction ---
  IF (iB(2,1) .EQ. 1 .AND. BC_2L_global .LT. 0) THEN
    nn2 = nn2
  ELSE IF (BC_2L .LE. 0) THEN
    nn2 = nn2 + 1
  END IF

  !--- compute number of elements ---
  ! minus 2 to avoid duplicate elements across process boundaries
  fn_val = (nn1-1)*(nn2-1)*2

  END FUNCTION total_n_local_tri_elements



  FUNCTION total_n_local_quad_elements() RESULT(fn_val)

  IMPLICIT NONE

  INTEGER                :: fn_val
  INTEGER                :: nn1, nn2, nn3

  IF (dimens .EQ. 3) THEN
    WRITE(*,*)"WARNING: Quad elements only supported for 2D simulations. Aborting ..."
    CALL MPI_FINALIZE(merror)
    STOP
  END IF

  CALL get_block_dims(nn1,nn2,nn3,.FALSE.,.FALSE.,-1)

  !--- correct for rows of ghost vertices (one on each side) ---
  !--- 1-direction ---
  IF (iB(1,1) .EQ. 1 .AND. BC_1L_global .LT. 0) THEN
    nn1 = nn1
  ELSE IF (BC_1L .LE. 0) THEN
    nn1 = nn1 + 1
  END IF
  !--- 2-direction ---
  IF (iB(2,1) .EQ. 1 .AND. BC_2L_global .LT. 0) THEN
    nn2 = nn2
  ELSE IF (BC_2L .LE. 0) THEN
    nn2 = nn2 + 1
  END IF

  !--- compute number of elements ---
  ! minus 2 to avoid duplicate elements across process boundaries
  fn_val = (nn1-1)*(nn2-1)

  END FUNCTION total_n_local_quad_elements




  FUNCTION total_n_local_tet_elements() RESULT(fn_val)

  IMPLICIT NONE

  INTEGER                :: fn_val
  INTEGER                :: nn1, nn2, nn3

  IF (dimens .EQ. 2) THEN
    WRITE(*,*)"WARNING: Tet elements only supported for 3D simulations. Aborting ..."
    CALL MPI_FINALIZE(merror)
    STOP
  END IF


  CALL get_block_dims(nn1,nn2,nn3,.FALSE.,.FALSE.,-1)

  !--- correct for rows of ghost vertices (one on each side) ---
  !--- 1-direction
  IF (iB(1,1) .EQ. 1 .AND. BC_1L_global .LT. 0) THEN
    nn1 = nn1
  ELSE IF (BC_1L .LE. 0) THEN
    nn1 = nn1 + 1
  END IF
  !--- 2-direction
  IF (iB(2,1) .EQ. 1 .AND. BC_2L_global .LT. 0) THEN
    nn2 = nn2
  ELSE IF (BC_2L .LE. 0) THEN
    nn2 = nn2 + 1
  END IF
  !--- 3-direction
  IF (iB(3,1) .EQ. 1 .AND. BC_3L_global .LT. 0) THEN
    nn3 = nn3
  ELSE IF (BC_3L .LE. 0) THEN
    nn3 = nn3 + 1
  END IF

  !--- compute number of elements ---
  ! minus 2 to avoid duplicate elements across process boundaries
  fn_val = (nn1-1)*(nn2-1)*(nn3-1)*6


  END FUNCTION total_n_local_tet_elements




  FUNCTION total_n_local_hex_elements() RESULT(fn_val)

  IMPLICIT NONE

  INTEGER                :: fn_val
  INTEGER                :: nn1, nn2, nn3

  IF (dimens .EQ. 2) THEN
    WRITE(*,*)"WARNING: Hex elements only supported for 3D simulations. Aborting ..."
    CALL MPI_FINALIZE(merror)
    STOP
  END IF


  CALL get_block_dims(nn1,nn2,nn3,.FALSE.,.FALSE.,-1)

  !--- correct for rows of ghost vertices (one on each side) ---
  !--- 1-direction
  IF (iB(1,1) .EQ. 1 .AND. BC_1L_global .LT. 0) THEN
    nn1 = nn1
  ELSE IF (BC_1L .LE. 0) THEN
    nn1 = nn1 + 1
  END IF
  !--- 2-direction
  IF (iB(2,1) .EQ. 1 .AND. BC_2L_global .LT. 0) THEN
    nn2 = nn2
  ELSE IF (BC_2L .LE. 0) THEN
    nn2 = nn2 + 1
  END IF
  !--- 3-direction
  IF (iB(3,1) .EQ. 1 .AND. BC_3L_global .LT. 0) THEN
    nn3 = nn3
  ELSE IF (BC_3L .LE. 0) THEN
    nn3 = nn3 + 1
  END IF

  !--- compute number of elements ---
  ! minus 2 to avoid duplicate elements across process boundaries
  fn_val = (nn1-1)*(nn2-1)*(nn3-1)


  END FUNCTION total_n_local_hex_elements





  FUNCTION total_n_global_tri_elements() RESULT(fn_val)

  IMPLICIT NONE

  INTEGER                :: fn_val

  fn_val = SUM(block_elem_sizes)

  END FUNCTION total_n_global_tri_elements

  FUNCTION total_n_global_quad_elements() RESULT(fn_val)

  IMPLICIT NONE

  INTEGER                :: fn_val

  fn_val = SUM(block_elem_sizes)

  END FUNCTION total_n_global_quad_elements

  FUNCTION total_n_global_tet_elements() RESULT(fn_val)

  IMPLICIT NONE

  INTEGER                :: fn_val

  fn_val = SUM(block_elem_sizes)

  END FUNCTION total_n_global_tet_elements

  FUNCTION total_n_global_hex_elements() RESULT(fn_val)

  IMPLICIT NONE

  INTEGER                :: fn_val

  fn_val = SUM(block_elem_sizes)

  END FUNCTION total_n_global_hex_elements


  !> New and updated, much slicker routine that makes use of the start and end
  !! indices
  SUBROUTINE get_block_dims(nn1, nn2, nn3, vel_grid_yes, boundary_yes, dir)

  IMPLICIT NONE

  INTEGER, INTENT(OUT)   ::  nn1, nn2, nn3
  LOGICAL, INTENT(IN )   ::  vel_grid_yes, boundary_yes
  INTEGER, INTENT(IN )   ::  dir
  INTEGER                ::  st1, st2, st3, nd1, nd2, nd3

  CALL get_start_and_end_indices(st1, st2, st3, nd1, nd2, nd3, vel_grid_yes, boundary_yes, dir)

  nn1 = nd1 - st1 + 1
  nn2 = nd2 - st2 + 1
  nn3 = nd3 - st3 + 1

  END SUBROUTINE get_block_dims




  SUBROUTINE get_start_and_end_indices(st1, st2, st3, nd1, nd2, nd3, vel_grid_yes, boundary_yes, dir)

  IMPLICIT NONE

  INTEGER, INTENT(OUT)   ::  st1, st2, st3, nd1, nd2, nd3
  LOGICAL, INTENT(IN )   ::  vel_grid_yes, boundary_yes
  INTEGER, INTENT(IN )   ::  dir
  
  IF (.NOT. vel_grid_yes) THEN
    !--- pressure grid indices ----------------------------------------------------------------------
    ! bbecsek: 
    ! if we want to convert the pressure grid indices into nodal IDs we have to
    ! consider that 
    ! if  BC_1U <= 0:  i = [1, N1-1]
    ! if  BC_1U >  0:  i = [1, N1  ]
    ! if  BC_2U <= 0:  j = [1, N2-1]
    ! if  BC_2U >  0:  j = [1, N2  ]
    ! if  BC_3U <= 0:  k = [1, N3-1]
    ! if  BC_3U >  0:  k = [1, N3  ]
    ! nodal IDs run from 0 to nn1xnn2xnn3-1
    st1 = S1p
    st2 = S2p
    st3 = S3p

    nd1 = N1p
    nd2 = N2p
    nd3 = N3p

  ELSE
    !--- velocity grid indices ----------------------------------------------------------------------
    ! bbecsek:
    ! if we want to convert the velocity grid indices into nodal IDs we have to
    ! consider that
    ! if we exclude boundary points, i.e. use sii and nii
    !   i = [1, N1-1]
    !   j = [1, N2-1]
    !   k = [1, N3-1]
    ! if we include boundary points, i.e. use siib and niib
    ! if  BC_1L <=0 and BC_1U <= 0:  i = [1, N1-1]
    ! if  BC_1L <=0 and BC_1U >  0:  i = [1, N1  ]
    ! if  BC_1L > 0 and BC_1U <= 0:  i = [0, N1-1]
    ! if  BC_1L > 0 and BC_1U >  0:  i = [0, N1  ]
    ! ... and analogously for the other spacial directions
    ! nodal IDs run from 0 to nn1xnn2xnn3-1
    IF (.NOT. boundary_yes) THEN
      IF (dir .EQ. 1) THEN
        st1 = S11
        st2 = S21 
        st3 = S31

        nd1 = N11
        nd2 = N21
        nd3 = N31
      ELSE IF (dir .EQ. 2) THEN
        st1 = S12
        st2 = S22
        st3 = S32

        nd1 = N12
        nd2 = N22
        nd3 = N32
      ELSE IF (dir .EQ. 3) THEN
        st1 = S13
        st2 = S23
        st3 = S33

        nd1 = N13
        nd2 = N23
        nd3 = N33
      ELSE
        WRITE(*,*)'WARNING: Invalid direction specified for get_block_dims. Aborting ...'
        CALL MPI_FINALIZE(merror)
        STOP
      END IF
    ELSE
      IF (dir .EQ. 1) THEN
        st1 = S11B
        st2 = S21B
        st3 = S31B

        nd1 = N11B
        nd2 = N21B
        nd3 = N31B
      ELSE IF (dir .EQ. 2) THEN
        st1 = S12B
        st2 = S22B
        st3 = S32B

        nd1 = N12B
        nd2 = N22B
        nd3 = N32B
      ELSE IF (dir .EQ. 3) THEN
        st1 = S13B
        st2 = S23B
        st3 = S33B

        nd1 = N13B
        nd2 = N23B
        nd3 = N33B
      ELSE
        WRITE(*,*)'WARNING: Invalid direction specified for get_block_dims. Aborting ...'
        CALL MPI_FINALIZE(merror)
        STOP
      END IF
    END IF

  END IF
  
  


  END SUBROUTINE get_start_and_end_indices




  SUBROUTINE get_boundary_block_dims(nn1l, nn1u, nn2l, nn2u, nn3l, nn3u, vel_grid_yes, boundary_yes, dir)

  IMPLICIT NONE

  INTEGER, INTENT(IN )   ::  dir
  LOGICAL, INTENT(IN )   ::  vel_grid_yes, boundary_yes
  
  INTEGER, INTENT(OUT)   ::  nn1l, nn1u, nn2l, nn2u, nn3l, nn3u

  IF (.NOT. vel_grid_yes) THEN
    !--- 1-direction ---
    nn1l = N1 - 1
    IF (BC_1U_global .GT. 0) THEN
      nn1u = N1
    ELSE 
      nn1u = N1 - 1
    END IF
    !--- 2-direction ---
    nn2l = N2 - 1
    IF (BC_2U_global .GT. 0) THEN
      nn2u = N2
    ELSE 
      nn2u = N2 - 1
    END IF
    !--- 3-direction ---
    nn3l = N3 - 1
    IF (BC_3U_global .GT. 0) THEN
      nn3u = N3
    ELSE 
      nn3u = N3 - 1
    END IF
    !-------------------
  ELSE
    !--- 1-direction ---
    !--- lowest block ---
    IF (BC_1L_global .GT. 0 .AND. dir .EQ. 1) THEN
      nn1l = N1
    ELSE
      nn1l = N1 - 1
    END IF
    !--- topmost block ---
    IF ((.NOT. boundary_yes .AND. BC_1U_global .GT. 0 .AND. dir .EQ. 1) .OR. (BC_1U_global .GT. 0 .AND. dir .NE. 1)) THEN
      nn1u = N1
    ELSE
      nn1u = N1 - 1
    END IF
    !--- 2-direction ---
    !--- lowest block ---
    IF (BC_2L_global .GT. 0 .AND. dir .EQ. 2) THEN
      nn2l = N2
    ELSE
      nn2l = N2 - 1
    END IF
    !--- topmost block ---
    IF ((.NOT. boundary_yes .AND. BC_2U_global .GT. 0 .AND. dir .EQ. 2) .OR. (BC_2U_global .GT. 0 .AND. dir .NE. 2)) THEN
      nn2u = N2
    ELSE
      nn2u = N2 - 1
    END IF
    !--- 3-direction ---
    !--- lowest block ---
    IF (BC_3L_global .GT. 0 .AND. dir .EQ. 3) THEN
      nn3l = N3
    ELSE
      nn3l = N3 - 1
    END IF
    !--- topmost block ---
    IF ((.NOT. boundary_yes .AND. BC_3U_global .GT. 0 .AND. dir .EQ. 3) .OR. (BC_3U_global .GT. 0 .AND. dir .NE. 3)) THEN
      nn3u = N3
    ELSE
      nn3u = N3 - 1
    END IF
    !---------------------
  END IF

  END SUBROUTINE get_boundary_block_dims


  !> This routine is called for each node on each process 
  !! i.e. communication may be very expensive!!
  !! maybe a better way is to save each block's # of nodes in a container 
  !! with a routine that is only called once.
  FUNCTION global_node_difference(to_block, from_block) RESULT(fn_val)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: to_block, from_block
  INTEGER                :: fn_val
  INTEGER                :: block
  INTEGER                :: nn1, nn2, nn3
  INTEGER                :: size_recv, size_send
  INTEGER                :: ib1, ib2, ib3

  IF (to_block .GE. from_block) THEN
    fn_val = SIGN(SUM(block_sizes(from_block:to_block-1)), to_block-from_block)
  ELSE
    fn_val = SIGN(SUM(block_sizes(to_block:from_block-1)), to_block-from_block)
  END IF

  END FUNCTION global_node_difference



  SUBROUTINE compute_and_store_global_block_sizes(vel_grid_yes, boundary_yes, dir)

  IMPLICIT NONE

  LOGICAL, INTENT(IN)    ::  vel_grid_yes, boundary_yes
  INTEGER, INTENT(IN)    ::  dir
  INTEGER                ::  nn1, nn2, nn3
  INTEGER                ::  current_elem_size
  INTEGER                ::  current_indices(1:6)
  INTEGER                ::  st1, st2, st3, nd1, nd2, nd3
  INTEGER                ::  i, total_blocks

  block_indices    = 0 ! these variable must be global
  block_sizes      = 0
  block_elem_sizes = 0

  CALL get_start_and_end_indices(st1, st2, st3, nd1, nd2, nd3, vel_grid_yes, boundary_yes, dir)

  current_indices = (/st1, st2, st3, nd1, nd2, nd3/)

  IF (dimens .LT. 3) THEN
    IF (elem_type .EQ. 1) current_elem_size = total_n_local_tri_elements()
    IF (elem_type .EQ. 2) current_elem_size = total_n_local_quad_elements()
  ELSE
    IF (elem_type .EQ. 1) current_elem_size = total_n_local_tet_elements()
    IF (elem_type .EQ. 2) current_elem_size = total_n_local_hex_elements()
  END IF

  CALL MPI_ALLGATHER(current_indices  , 6, MPI_INTEGER, block_indices(0)   , 6, MPI_INTEGER, COMM_CART, merror)
  CALL MPI_ALLGATHER(current_elem_size, 1, MPI_INTEGER, block_elem_sizes(0), 1, MPI_INTEGER, COMM_CART, merror)

  DO i = 0, NB1*NB2*NB3-1
    block_sizes(i) = (block_indices(i*6 + 3) - block_indices(i*6)     + 1) &
                   * (block_indices(i*6 + 4) - block_indices(i*6 + 1) + 1) &
                   * (block_indices(i*6 + 5) - block_indices(i*6 + 2) + 1)
  END DO

 
  END SUBROUTINE compute_and_store_global_block_sizes






  SUBROUTINE interpolate_force_pre_vel(dir)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: dir
  REAL                   :: inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  inter = 0.

  IF (dir .EQ. 1) THEN
    ! exchange across boundaries on pressure grid
    ! also updates ghost cells
    CALL exchange(1,0,fd(b1L,b2L,b3L,1))
    CALL exchange(2,0,fd(b1L,b2L,b3L,1))
    CALL exchange(3,0,fd(b1L,b2L,b3L,1))
    ! Force comp 1 onto vel 1 grid
    CALL interpolate2_pre_vel(.FALSE., 1, fd(b1L,b2L,b3L,1), inter)
    fd(:,:,:,1) = inter(:,:,:)
    ! exchange across boundaries on vel grid
    CALL exchange(1,1,fd(b1L,b2L,b3L,1))
    CALL exchange(2,1,fd(b1L,b2L,b3L,1))
    CALL exchange(3,1,fd(b1L,b2L,b3L,1))
  END IF

  IF (dir .EQ. 2) THEN
    ! exchange across boundaries on pressure grid
    ! also updates ghost cells
    CALL exchange(1,0,fd(b1L,b2L,b3L,2))
    CALL exchange(2,0,fd(b1L,b2L,b3L,2))
    CALL exchange(3,0,fd(b1L,b2L,b3L,2))
    ! Force comp 2 onto vel 2 grid
    CALL interpolate2_pre_vel(.FALSE., 2, fd(b1L,b2L,b3L,2), inter)
    fd(:,:,:,2) = inter(:,:,:)
    ! exchange across boundaries on vel grid
    CALL exchange(1,2,fd(b1L,b2L,b3L,2))
    CALL exchange(2,2,fd(b1L,b2L,b3L,2))
    CALL exchange(3,2,fd(b1L,b2L,b3L,2))
  END IF

  IF (dir .EQ. 3) THEN
    ! exchange across boundaries on pressure grid
    ! also updates ghost cells
    CALL exchange(1,0,fd(b1L,b2L,b3L,3))
    CALL exchange(2,0,fd(b1L,b2L,b3L,3))
    CALL exchange(3,0,fd(b1L,b2L,b3L,3))
    ! Force comp 3 onto vel 3 grid
    CALL interpolate2_pre_vel(.FALSE., 3, fd(b1L,b2L,b3L,3), inter)
    fd(:,:,:,3) = inter(:,:,:)
    ! exchange across boundaries on vel grid
    CALL exchange(1,3,fd(b1L,b2L,b3L,3))
    CALL exchange(2,3,fd(b1L,b2L,b3L,3))
    CALL exchange(3,3,fd(b1L,b2L,b3L,3))
  END IF

  END SUBROUTINE interpolate_force_pre_vel



  SUBROUTINE residual2volume_force(dir)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: dir
  INTEGER                :: i,j,k
  REAL                   :: vol

  ! At the moment it is all done on the pressure grid, so before interpolating
  ! the force onto the velocity grids!

  IF (dimens .EQ. 3) THEN
    DO k = S3p, N3p
      DO j = S2p, N2p
        DO i = S1p, N1p
          IF (elem_type .EQ. 1) THEN ! TET4
            vol = (dx1p(i  ) * dx2p(j  ) * dx3p(k  ))/ 4. &
                + (dx1p(i+1) * dx2p(j  ) * dx3p(k  ))/12. &
                + (dx1p(i  ) * dx2p(j+1) * dx3p(k  ))/12. &
                + (dx1p(i  ) * dx2p(j  ) * dx3p(k+1))/12. &
                + (dx1p(i+1) * dx2p(j+1) * dx3p(k  ))/12. &
                + (dx1p(i+1) * dx2p(j  ) * dx3p(k+1))/12. &
                + (dx1p(i  ) * dx2p(j+1) * dx3p(k+1))/12. &
                + (dx1p(i+1) * dx2p(j+1) * dx3p(k+1))/4.
          ELSE IF (elem_type .EQ. 2) THEN ! HEX8
            vol = (dx1p(i  ) * dx2p(j  ) * dx3p(k  ))/8. &
                + (dx1p(i+1) * dx2p(j  ) * dx3p(k  ))/8. &
                + (dx1p(i  ) * dx2p(j+1) * dx3p(k  ))/8. &
                + (dx1p(i  ) * dx2p(j  ) * dx3p(k+1))/8. &
                + (dx1p(i+1) * dx2p(j+1) * dx3p(k  ))/8. &
                + (dx1p(i+1) * dx2p(j  ) * dx3p(k+1))/8. &
                + (dx1p(i  ) * dx2p(j+1) * dx3p(k+1))/8. &
                + (dx1p(i+1) * dx2p(j+1) * dx3p(k+1))/8.
          ELSE
            WRITE(*,*) 'WARNING: Error in computing volume force (elem_type)'
          END IF
          IF (dir .EQ. 1) fd(i,j,k,1) = fd(i,j,k,1) / (vol * L_ref**3.)
          IF (dir .EQ. 2) fd(i,j,k,2) = fd(i,j,k,2) / (vol * L_ref**3.)
          IF (dir .EQ. 3) fd(i,j,k,3) = fd(i,j,k,3) / (vol * L_ref**3.)
        END DO
      END DO
    END DO
  ELSE
    k = 1
      DO j = S2p, N2p
        DO i = S1p, N1p
          IF (elem_type .EQ. 1) THEN ! TRI3
            vol = (dx1p(i  ) * dx2p(j  ))/3. &
                + (dx1p(i+1) * dx2p(j  ))/6. &
                + (dx1p(i  ) * dx2p(j+1))/6. &
                + (dx1p(i+1) * dx2p(j+1))/3.
          ELSE IF (elem_type .EQ. 2) THEN ! QUAD4
            vol = (dx1p(i  ) * dx2p(j  ))/4. &
                + (dx1p(i+1) * dx2p(j  ))/4. &
                + (dx1p(i  ) * dx2p(j+1))/4. &
                + (dx1p(i+1) * dx2p(j+1))/4.
          ELSE
            WRITE(*,*) 'WARNING: Error in computing volume force (elem_type)'
          END IF
          IF (dir .EQ. 1) fd(i,j,k,1) = fd(i,j,k,1) / (vol * L_ref**2.)
          IF (dir .EQ. 2) fd(i,j,k,2) = fd(i,j,k,2) / (vol * L_ref**2.)
          IF (dir .EQ. 3) WRITE(*,*) 'WARNING: There is no 3rd direction in a 2D force.'
        END DO
      END DO
    ! end k
  END IF


  END SUBROUTINE residual2volume_force


  !SUBROUTINE set_pointers_to_non_allocatable_arrays()

  !IMPLICIT NONE

  !!===========================================================================================================
  !!=== mod_vars ==============================================================================================
  !!===========================================================================================================
  !!===========================================================================================================
  !!=== physikalische Parameter ===============================================================================
  !!===========================================================================================================
  !!--- RK coefficients ---
  !arkptr = C_LOC(ark_c(1))
  !brkptr = C_LOC(brk_c(1))
  !stride_largeptr = C_LOC(stride_large(1))
  !stride_medptr = C_LOC(stride_med(1))
  !stride_smallptr = C_LOC(stride_small(1))
  !!============================================================================================================
  !!=== iteration parameters ===================================================================================
  !!============================================================================================================
  !init_preptr = C_LOC(init_pre(1))
  !init_velptr = C_LOC(init_vel(1))
  !precoffset0ptr = C_LOC(precOffset0(1))
  !precratio0ptr = C_LOC(precRatio0(1))
  !impl_dirptr = C_LOC(impl_dir(1))

  !!--- Addt'l C pointers --------------------------------------------------------------------------------------
  !ibptr = C_LOC(iB(1,1))
  !outletptr = C_LOC(outlet(1,1,1))

  !!===========================================================================================================
  !!=== usr_vars ==============================================================================================
  !!===========================================================================================================
  !fringe_centerptr = C_LOC(fringe_center(1))
  !wk_flow_centerptr = C_LOC(WK_flow_center(1))
  !wk_preptr = C_LOC(WK_pre(1))

  !END SUBROUTINE set_pointers_to_non_allocatable_arrays




  FUNCTION rhs_L2_norm() RESULT(fn_val)

  IMPLICIT NONE

  REAL                   ::  fn_val
  INTEGER                ::  i,j,k,ii,jj,kk
  REAL                   ::  rhs1, rhs2, rhs3
  REAL                   ::  L2_norm_rhs, L2_norm_rhs_global

  L2_norm_rhs = 0.
  rhs1 = 0.; rhs2 = 0.; rhs3 = 0.

  DO k = S3p, N3p
    DO j = S2p, N2p
      DO i = S1p, N1p
        rhs1 = cIup(d1L,i)*rhs(i+d1L,j,k,1)
        DO ii = d1L+1,d1U
          rhs1 = rhs1 + cIup(ii,i)*rhs(i+ii,j,k,1)
        END DO
        rhs2 = cIvp(d2L,j)*rhs(i,j+d2L,k,2)
        DO jj = d2L+1,d2U
          rhs2 = rhs2 + cIvp(jj,j)*rhs(i,j+jj,k,2)
        END DO
        rhs3 = cIwp(d3L,k)*rhs(i,j,k+d3L,3)
        DO kk = d3L+1,d3U
          rhs3 = rhs3 + cIwp(kk,k)*rhs(i,j,k+kk,3)
        END DO
        L2_norm_rhs = L2_norm_rhs + res(i,j,k)*(rhs1**2. + rhs2**2. + rhs3**2.)
      END DO
    END DO
  END DO

  CALL MPI_ALLREDUCE(L2_norm_rhs,L2_norm_rhs_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  fn_val = SQRT(L2_norm_rhs_global)



  END FUNCTION rhs_L2_norm



  FUNCTION force_L2_norm() RESULT(fn_val)

  IMPLICIT NONE

  REAL                   ::  fn_val
  INTEGER                ::  i,j,k,ii,jj,kk
  REAL                   ::  fd1, fd2, fd3
  REAL                   ::  L2_norm_fd, L2_norm_fd_global

  L2_norm_fd = 0.
  fd1 = 0.; fd2 = 0.; fd3 = 0.

  ! we interpolate back to the pressure grid here
  ! since the force density is moved to the velocity
  ! grid right after the transfer
  DO k = S3p, N3p
    DO j = S2p, N2p
      DO i = S1p, N1p
        fd1 = cIup(d1L,i)*fd(i+d1L,j,k,1)
        DO ii = d1L+1,d1U
          fd1 = fd1 + cIup(ii,i)*fd(i+ii,j,k,1)
        END DO
        fd2 = cIvp(d2L,j)*fd(i,j+d2L,k,2)
        DO jj = d2L+1,d2U
          fd2 = fd2 + cIvp(jj,j)*fd(i,j+jj,k,2)
        END DO
        fd3 = cIwp(d3L,k)*fd(i,j,k+d3L,3)
        DO kk = d3L+1,d3U
          fd3 = fd3 + cIwp(kk,k)*fd(i,j,k+kk,3)
        END DO
        L2_norm_fd = L2_norm_fd + res(i,j,k)*(fd1**2. + fd2**2. + fd3**2.)
      END DO
    END DO
  END DO

  CALL MPI_ALLREDUCE(L2_norm_fd,L2_norm_fd_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  fn_val = SQRT(L2_norm_fd_global)*L_ref/(rho_fluid*U_ref**2.)


  END FUNCTION force_L2_norm


  FUNCTION vel_L2_norm() RESULT(fn_val)

  IMPLICIT NONE

  REAL                 ::  fn_val
  INTEGER              ::  i,j,k
  REAL                 ::  L2_norm_vel, L2_norm_vel_global

  L2_norm_vel = 0.

  DO k = S3p, N3p
    DO j = S2p, N2p
      DO i = S1p, N1p
        L2_norm_vel = L2_norm_vel + res(i,j,k)*(work1(i,j,k)**2. + work2(i,j,k)**2. + work3(i,j,k)**2.)
      END DO
    END DO
  END DO
    
  CALL MPI_ALLREDUCE(L2_norm_vel,L2_norm_vel_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  fn_val = SQRT(L2_norm_vel_global)

  END FUNCTION vel_L2_norm

  !> subroutine to reset the respective force component to 0
  !! will be called inside the transfer routine prior to 
  !! setting new data. Need to pass the component because 
  !! transferring is done on a per-component basis
  SUBROUTINE reset_force_density_component(dir)

  IMPLICIT NONE

  INTEGER, INTENT(IN)   ::  dir

  ! Do the reset on the fluid grid since these are the last 
  ! grids that the values were stored on
  SELECT CASE (dir)
    CASE (1)
      fd(S11B:N11B,S21B:N21B,S31B:N31B,1) = 0.
    CASE (2)
      fd(S12B:N12B,S22B:N22B,S32B:N32B,2) = 0.
    CASE (3)
      fd(S13B:N13B,S23B:N23B,S33B:N33B,3) = 0.
    CASE DEFAULT
      write(*,*) "WARNING: Something went wrong while resetting the force density component:",dir
  END SELECT 

  END SUBROUTINE reset_force_density_component


  SUBROUTINE open_log_iterations_fortran()

  IMPLICIT NONE

  OPEN(10,FILE='log_iterations.txt',STATUS='OLD',ACCESS='APPEND')

  END SUBROUTINE open_log_iterations_fortran

  SUBROUTINE close_log_iterations_fortran()

  IMPLICIT NONE

  CLOSE(10)

  END SUBROUTINE close_log_iterations_fortran


  SUBROUTINE save_old_velocity()

  IMPLICIT NONE

  INTEGER          ::  m
  
  DO m = 1,dimens
    IF(m .EQ. 1) vel_old(S11B:N11B,S21B:N21B,S31B:N31B,1) = vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
    IF(m .EQ. 2) vel_old(S12B:N12B,S22B:N22B,S32B:N32B,2) = vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
    IF(m .EQ. 3) vel_old(S13B:N13B,S23B:N23B,S33B:N33B,3) = vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
  END DO

  END SUBROUTINE save_old_velocity

  SUBROUTINE restore_old_velocity()

  IMPLICIT NONE

  INTEGER          ::  m

  DO m = 1,dimens
    IF(m .EQ. 1) vel(S11B:N11B,S21B:N21B,S31B:N31B,1) = vel_old(S11B:N11B,S21B:N21B,S31B:N31B,1)
    IF(m .EQ. 2) vel(S12B:N12B,S22B:N22B,S32B:N32B,2) = vel_old(S12B:N12B,S22B:N22B,S32B:N32B,2)
    IF(m .EQ. 3) vel(S13B:N13B,S23B:N23B,S33B:N33B,3) = vel_old(S13B:N13B,S23B:N23B,S33B:N33B,3)
  END DO

  END SUBROUTINE restore_old_velocity

END MODULE usr_func
