!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                              							     *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!*************************************************************************************************************

!> @file usr_config.f90
!! File holding subroutine configuration. User-defined parameter values can be set here.

!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> This subroutine is called after MPI and HDF5 have been initialized and IMPACT checked for 'queue.txt'.
  !! It uses mod_dims, mod_vars, usr_vars, usr_func.
  SUBROUTINE configuration
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE usr_vars
  USE usr_func
  
  IMPLICIT NONE
  
  !bbecsek
  CHARACTER(len=80)        ::  text
  CHARACTER(len=80)        ::  dummy
  INTEGER                  ::  ios

  !************************** READ config.txt FILE **********************************************************!
  
  OPEN(10,FILE=TRIM('config.txt'),FORM='formatted',ACTION='read',STATUS='old',IOSTAT=ios)

  IF (ios /= 0) THEN
    IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open config.txt file!'
    CALL MPI_FINALIZE(merror)
    STOP
  END IF

  DO
    READ(10,FMT=*,IOSTAT=ios) text
    IF (text == 'BEGIN') EXIT
  END DO

  DO
    READ(10,FMT='(a)',IOSTAT=ios) text
    
    IF (ios /= 0) THEN
      IF (rank == 0) WRITE(*,*) 'ERROR! Cannot read the config.txt file!'
      CALL MPI_FINALIZE(merror)
      STOP
    END IF

    IF (INDEX(text,'#') == 1) CYCLE

  !***********************************************************************************************************
  !--- all specifications of the configuration ---
  ! note: - mandatory specifications
  !       - additional, user-defined specifications (cf. end of this subroutine)
  !
  
  
  !===========================================================================================================
  !=== job control ===========================================================================================
  !===========================================================================================================
  !--- tasks ---
  ! task = 1 (time integration)
  ! task = 2 (post-processing)
  ! task = 3 (Eigenvalue problems, under construction)
  ! task = 4 (other matrix analysis, under construction)
  !
  !task = 1
    IF (INDEX(text,'task'                        ) == 1) READ(UNIT=text,FMT=*) dummy, task
  
  !--- restart number ---
  !restart = 0
    IF (INDEX(text,'restart_from'                ) == 1) READ(UNIT=text,FMT=*) dummy, restart
  
  
  
  !===========================================================================================================
  !=== resolution ============================================================================================
  !===========================================================================================================
  !--- max. number of time steps (independend of restart count) ---
  ! note: - time integration is terminated when "n_timesteps" or "time_end" is reached
  !
  !n_timesteps = 200000000
    IF (INDEX(text,'n_timesteps'                 ) == 1) READ(UNIT=text,FMT=*) dummy, n_timesteps
  
#ifdef ALLOC
  !--- number of grid points (global) ---
  ! note: - multigrid works most efficiently if Mi = a*2**q+1 and a is a "small" integer
  !       - Mi must be sufficiently large such that at least one central finite difference stencil fits in
  !         (will be improved in future releases)
  !       - for 2D simulations set M3 = 2 (all terms in the third direction are switched off)
  !
  !M1 = 8*2**6+1
    IF (INDEX(text,'M1'                          ) == 1) READ(UNIT=text,FMT=*) dummy, M1
  !M2 = 1*2**6+1
    IF (INDEX(text,'M2'                          ) == 1) READ(UNIT=text,FMT=*) dummy, M2
  !M3 = 2 
    IF (INDEX(text,'M3'                          ) == 1) READ(UNIT=text,FMT=*) dummy, M3
  
  !--- numbers of processor blocks ---
  ! note: - multigrid requires that MOD(Mi-1,NBi) = 0
  !       - 2D simulations (M3 = 2) require NB3 = 1
  !
  !NB1 = 32
    IF (INDEX(text,'NB1'                         ) == 1) READ(UNIT=text,FMT=*) dummy, NB1
  !NB2 = 1
    IF (INDEX(text,'NB2'                         ) == 1) READ(UNIT=text,FMT=*) dummy, NB2
  !NB3 = 1
    IF (INDEX(text,'NB3'                         ) == 1) READ(UNIT=text,FMT=*) dummy, NB3
#endif
  
  
  !===========================================================================================================
  !=== non-dimensional numbers ===============================================================================
  !===========================================================================================================

  !--- Reynolds number ---
  !Re = 1500. 
    IF (INDEX(text,'Re'                          ) == 1) READ(UNIT=text,FMT=*) dummy, Re
  
  
  !===========================================================================================================
  !=== extents ===============================================================================================
  !===========================================================================================================
  !--- start time ---
  !time_start = 0.
    IF (INDEX(text,'time_start'                  ) == 1) READ(UNIT=text,FMT=*) dummy, time_start
  
  !--- end time ---
  ! note: - time integration is terminated when "n_timesteps" or "time_end" is reached
  !
  !time_end = 20.
    IF (INDEX(text,'time_end'                    ) == 1) READ(UNIT=text,FMT=*) dummy, time_end
  
  !--- extents of the physical domain ---
  !L1 = 8.
    IF (INDEX(text,'L1'                          ) == 1) READ(UNIT=text,FMT=*) dummy, L1
  !L2 = 1.
    IF (INDEX(text,'L2'                          ) == 1) READ(UNIT=text,FMT=*) dummy, L2
  !L3 = 1.
    IF (INDEX(text,'L3'                          ) == 1) READ(UNIT=text,FMT=*) dummy, L3
  
    IF (INDEX(text,'refine_2dir_yes'             ) == 1) READ(UNIT=text,FMT=*) dummy, refine_2dir_yes

    IF (INDEX(text,'refine_2dir_xi'             ) == 1) READ(UNIT=text,FMT=*) dummy, refine_2dir_xi

  
  
  !===========================================================================================================
  !=== configuration =========================================================================================
  !===========================================================================================================
  !--- velocity boundary conditions (BC) ---
  ! note: - concentration boundary conditions are controlled via "isopycnal" and may depend on the velocity
  !       - for advective boundary conditions you have to specify Dirichlet boundary conditions
  !                           _
  ! symmetry BC:    BC = -2    |
  ! periodicity BC: BC = -1   _|- symmetric, central FD stencils
  ! Dirichlet BC:   BC =  1    |
  ! Neumann BC*:    BC =  2    |- non-symmetric, skew FD stencils
  ! Robin BC*:      BC =  3   _|
  ! *not yet implemented for velocity
  !
  !    orientation of boundary normal
  !    |lower/upper boundary
  !    ||
  ! BC_1L_global
  !
  !BC_1L_global = -1
    IF (INDEX(text,'BC_1L_global'                ) == 1) READ(UNIT=text,FMT=*) dummy, BC_1L_global
  !BC_1U_global = -1
    IF (INDEX(text,'BC_1U_global'                ) == 1) READ(UNIT=text,FMT=*) dummy, BC_1U_global

  !IF (1==2) THEN
  !  BC_1L_global = 1
  !  BC_1U_global = 1
  !END IF
  
  !BC_2L_global = 1
    IF (INDEX(text,'BC_2L_global'                ) == 1) READ(UNIT=text,FMT=*) dummy, BC_2L_global
  !BC_2U_global = 1
    IF (INDEX(text,'BC_2U_global'                ) == 1) READ(UNIT=text,FMT=*) dummy, BC_2U_global
  
  !BC_3L_global = -1
    IF (INDEX(text,'BC_3L_global'                ) == 1) READ(UNIT=text,FMT=*) dummy, BC_3L_global
  !BC_3U_global = -1
    IF (INDEX(text,'BC_3U_global'                ) == 1) READ(UNIT=text,FMT=*) dummy, BC_3U_global
  
  !--- advective velocity boundaries ---
  ! note: - requires BC_global = 1
  !
  !         orientation of the boundary normal
  !         |   lower/upper boundary
  !         |   |   velocity component
  !         |   |   |
  ! outlet(1:3,1:2,1:3)
  !
  !outlet          = .FALSE.
    IF (INDEX(text,'outlet11'                    ) == 1) READ(UNIT=text,FMT=*) dummy, outlet(1,1,1:3)
    IF (INDEX(text,'outlet12'                    ) == 1) READ(UNIT=text,FMT=*) dummy, outlet(1,2,1:3)
    IF (INDEX(text,'outlet21'                    ) == 1) READ(UNIT=text,FMT=*) dummy, outlet(2,1,1:3)
    IF (INDEX(text,'outlet22'                    ) == 1) READ(UNIT=text,FMT=*) dummy, outlet(2,2,1:3)
    IF (INDEX(text,'outlet31'                    ) == 1) READ(UNIT=text,FMT=*) dummy, outlet(3,1,1:3)
    IF (INDEX(text,'outlet32'                    ) == 1) READ(UNIT=text,FMT=*) dummy, outlet(3,2,1:3)

  !outlet(1,1,1:3) = .TRUE.
  !outlet(2,2,1:3) = .TRUE.
  
  !IF (1==2) THEN
  !  outlet(1,2,1:3) = .TRUE.
  !END IF
  
  
  !--- bulk flow forcing mode ---
  ! (under construction, currently limited to channel flows)
  ! note: - fixes bulk flow velocity to 2./3. in periodic channels
  !
  ! forcing_mode = 0 (no forcing)
  ! forcing_mode = 1 (explicit forcing in time)
  ! forcing_mode = 2 (implicit forcing in time)
  !
  !forcing_mode = 0
    IF (INDEX(text,'forcing_mode'                ) == 1) READ(UNIT=text,FMT=*) dummy, forcing_mode
  
  !--- direction of bulk flow (relevant only for bulk flow forcing) ---
  ! (under construction, currently limited to channel flows)
  !
  !bulkflow_dir = 1
    IF (INDEX(text,'bulkflow_dir'                ) == 1) READ(UNIT=text,FMT=*) dummy, bulkflow_dir
  
  !--- bulk flow velocity ---
  !vel_bulk = 2./3.
    IF (INDEX(text,'vel_bulk'                    ) == 1) READ(UNIT=text,FMT=*) dummy, vel_bulk
  
  
  !===========================================================================================================
  !=== time discretization ===================================================================================
  !===========================================================================================================
  !--- time integration scheme ---
  ! timeint_mode = 0 (CN-RK3)
  ! timeint_mode = 1 (RK3-O3)
  !
  !timeint_mode = 1
    IF (INDEX(text,'timeint_mode'                ) == 1) READ(UNIT=text,FMT=*) dummy, timeint_mode
  
  !--- Crank-Nicolson factor ---
  ! note: - must be 0. < thetaL < 1.
  !
  ! theta = 0. (fully explicit time integration of viscous terms)
  ! theta = 1. (fully implicit time integration of viscous terms)
  !
  !thetaL = 0.5
    IF (INDEX(text,'thetaL'                      ) == 1) READ(UNIT=text,FMT=*) dummy, thetaL
  
  !--- Euler flow (no physical viscosity) ---
  !Euler_yes = .FALSE.
    IF (INDEX(text,'Euler_yes'                   ) == 1) READ(UNIT=text,FMT=*) dummy, Euler_yes
  
  !--- Stokes flow (no fluid advection in the domain) ---
  !Stokes_yes = .FALSE.
    IF (INDEX(text,'Stokes_yes'                  ) == 1) READ(UNIT=text,FMT=*) dummy, Stokes_yes
  
  !--- two-step time integration (relevant only for CN-RK3) ---
  !twostep_yes = .FALSE.
    IF (INDEX(text,'twostep_yes'                 ) == 1) READ(UNIT=text,FMT=*) dummy, twostep_yes
  
  !--- CFL number ---
  ! note: - normalized with maximum stable number, i.e. should be 0 < CFL < 1
  !       - applies only to advective and viscous terms of velocity and concentrations but not on other user-
  !         defined volume forces
  !       - full RK3 and CN-RK3 stability domains are taken into account
  !       - symbols of the spatial discretization are taken into account 
  !
  !CFL = 0.75
    IF (INDEX(text,'CFL'                         ) == 1) READ(UNIT=text,FMT=*) dummy, CFL

  !--- max. time step size ---
  !dtime_max = 5.*10.**(-1)
    IF (INDEX(text,'dtime_max'                   ) == 1) READ(UNIT=text,FMT=*) dummy, dtime_max
  
  !--- max. time step size at beginning of time integration ---
  ! note: - "dtime" is limited by "dtime0" for the first "Int_dtime" time steps
  !
  !dtime0 = 10.**(-8)
    IF (INDEX(text,'dtime0'                      ) == 1) READ(UNIT=text,FMT=*) dummy, dtime0
  
  !--- number of time steps after which time step size is recomputed ---
  !Int_dtime = 500
    IF (INDEX(text,'Int_dtime'                   ) == 1) READ(UNIT=text,FMT=*) dummy, Int_dtime
  
  
  
  !===========================================================================================================
  !=== spatial discretization ================================================================================
  !===========================================================================================================
  !--- upwinding for advective terms of velocity momentum equation ---
  ! note: - option not available for compact finite differences (cf. comp_conv_yes)
  !
  !upwind_yes = .TRUE.
    IF (INDEX(text,'upwind_yes'                  ) == 1) READ(UNIT=text,FMT=*) dummy, upwind_yes
  
  
  !--- computation of finite difference coefficients (relevant only for explicit finite differences) ---
  ! note: - affects the accuracy and temporal stability of the discretization
  !
  ! mapping_yes = .TRUE.  (computation on computational, equidistant grid and subsequent mapping on
  !                        physical grid)
  ! mapping_yes = .FALSE. (direct computation on physical grid)
  !
  !mapping_yes = .TRUE.  
    IF (INDEX(text,'mapping_yes'                 ) == 1) READ(UNIT=text,FMT=*) dummy, mapping_yes
  
  
  !===========================================================================================================
  !=== iterative solver ======================================================================================
  !===========================================================================================================
  !--- infinity norm of velocity divergence ---
  ! note: - cf. comments on "weighting_yes"
  !
  !epsU = 10.**(-6)
    IF (INDEX(text,'epsU'                        ) == 1) READ(UNIT=text,FMT=*) dummy, epsU
  
  ! weighting_yes = .TRUE.  (div(vel) is normalized with the local grid spacing)
  ! weighting_yes = .FALSE. (div(vel) is not normalized)
  !
  !weighting_yes = .FALSE.
    IF (INDEX(text,'weighting_yes'               ) == 1) READ(UNIT=text,FMT=*) dummy, weighting_yes
  
  !--- number of time steps after which the mean pressure is set to zero ---
  ! note: - the mean pressure level is not implicitly fixed in the pressure solver
  !
  !Int_lev_pre = 1
    IF (INDEX(text,'Int_lev_pre'                 ) == 1) READ(UNIT=text,FMT=*) dummy, Int_lev_pre
  
  !--- preconditioners ---
  ! precond            = 0: no preconditioner, i.e. only primary solvers (Richardson iteration or BiCGstab)
  ! precond_outer      = 1: Richardson iteration with Laplace preconditioner
  ! precond_outer      = 2: Richardson iteration with commuation-based preconditioner
  ! precond_Poisson    = 1: BiCGstab with V-cycle multigrid
  ! precond_Poisson    = 2: BiCGstab with F-cycle multigrid
  ! precond_Helmh_vel  = 1: BiCGstab with V-cycle multigrid
  ! precond_Helmh_conc = 1: BiCGstab with V-cycle multigrid
  !
  ! from top: outer pressure iteration, Helmholtz equation (velocity), Poisson equation,
  !           Helmholtz equations (concentrations)
  !
  !precond_outer      = 2
    IF (INDEX(text,'precond_outer'               ) == 1) READ(UNIT=text,FMT=*) dummy, precond_outer
  !precond_Poisson    = 1
    IF (INDEX(text,'precond_Poisson'             ) == 1) READ(UNIT=text,FMT=*) dummy, precond_Poisson
  !precond_Helmh_vel  = 0
    IF (INDEX(text,'precond_Helmh_vel'           ) == 1) READ(UNIT=text,FMT=*) dummy, precond_Helmh_vel
  
  !--- max. numbers of iterations ---
  ! from top: outer pressure iteration, Helmholtz equation (velocity), Poisson equation,
  !           Helmholtz equations (concentrations)
  !
  !n_it_outer      = 10
    IF (INDEX(text,'n_it_outer'                  ) == 1) READ(UNIT=text,FMT=*) dummy, n_it_outer
  !n_it_Poisson    = 10
    IF (INDEX(text,'n_it_Poisson'                ) == 1) READ(UNIT=text,FMT=*) dummy, n_it_Poisson
  !n_it_Helmh_vel  = 10
    IF (INDEX(text,'n_it_Helmh_vel'              ) == 1) READ(UNIT=text,FMT=*) dummy, n_it_Helmh_vel
  
  !--- initialization before implicit problems are solved in a Runge-Kutta sub-step ---
  ! note: - initialization of velocity and concentration applies only to CN-RK3 time integration
  !
  ! from top: pressure, velocity, concentrations
  !
  !init_pre (1:RK_steps) = .FALSE.
    IF (INDEX(text,'init_pre'                    ) == 1) READ(UNIT=text,FMT=*) dummy, init_pre(1:RK_steps)
  !init_vel (1:RK_steps) = .FALSE.
    IF (INDEX(text,'init_vel'                    ) == 1) READ(UNIT=text,FMT=*) dummy, init_vel(1:RK_steps)
  
  !--- accuracy of preconditioner problems in outer pressure iteration (relevant only for CN-RK3) ---
  ! note: - sets the solution accuracy of the inner Poisson problems relative to the expected accuracy of the
  !         outer pressure iteration cycle
  !
  !precOffset0(1:RK_steps) = 0.5
    IF (INDEX(text,'precOffset0'                 ) == 1) READ(UNIT=text,FMT=*) dummy, precOffset0(1:RK_steps)
  
  !--- expected convervengence ratios of the first time step of a restart (relevant only for CN-RK3) ---
  !precRatio0 (1:RK_steps) = 10.**(-4)
    IF (INDEX(text,'precRatio0'                  ) == 1) READ(UNIT=text,FMT=*) dummy, precRatio0(1:RK_steps)
  
  !--- use flux corrections on fine/coarse grids ---
  ! note: - all pressure-Poisson problems are singular (including all multigrid levels)
  !       - to guarantee at least one solution, each right-hand side of a pressure-Poisson problem must be
  !         fully covered by the column space of the corresponding pressure matrix, i.e. each right-hand side
  !         must be orthogonal to the left null space of the corresponding pressure matrix
  !       - setting these variables enforces corrections of the fluxes over the boundaries in each sub-time
  !         step such that the corresponding pressure-Poisson problems have at least one solution
  !       - "nullspace_coarse_yes" (for the multigrid problems) is normally not necessary
  !
  !nullspace_yes        = .TRUE.
    IF (INDEX(text,'nullspace_yes'               ) == 1) READ(UNIT=text,FMT=*) dummy, nullspace_yes
  !nullspace_coarse_yes = .FALSE.
    IF (INDEX(text,'nullspace_coars_yes'         ) == 1) READ(UNIT=text,FMT=*) dummy, nullspace_coarse_yes
  
  !--- read null space base vectors instead of recomputing them ---
  !read_nullspace_yes = .FALSE.
    IF (INDEX(text,'read_nullspace_yes'          ) == 1) READ(UNIT=text,FMT=*) dummy, read_nullspace_yes
  
  !--- specify if the fluxes are corrected by means of orthogonal projection ---
  ! note: - a vector of the left null space is used as correction vector
  !       - leads to minimal L2-norms of the flux corrections
  !       - fluxes are corrected on the entire boundary (excluding symmetry boundaries)
  !       - if orthogonal projection is not used you need to specify the correction vector "th" in
  !         "usr_initcond.f90"
  !
  !nullspace_ortho_yes = .TRUE.
    IF (INDEX(text,'nullspace_ortho_yes'         ) == 1) READ(UNIT=text,FMT=*) dummy, nullspace_ortho_yes
  
  
  
  !===========================================================================================================
  !=== multigrid =============================================================================================
  !===========================================================================================================
  !--- numbers of relaxation sweeps ---
  ! from top: restriction, prolongation, coarsest grid
  !
  !n_relax_down   = 4
    IF (INDEX(text,'n_relax_down'                ) == 1) READ(UNIT=text,FMT=*) dummy, n_relax_down
  !n_relax_up     = 4
    IF (INDEX(text,'n_relax_up'                  ) == 1) READ(UNIT=text,FMT=*) dummy, n_relax_up
  !n_relax_bottom = 4
    IF (INDEX(text,'n_relax_bottom'              ) == 1) READ(UNIT=text,FMT=*) dummy, n_relax_bottom
  
  !--- directions of (alternating) line-relaxations ---
  ! note: - line-relaxation is necessary if the grid is strongly unisotropic, i.e. when the grid spacings
  !         differ strongly in the different spatial directions
  !
  ! impl_dir = 0 ((partially) red-black-ordered point-block relaxation)
  ! impl_dir = 1 (lexicographically-ordered line-relaxation)
  !
  !impl_dir(1:3) = 0
    IF (INDEX(text,'impl_dir'                    ) == 1) READ(UNIT=text,FMT=*) dummy, impl_dir
  
  !--- Jacobi smoothing instead of Gauss-Seidel ---
  !Jacobi_yes = .FALSE.
    IF (INDEX(text,'Jacobi_yes'                  ) == 1) READ(UNIT=text,FMT=*) dummy, Jacobi_yes
  
  !--- max. number of coarse grid levels ---
  ! note: - must be 1 <= n_grids_limit <= n_grids_max
  !
  !n_grids_limit = n_grids_max
    IF (INDEX(text,'n_grids_limit'               ) == 1) READ(UNIT=text,FMT=*) dummy, n_grids_limit 
  
  
  !===========================================================================================================
  !=== output ================================================================================================
  !===========================================================================================================
  !--- write standard output ---
  !write_stout_yes = .TRUE.
    IF (INDEX(text,'write_stout_yes'             ) == 1) READ(UNIT=text,FMT=*) dummy, write_stout_yes
  
  !--- write log file for iterations ---
  !log_iteration_yes = .TRUE.
    IF (INDEX(text,'log_iteration_yes'           ) == 1) READ(UNIT=text,FMT=*) dummy, log_iteration_yes
  
  !--- write restart files ---
  !write_restart_yes = .FALSE.
    IF (INDEX(text,'write_restart_yes'           ) == 1) READ(UNIT=text,FMT=*) dummy, write_restart_yes
  
  !--- time interval of field output ---
  !dtime_out_vect = 0.02
    IF (INDEX(text,'dtime_out_vect'              ) == 1) READ(UNIT=text,FMT=*) dummy, dtime_out_vect
  
  !--- time interval of other ouput (statistics) ---
  !dtime_out_scal = 1.
    IF (INDEX(text,'dtime_out_scal'              ) == 1) READ(UNIT=text,FMT=*) dummy, dtime_out_scal
  
  !--- time interval of other ouput (kalman filter) ---
  !dtime_out_kalm = 1.
    IF (INDEX(text,'dtime_out_kalm'              ) == 1) READ(UNIT=text,FMT=*) dummy, dtime_out_kalm
  
  !--- compute and write lambda2 fields from velocity ---
  !write_lambda2_yes = .TRUE.
    IF (INDEX(text,'write_lambda2_yes'           ) == 1) READ(UNIT=text,FMT=*) dummy, write_lambda2_yes
 
  !--- bbecsek: ---
  !--- output force density received from FSI coupling (Moonolith) ---
    IF (INDEX(text,'write_force_yes'             ) == 1) READ(UNIT=text,FMT=*) dummy, write_force_yes
  ! --- writing covariance into xdmf file  
    IF (INDEX(text,'write_covariance_yes'        ) == 1) READ(UNIT=text,FMT=*) dummy, write_covariance_yes
  !--- output .xmf files for reading in Paraview and Visit ---
    IF (INDEX(text,'write_xdmf_yes'             ) == 1) READ(UNIT=text,FMT=*) dummy, write_xdmf_yes

  !--- scale output with U_ref and L_ref to make them dimensional ---
    IF (INDEX(text,'scale_output_yes'           ) == 1) READ(UNIT=text,FMT=*) dummy, scale_output_yes

  !--- coarsening factors for 3D fields ---
  ! note: - helpful to save disc space
  !       - 2D data are written without strides
  !       - stride = 0 indicates no output
  !
  !stride_large(1:3) = 1
    IF (INDEX(text,'stride_large'                ) == 1) READ(UNIT=text,FMT=*) dummy, stride_large(1:3)
  !stride_med  (1:3) = 2
    IF (INDEX(text,'stride_med'                  ) == 1) READ(UNIT=text,FMT=*) dummy, stride_med(1:3)
  !stride_small(1:3) = 0
    IF (INDEX(text,'stride_small'                ) == 1) READ(UNIT=text,FMT=*) dummy, stride_small(1:3)
  
  !--- write debugging files ---
  ! (deprecated, will vanish in future releases)
  !
  !write_test_yes = .FALSE.
    IF (INDEX(text,'write_test_yes'              ) == 1) READ(UNIT=text,FMT=*) dummy, write_test_yes
  
  
  
  !===========================================================================================================
  !=== additional, user-defined input ========================================================================
  !===========================================================================================================
  !--- pi ---
  !pi = 2.*ABS(ACOS(0.)) !--- hard-coded
  
  !--- direction of wall normal ---
  !wallnorm_dir = 1
    IF (INDEX(text,'wallnorm_dir'                ) == 1) READ(UNIT=text,FMT=*) dummy, wallnorm_dir
  
  !--- max. number of wavenumbers of discrete Fourier transform ---
  !amax = 4
    IF (INDEX(text,'amax'                        ) == 1) READ(UNIT=text,FMT=*) dummy, amax
  !bmax = 4
    IF (INDEX(text,'bmax'                        ) == 1) READ(UNIT=text,FMT=*) dummy, bmax
  
  !--- origin of coordinate system ---
  !y1_origin = 0.
    IF (INDEX(text,'y1_origin'                   ) == 1) READ(UNIT=text,FMT=*) dummy, y1_origin
  !y2_origin = 0.
    IF (INDEX(text,'y2_origin'                   ) == 1) READ(UNIT=text,FMT=*) dummy, y2_origin
  !y3_origin = 0.
    IF (INDEX(text,'y3_origin'                   ) == 1) READ(UNIT=text,FMT=*) dummy, y3_origin

  !===========================================================================================================
  !=== bbecsek: Fringe Region Forcing ========================================================================
  !===========================================================================================================
  !--- Switch on/off fringe forcing ---
    IF (INDEX(text,'fringe_yes'                  ) == 1) READ(UNIT=text,FMT=*) dummy, fringe_yes

  !--- Starting location of fringe forcing ---
    IF (INDEX(text,'fringe_start'                ) == 1) READ(UNIT=text,FMT=*) dummy, fringe_start

  !--- Ending location of fringe forcing ---
    IF (INDEX(text,'fringe_end'                  ) == 1) READ(UNIT=text,FMT=*) dummy, fringe_end

  !--- Fringe rise interval ---
    IF (INDEX(text,'fringe_rise'                 ) == 1) READ(UNIT=text,FMT=*) dummy, fringe_rise

  !--- Fringe fall interval ---
    IF (INDEX(text,'fringe_fall'                 ) == 1) READ(UNIT=text,FMT=*) dummy, fringe_fall

  !--- Effective direction of fringe region ---
    IF (INDEX(text,'fringe_dir'                  ) == 1) READ(UNIT=text,FMT=*) dummy, fringe_dir
  
  !--- Fringe region forcing amplitude ---
    IF (INDEX(text,'fringe_amp'                  ) == 1) READ(UNIT=text,FMT=*) dummy, fringe_amp
 
  !--- Fringe region perpendicular center ---
    IF (INDEX(text,'fringe_center'               ) == 1) READ(UNIT=text,FMT=*) dummy, fringe_center(1:3)

  !--- Fringe region perpendicular radius ---
    IF (INDEX(text,'fringe_radius'               ) == 1) READ(UNIT=text,FMT=*) dummy, fringe_radius

  !--- Whether to ramp up forcing gradually ---
    IF (INDEX(text,'step_yes'                    ) == 1) READ(UNIT=text,FMT=*) dummy, step_yes

  !--- Interval of ramp-up ---
    IF (INDEX(text,'t_rampup'                    ) == 1) READ(UNIT=text,FMT=*) dummy, t_rampup

  !--- Parabolic shape of forcing ---
    IF (INDEX(text,'parabola_yes'                ) == 1) READ(UNIT=text,FMT=*) dummy, parabola_yes

  !===========================================================================================================
  !=== bbecsek: Immersed Boundary Method =====================================================================
  !===========================================================================================================
  !--- Quadi FEM representation of IMPACT grid ---
  !elem_type = 1
  ! 2D: 1 -> TRI, 2 -> QUAD
  ! 3D: 1 -> TET, 2 -> HEX
    IF (INDEX(text,'elem_type'                   ) == 1) READ(UNIT=text,FMT=*) dummy, elem_type

 !--- dynamic viscosity of blood (or other fluid) [Pa*s] Literature Value: ~0.01
  !mu_fluid = 0.01
    IF (INDEX(text,'mu_fluid'                    ) == 1) READ(UNIT=text,FMT=*) dummy, mu_fluid

  !--- bulk density of blood (or other fluid) [kg/m3]
  !rho_fluid = 1.
    IF (INDEX(text,'rho_fluid'                   ) == 1) READ(UNIT=text,FMT=*) dummy, rho_fluid

  !--- reference length scale --- [m]
  !L_ref = 1.
    IF (INDEX(text,'L_ref'                       ) == 1) READ(UNIT=text,FMT=*) dummy, L_ref

  !--- reference velocity --- [m/s]
  !U_ref = 1.
    IF (INDEX(text,'U_ref'                       ) == 1) READ(UNIT=text,FMT=*) dummy, U_ref

  !==========================================================================================================
  !=== bbecsek: Lumped Parameter Windkessel Loading Model ===================================================
  !==========================================================================================================
  !--- Switch on/off WK ---
    IF (INDEX(text,'WK_yes'                      ) == 1) READ(UNIT=text,FMT=*) dummy, WK_yes
  
  !--- Peripheral resistance --- [Pa*s/m3]
    IF (INDEX(text,'R_p'                         ) == 1) READ(UNIT=text,FMT=*) dummy, R_p

  !--- Characteristic resistance --- [Pa*s/m3]
    IF (INDEX(text,'R_c'                         ) == 1) READ(UNIT=text,FMT=*) dummy, R_c

  !--- Arterial compliance --- [m3/Pa]
    IF (INDEX(text,'C_art'                       ) == 1) READ(UNIT=text,FMT=*) dummy, C_art

  !--- Inertance --- [Pa * s2/m3]
    IF (INDEX(text,'L_art'                     ) == 1) READ(UNIT=text,FMT=*) dummy, L_art

  !--- Windkessel model type 3 == three-element, 4 == four-element ---
    IF (INDEX(text,'WK_type'                     ) == 1) READ(UNIT=text,FMT=*) dummy, WK_type

  !--- Direction in which the flow for Windkessel model is computed ---
    IF (INDEX(text,'WK_flow_dir'                 ) == 1) READ(UNIT=text,FMT=*) dummy, WK_flow_dir

  !--- Position in flow direction at which the flow is computed ---
    IF (INDEX(text,'WK_flow_pos'                 ) == 1) READ(UNIT=text,FMT=*) dummy, WK_flow_pos

  !--- Center of the circle (3D) or line (2D) through which the flow is computed ---
    IF (INDEX(text,'WK_flow_center'              ) == 1) READ(UNIT=text,FMT=*) dummy, WK_flow_center(1:3)

  !--- Radius of the circle (3D) or half-lenght of line (2D) ---
    IF (INDEX(text,'WK_flow_radius'              ) == 1) READ(UNIT=text,FMT=*) dummy, WK_flow_radius

  !--- Initial conditions for the outlet pressure ---
    IF (INDEX(text,'WK_pre'                      ) == 1) READ(UNIT=text,FMT=*) dummy, WK_pre(1:3)

  !--- Winkessel fringe region start location ---
    IF (INDEX(text,'WK_frge_start'             ) == 1) READ(UNIT=text,FMT=*) dummy, WK_frge_start

  !--- Windkessel fringe region end location ---
    IF (INDEX(text,'WK_frge_end'               ) == 1) READ(UNIT=text,FMT=*) dummy, WK_frge_end

  !--- Windkessel fringe region rise interval ---
    IF (INDEX(text,'WK_frge_rise'              ) == 1) READ(UNIT=text,FMT=*) dummy, WK_frge_rise

  !--- Windkessel fringe region fall interval ---
    IF (INDEX(text,'WK_frge_fall'              ) == 1) READ(UNIT=text,FMT=*) dummy, WK_frge_fall

  !--- Windkessel fringe region amplitude ---
    IF (INDEX(text,'WK_frge_amp'               ) == 1) READ(UNIT=text,FMT=*) dummy, WK_frge_amp

    IF (INDEX(text,'vel_initcond_file_yes'     ) == 1) READ(UNIT=text,FMT=*) dummy, vel_initcond_file_yes


    IF (INDEX(text,'END'                         ) == 1) EXIT
    
  END DO

  CLOSE(10)

  
END SUBROUTINE configuration
  
