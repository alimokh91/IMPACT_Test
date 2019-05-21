!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

MODULE usr_vars
  
  USE mod_dims
  USE mod_vars
  USE ISO_C_BINDING !bbecsek
  
  IMPLICIT NONE
  
  !--- specification of additional, user-defined variables ---
  !
  REAL                   ::  pi
  INTEGER                ::  wallnorm_dir
  INTEGER                ::  amax, bmax
  
  REAL                   ::  energy_visc
  REAL                   ::  diss_viscInt_old
  
  !--- origin of coordinate system ---------------------------------------------------------------------------
  REAL                   ::  y1_origin
  REAL                   ::  y2_origin
  REAL                   ::  y3_origin

  !--- quasi-FEM options ---
  INTEGER                ::  elem_type !> 2D: 1 -> TRI, 2 -> QUAD
                                       !! 3D: 1 -> TET, 2 -> HEX

  !--- Fringe forcing region specs ---------------------------------------------------------------------------
  LOGICAL                ::  fringe_yes
  REAL                   ::  fringe_start
  REAL                   ::  fringe_end
  REAL                   ::  fringe_rise
  REAL                   ::  fringe_fall
  INTEGER                ::  fringe_dir
  REAL, TARGET           ::  fringe_center(1:3)
  REAL                   ::  fringe_amp

  REAL                   ::  fringe_radius
  LOGICAL                ::  step_yes
  REAL                   ::  t_rampup

  LOGICAL                ::  parabola_yes
  LOGICAL                ::  D3_channel_parab_yes
  LOGICAL                ::  pressure_yes
  REAL                   ::  p_freq !Hz
  REAL                   ::  p_amp !pascal
  !--- Windkessel specifications ------------------------------------------------------------------------------
  REAL                   ::  Q_AV
  REAL                   ::  dQ_AV_dt
  REAL                   ::  Q_AV_0, Q_AV_1, Q_AV_2, Q_AV_3 ! helper variables
  LOGICAL                ::  WK_yes
  INTEGER                ::  WK_type
  INTEGER                ::  WK_flow_dir
  REAL                   ::  WK_flow_pos
  REAL, TARGET           ::  WK_flow_center(1:3)
  REAL                   ::  WK_flow_radius
  REAL                   ::  WK_predot(1:3)
  REAL                   ::  R_c   ! characteristic resistance 
  REAL                   ::  R_p   ! peripheral resistance
  REAL                   ::  C_art ! arterial compliance
  REAL                   ::  L_art ! inertance

  REAL                   ::  WK_frge_start
  REAL                   ::  WK_frge_end
  REAL                   ::  WK_frge_rise
  REAL                   ::  WK_frge_fall
  REAL                   ::  WK_frge_amp

  REAL, TARGET           ::  WK_pre(1:3)
#ifndef FTOPY
  !--- FORTAN-C interface ---
  TYPE(C_PTR), bind(C, name='_fringe_center_c') :: fringe_centerptr
  TYPE(C_PTR), bind(C,name='_wk_flow_center_c') :: wk_flow_centerptr
  TYPE(C_PTR), bind(C,name='_wk_pre_c') :: wk_preptr
#endif

  INTEGER :: n_data,n_data_tot,data_shift
	REAL, ALLOCATABLE :: write_mean_gbl(:,:,:,:)

  TYPE stats_t
     INTEGER :: i_data,m
     INTEGER, POINTER :: x(:),y(:),z(:)  !i,j,k of the grid node for each m component
     REAL, POINTER :: mean_xyz(:),covar_xyz(:),mean_xyzt(:),covar_xyzt(:) !stats
     REAL, POINTER :: wgt(:) !stats
     TYPE(stats_t), POINTER :: next
  END TYPE stats_t

  TYPE stats_group_t
     INTEGER :: group_id
     INTEGER :: n_data,n_data_tot,data_shift
     TYPE(stats_t), pointer :: stats_first
     TYPE(stats_group_t), POINTER :: next
  END TYPE stats_group_t

  TYPE(stats_group_t), pointer :: stats_group_first

  TYPE kalman_t
     INTEGER :: i_data,m
     INTEGER, POINTER :: x(:),y(:),z(:)  !i,j,k of the grid node for each m component
     REAL, POINTER :: mean(:,:),covar(:,:) !kalman
     REAL, POINTER :: muf(:),pf(:,:),obs_data(:),obs_covar(:,:),obs_oper(:,:),K(:,:)
     TYPE(kalman_t), POINTER :: next
  END TYPE kalman_t

  TYPE(kalman_t), pointer :: kalman_first

END MODULE usr_vars
