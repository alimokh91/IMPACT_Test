!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014 										                                                                         *
!* added ISO C bindings, Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch          *
!* February 2016                                                                                             *
!*************************************************************************************************************

!*************************************************************************************************************
!* Note: (bbecsek 110216)                                                                                    *
!*   New variables of PROTECTED type had to be added in order to make them accessible from C.                *
!*   PARAMETER type values cannot be accessed from C.                                                        *
!*   Other variables are not assigned ISO_C_BINDINGS because we are only using GCC compilers, and thus       *
!*   we do not need compiler independent access.                                                             *
!*   Arrays are accessed through pointers from the binding.                                                  *
!*************************************************************************************************************

!> @file mod_vars.f90
!! File containing variable declarations

!> All the variables are declared here. It uses the module mod_vars
MODULE mod_vars
  
  USE mod_dims
  USE ISO_C_BINDING !bbecsek

  IMPLICIT NONE
  

#ifdef ALLOC
  
  !===========================================================================================================
  !=== Raemliche Dimensionen =================================================================================
  !===========================================================================================================
  ! 2D wird ueber M3==2 eingeschaltet
  INTEGER                ::  dimens
  
  
  !===========================================================================================================
  !=== Domain- und Blockspezifikationen ======================================================================
  !===========================================================================================================
  !--- zulaessige Blockgroessen ------------------------------------------------------------------------------
  ! n  x 1   2   4   8   16   32   64   128   256  512   ...
  !--------------------------------------------------------------
  ! n2 = 2,  4,  8, 16,  32,  64, 128,  256,  512, 1024, ...
  ! n3 = 3,  6, 12, 24,  48,  96, 192,  384,  768, 1536, ...
  ! n5 = 5, 10, 20, 40,  80, 160, 320,  640, 1280, 2560, ...
  ! n7 = 7, 14, 28, 56, 112, 224, 448,  896, 1792, 3584, ...
  ! n9 = 9, 18, 36, 72, 144, 288, 576, 1152, 2304, 4608, ...
  ! ...
  INTEGER                ::  N1 
  INTEGER                ::  N2
  INTEGER                ::  N3
  
  !--- Anzahl grobe Gitter (Multigrid) -----------------------------------------------------------------------
  INTEGER, PARAMETER     ::  n_grids_max = 15
  INTEGER                ::  n_grids, n_grids_limit
#ifndef FTOPY
  !--- C-FORTRAN interface ---
  INTEGER, PROTECTED, bind(C, name='_n_grids_max_c') :: n_grids_max_c = n_grids_max
#endif

  
#ifndef FTOPY
  !--- Dimensionen -------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  dim_ncb1c = SIZE(ncb1c)
  INTEGER, PARAMETER     ::  dim_ncb1g = SIZE(ncb1g)
  INTEGER, PARAMETER     ::  dim_ncb1d = SIZE(ncb1d)
  
  INTEGER, PARAMETER     ::  dim_ncb2c = SIZE(ncb2c)
  INTEGER, PARAMETER     ::  dim_ncb2g = SIZE(ncb2g)
  INTEGER, PARAMETER     ::  dim_ncb2d = SIZE(ncb2d)
  
  INTEGER, PARAMETER     ::  dim_ncb3c = SIZE(ncb3c)
  INTEGER, PARAMETER     ::  dim_ncb3g = SIZE(ncb3g)
  INTEGER, PARAMETER     ::  dim_ncb3d = SIZE(ncb3d)

  !--- C-Fortran interface ---
  INTEGER, PROTECTED, bind(C, name='_dim_ncb1c_c') :: dim_ncb1c_c = dim_ncb1c
  INTEGER, PROTECTED, bind(C, name='_dim_ncb1g_c') :: dim_ncb1g_c = dim_ncb1g
  INTEGER, PROTECTED, bind(C, name='_dim_ncb1d_c') :: dim_ncb1d_c = dim_ncb1d
      
  INTEGER, PROTECTED, bind(C, name='_dim_ncb2c_c') :: dim_ncb2c_c = dim_ncb2c
  INTEGER, PROTECTED, bind(C, name='_dim_ncb2g_c') :: dim_ncb2g_c = dim_ncb2g
  INTEGER, PROTECTED, bind(C, name='_dim_ncb2d_c') :: dim_ncb2d_c = dim_ncb2d
         
  INTEGER, PROTECTED, bind(C, name='_dim_ncb3c_c') :: dim_ncb3c_c = dim_ncb3c
  INTEGER, PROTECTED, bind(C, name='_dim_ncb3g_c') :: dim_ncb3g_c = dim_ncb3d
  INTEGER, PROTECTED, bind(C, name='_dim_ncb3d_c') :: dim_ncb3d_c = dim_ncb3g

  !--- Anzahl Stencil-Koeffizienten (Feld) -------------------------------------------------------------------
  ! Anzahl der Koeffizienten im Feld (zentrale Differenzen angenommen):
  INTEGER, PARAMETER     ::  nc1c = ncb1c(dim_ncb1c)
  INTEGER, PARAMETER     ::  nc1s = ncb1g(dim_ncb1g)
  
  INTEGER, PARAMETER     ::  nc2c = ncb2c(dim_ncb2c)
  INTEGER, PARAMETER     ::  nc2s = ncb2g(dim_ncb2g)
  
  INTEGER, PARAMETER     ::  nc3c = ncb3c(dim_ncb3c)
  INTEGER, PARAMETER     ::  nc3s = ncb3g(dim_ncb3g)

  !--- C-Fortran interface ---
  INTEGER, PROTECTED, bind(C, name='_nc1c_c') :: nc1c_c = nc1c
  INTEGER, PROTECTED, bind(C, name='_nc1s_c') :: nc1s_c = nc1s
  
  INTEGER, PROTECTED, bind(C, name='_nc2c_c') :: nc2c_c = nc2c
  INTEGER, PROTECTED, bind(C, name='_nc2s_c') :: nc2s_c = nc2s
  
  INTEGER, PROTECTED, bind(C, name='_nc3c_c') :: nc3c_c = nc3c
  INTEGER, PROTECTED, bind(C, name='_nc3s_c') :: nc3s_c = nc3s

  !===========================================================================================================
  !=== Intervallgrenzen der Differenzen-Koeffizienten-Arrays =================================================
  !===========================================================================================================
  !--- zentral -----------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  b1U = nc1s/2
  INTEGER, PARAMETER     ::  b2U = nc2s/2
  INTEGER, PARAMETER     ::  b3U = nc3s/2
  
  INTEGER, PARAMETER     ::  b1L = -b1U
  INTEGER, PARAMETER     ::  b2L = -b2U
  INTEGER, PARAMETER     ::  b3L = -b3U

  !--- C-FORTRAN interface ---
  INTEGER, PROTECTED, bind(C, name='_b1u_c') :: b1u_c = b1U
  INTEGER, PROTECTED, bind(C, name='_b2u_c') :: b2u_c = b2U
  INTEGER, PROTECTED, bind(C, name='_b3u_c') :: b3u_c = b3U

  INTEGER, PROTECTED, bind(C, name='_b1l_c') :: b1l_c = b1L
  INTEGER, PROTECTED, bind(C, name='_b2l_c') :: b2l_c = b2L
  INTEGER, PROTECTED, bind(C, name='_b3l_c') :: b3l_c = b3L
 
  !--- upwind (nicht-linear) ---------------------------------------------------------------------------------
  ! (aktuell wird nicht zwischen auf- und abwärtsgerichteten Stencils unterschieden, um auch am Rand arbeiten
  !  zu können, wo KEINE upwind-Differenzierung verwendet wird)
  INTEGER, PARAMETER     ::  n1L = b1L
  INTEGER, PARAMETER     ::  n2L = b2L
  INTEGER, PARAMETER     ::  n3L = b3L
  
  INTEGER, PARAMETER     ::  n1U = b1U
  INTEGER, PARAMETER     ::  n2U = b2U
  INTEGER, PARAMETER     ::  n3U = b3U

  !--- C-FORTRAN interface ---
  INTEGER, PROTECTED, bind(C, name='_n1u_c') :: n1u_c = n1U
  INTEGER, PROTECTED, bind(C, name='_n2u_c') :: n2u_c = n2U
  INTEGER, PROTECTED, bind(C, name='_n3u_c') :: n3u_c = n3U

  INTEGER, PROTECTED, bind(C, name='_n1l_c') :: n1l_c = n1L
  INTEGER, PROTECTED, bind(C, name='_n2l_c') :: n2l_c = n2L
  INTEGER, PROTECTED, bind(C, name='_n3l_c') :: n3l_c = n3L
 
  !--- Divergenz ---------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  d1L = b1L
  INTEGER, PARAMETER     ::  d2L = b2L
  INTEGER, PARAMETER     ::  d3L = b3L
  
  INTEGER, PARAMETER     ::  d1U = b1U-1
  INTEGER, PARAMETER     ::  d2U = b2U-1
  INTEGER, PARAMETER     ::  d3U = b3U-1

  !--- C-FORTRAN interface ---
  INTEGER, PROTECTED, bind(C, name='_d1u_c') :: d1u_c = d1U
  INTEGER, PROTECTED, bind(C, name='_d2u_c') :: d2u_c = d2U
  INTEGER, PROTECTED, bind(C, name='_d3u_c') :: d3u_c = d3U

  INTEGER, PROTECTED, bind(C, name='_d1l_c') :: d1l_c = d1L
  INTEGER, PROTECTED, bind(C, name='_d2l_c') :: d2l_c = d2L
  INTEGER, PROTECTED, bind(C, name='_d3l_c') :: d3l_c = d3L
 
  !--- Gradient ----------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  g1L = b1L+1
  INTEGER, PARAMETER     ::  g2L = b2L+1
  INTEGER, PARAMETER     ::  g3L = b3L+1
  
  INTEGER, PARAMETER     ::  g1U = b1U
  INTEGER, PARAMETER     ::  g2U = b2U
  INTEGER, PARAMETER     ::  g3U = b3U

  !--- C-FORTRAN interface ---
  INTEGER, PROTECTED, bind(C, name='_g1u_c') :: g1u_c = g1U
  INTEGER, PROTECTED, bind(C, name='_g2u_c') :: g2u_c = g2U
  INTEGER, PROTECTED, bind(C, name='_g3u_c') :: g3u_c = g3U

  INTEGER, PROTECTED, bind(C, name='_g1l_c') :: g1l_c = g1L
  INTEGER, PROTECTED, bind(C, name='_g2l_c') :: g2l_c = g2L
  INTEGER, PROTECTED, bind(C, name='_g3l_c') :: g3l_c = g3L
 
#else
  INTEGER, PARAMETER     ::  dim_ncb1c = 4
  INTEGER, PARAMETER     ::  dim_ncb1g = 4
  INTEGER, PARAMETER     ::  dim_ncb1d = 4
  
  INTEGER, PARAMETER     ::  dim_ncb2c = 4
  INTEGER, PARAMETER     ::  dim_ncb2g = 4
  INTEGER, PARAMETER     ::  dim_ncb2d = 4
  
  INTEGER, PARAMETER     ::  dim_ncb3c = 4
  INTEGER, PARAMETER     ::  dim_ncb3g = 4
  INTEGER, PARAMETER     ::  dim_ncb3d = 4

  !===========================================================================================================
  !=== Intervallgrenzen der Differenzen-Koeffizienten-Arrays =================================================
  !===========================================================================================================
  !--- zentral -----------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  b1U = 3
  INTEGER, PARAMETER     ::  b2U = 3
  INTEGER, PARAMETER     ::  b3U = 3
  
  INTEGER, PARAMETER     ::  b1L = -3
  INTEGER, PARAMETER     ::  b2L = -3
  INTEGER, PARAMETER     ::  b3L = -3

  !--- upwind (nicht-linear) ---------------------------------------------------------------------------------
  ! (aktuell wird nicht zwischen auf- und abwärtsgerichteten Stencils unterschieden, um auch am Rand arbeiten
  !  zu können, wo KEINE upwind-Differenzierung verwendet wird)
  INTEGER, PARAMETER     ::  n1L = -3
  INTEGER, PARAMETER     ::  n2L = -3
  INTEGER, PARAMETER     ::  n3L = -3
  
  INTEGER, PARAMETER     ::  n1U = 3
  INTEGER, PARAMETER     ::  n2U = 3
  INTEGER, PARAMETER     ::  n3U = 3

  !--- Divergenz ---------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  d1L = -3
  INTEGER, PARAMETER     ::  d2L = -3
  INTEGER, PARAMETER     ::  d3L = -3
  
  INTEGER, PARAMETER     ::  d1U = 2
  INTEGER, PARAMETER     ::  d2U = 2
  INTEGER, PARAMETER     ::  d3U = 2

  !--- Gradient ----------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  g1L = -2
  INTEGER, PARAMETER     ::  g2L = -2
  INTEGER, PARAMETER     ::  g3L = -2
  
  INTEGER, PARAMETER     ::  g1U = 3
  INTEGER, PARAMETER     ::  g2U = 3
  INTEGER, PARAMETER     ::  g3U = 3

#endif
 
  
  !===========================================================================================================
  !=== Differenzen-Koeffizienten-Arrays ======================================================================
  !===========================================================================================================
  !--- 1. Ableitung (zentral) --------------------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  cp1  (:,:)
  REAL   , ALLOCATABLE   ::  cp2  (:,:)
  REAL   , ALLOCATABLE   ::  cp3  (:,:)
  
  REAL   , ALLOCATABLE   ::  cu1  (:,:)
  REAL   , ALLOCATABLE   ::  cv2  (:,:)
  REAL   , ALLOCATABLE   ::  cw3  (:,:)
  
  !--- 1. Ableitung (upwind) ---------------------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  cNp1D(:,:)
  REAL   , ALLOCATABLE   ::  cNp2D(:,:)
  REAL   , ALLOCATABLE   ::  cNp3D(:,:)
  
  REAL   , ALLOCATABLE   ::  cNp1U(:,:)
  REAL   , ALLOCATABLE   ::  cNp2U(:,:)
  REAL   , ALLOCATABLE   ::  cNp3U(:,:)
  
  REAL   , ALLOCATABLE   ::  cNu1D(:,:)
  REAL   , ALLOCATABLE   ::  cNv2D(:,:)
  REAL   , ALLOCATABLE   ::  cNw3D(:,:)
  
  REAL   , ALLOCATABLE   ::  cNu1U(:,:)
  REAL   , ALLOCATABLE   ::  cNv2U(:,:)
  REAL   , ALLOCATABLE   ::  cNw3U(:,:)
  
  !--- Divergenz ---------------------------------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  cDu1 (:,:)
  REAL   , ALLOCATABLE   ::  cDv2 (:,:)
  REAL   , ALLOCATABLE   ::  cDw3 (:,:)
  
  !--- Divergenz (transponiert) ------------------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  cDu1T(:,:)
  REAL   , ALLOCATABLE   ::  cDv2T(:,:)
  REAL   , ALLOCATABLE   ::  cDw3T(:,:)
  
  !--- Gradient ----------------------------------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  cGp1 (:,:)
  REAL   , ALLOCATABLE   ::  cGp2 (:,:)
  REAL   , ALLOCATABLE   ::  cGp3 (:,:)
  
  !--- Gradient (transponiert) -------------------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  cGp1T(:,:)
  REAL   , ALLOCATABLE   ::  cGp2T(:,:)
  REAL   , ALLOCATABLE   ::  cGp3T(:,:)
  
  !--- 2. Ableitung (zentral) --------------------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  cp11 (:,:)
  REAL   , ALLOCATABLE   ::  cp22 (:,:)
  REAL   , ALLOCATABLE   ::  cp33 (:,:)
  
  REAL   , ALLOCATABLE   ::  cu11 (:,:)
  REAL   , ALLOCATABLE   ::  cv22 (:,:)
  REAL   , ALLOCATABLE   ::  cw33 (:,:)
  
  !--- Interpolation ----------------------------------------------------------------------------------------- 
  REAL   , ALLOCATABLE   ::  cIpu(:,:)
  REAL   , ALLOCATABLE   ::  cIpv(:,:)
  REAL   , ALLOCATABLE   ::  cIpw(:,:)
  
  REAL   , ALLOCATABLE   ::  cIup(:,:)
  REAL   , ALLOCATABLE   ::  cIvp(:,:)
  REAL   , ALLOCATABLE   ::  cIwp(:,:)
  
  !--- Filter ------------------------------------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  cFp1(:,:)
  REAL   , ALLOCATABLE   ::  cFp2(:,:)
  REAL   , ALLOCATABLE   ::  cFp3(:,:)
  
  REAL   , ALLOCATABLE   ::  cFu1(:,:)
  REAL   , ALLOCATABLE   ::  cFv2(:,:)
  REAL   , ALLOCATABLE   ::  cFw3(:,:)
  
  !--- Integrator (nur für Druckgitter) ----------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  cInt1(:,:)
  REAL   , ALLOCATABLE   ::  cInt2(:,:)
  REAL   , ALLOCATABLE   ::  cInt3(:,:)  
  
  !--- 2. Ableitung (Multigrid) ------------------------------------------------------------------------------ 
  ! Anmerkung: Die Koeffizientensätze unterscheiden sich z.T. lediglich durch die Randbedingungen.
  REAL   , ALLOCATABLE   ::  cp11R(:,:,:)
  REAL   , ALLOCATABLE   ::  cp22R(:,:,:)
  REAL   , ALLOCATABLE   ::  cp33R(:,:,:)
  
  REAL   , ALLOCATABLE   ::  cu11R(:,:,:)
  REAL   , ALLOCATABLE   ::  cv22R(:,:,:)
  REAL   , ALLOCATABLE   ::  cw33R(:,:,:)
  
  REAL   , ALLOCATABLE   ::  cdg1 (:,:,:)
  REAL   , ALLOCATABLE   ::  cdg2 (:,:,:)
  REAL   , ALLOCATABLE   ::  cdg3 (:,:,:)
  
  !--- Interpolation (Multigrid) ----------------------------------------------------------------------------- 
  REAL   , ALLOCATABLE   ::  cI1(:,:,:)
  REAL   , ALLOCATABLE   ::  cI2(:,:,:)
  REAL   , ALLOCATABLE   ::  cI3(:,:,:)
  
  REAL   , ALLOCATABLE   ::  cIH1(:,:)
  REAL   , ALLOCATABLE   ::  cIH2(:,:)
  REAL   , ALLOCATABLE   ::  cIH3(:,:)
  
  !--- Restriktion (Multigrid) ------------------------------------------------------------------------------- 
  REAL   , ALLOCATABLE   ::  cR1 (:,:,:)
  REAL   , ALLOCATABLE   ::  cR2 (:,:,:)
  REAL   , ALLOCATABLE   ::  cR3 (:,:,:)
  
  REAL   , ALLOCATABLE   ::  cRest1(:,:,:) ! TEST!!!
  REAL   , ALLOCATABLE   ::  cRest2(:,:,:)
  REAL   , ALLOCATABLE   ::  cRest3(:,:,:)
  
  REAL   , ALLOCATABLE   ::  cRH1(:,:)
  REAL   , ALLOCATABLE   ::  cRH2(:,:)
  REAL   , ALLOCATABLE   ::  cRH3(:,:)
  
  
  !===========================================================================================================
  !=== Gitterspezifikationen =================================================================================
  !===========================================================================================================
  !--- physiklische Koordinaten (global) ---------------------------------------------------------------------
  REAL   , ALLOCATABLE, TARGET   ::  y1p(:), y1u(:)
  REAL   , ALLOCATABLE, TARGET   ::  y2p(:), y2v(:)
  REAL   , ALLOCATABLE, TARGET   ::  y3p(:), y3w(:)

#ifndef FTOPY
  !--- C-FORTRAN interface ---
  TYPE(C_PTR), bind(C, name='_y1p_c') :: y1pptr
  TYPE(C_PTR), bind(C, name='_y2p_c') :: y2pptr
  TYPE(C_PTR), bind(C, name='_y3p_c') :: y3pptr

  TYPE(C_PTR), bind(C, name='_y1u_c') :: y1uptr
  TYPE(C_PTR), bind(C, name='_y2v_c') :: y2vptr
  TYPE(C_PTR), bind(C, name='_y3w_c') :: y3wptr
#endif

  !--- physiklische Koordinaten (Block) ----------------------------------------------------------------------
  REAL   , ALLOCATABLE, TARGET   ::  x1p(:), x1u(:)
  REAL   , ALLOCATABLE, TARGET   ::  x2p(:), x2v(:)
  REAL   , ALLOCATABLE, TARGET   ::  x3p(:), x3w(:)

#ifndef FTOPY
  !--- C-FORTRAN interface --- 
  TYPE(C_PTR), bind(C, name='_x1p_c') :: x1pptr
  TYPE(C_PTR), bind(C, name='_x2p_c') :: x2pptr
  TYPE(C_PTR), bind(C, name='_x3p_c') :: x3pptr

  TYPE(C_PTR), bind(C, name='_x1u_c') :: x1uptr
  TYPE(C_PTR), bind(C, name='_x2v_c') :: x2vptr
  TYPE(C_PTR), bind(C, name='_x3w_c') :: x3wptr
#endif

  !--- physiklische Koordinaten (Block, Multigrid) -----------------------------------------------------------
  REAL   , ALLOCATABLE   ::  x1pR(:,:), x1uR(:,:)
  REAL   , ALLOCATABLE   ::  x2pR(:,:), x2vR(:,:)
  REAL   , ALLOCATABLE   ::  x3pR(:,:), x3wR(:,:)
  
  !--- Gitterweiten (global) ---------------------------------------------------------------------------------
  REAL   , ALLOCATABLE, TARGET   ::  dy1p(:), dy1u(:)
  REAL   , ALLOCATABLE, TARGET   ::  dy2p(:), dy2v(:)
  REAL   , ALLOCATABLE, TARGET   ::  dy3p(:), dy3w(:)

#ifndef FTOPY
  !--- C-FORTRAN interface --- 
  TYPE(C_PTR), bind(C, name='_dy1p_c') :: dy1pptr
  TYPE(C_PTR), bind(C, name='_dy2p_c') :: dy2pptr
  TYPE(C_PTR), bind(C, name='_dy3p_c') :: dy3pptr

  TYPE(C_PTR), bind(C, name='_dy1u_c') :: dy1uptr
  TYPE(C_PTR), bind(C, name='_dy2v_c') :: dy2vptr
  TYPE(C_PTR), bind(C, name='_dy3w_c') :: dy3wptr
#endif

  !--- Gitterweiten (Block) ----------------------------------------------------------------------------------
  REAL   , ALLOCATABLE, TARGET   ::  dx1p(:), dx1u(:)
  REAL   , ALLOCATABLE, TARGET   ::  dx2p(:), dx2v(:)
  REAL   , ALLOCATABLE, TARGET   ::  dx3p(:), dx3w(:)
  
  REAL   , ALLOCATABLE   ::  dx1DM(:), dx1pM(:), ddx1pM(:)
  REAL   , ALLOCATABLE   ::  dx2DM(:), dx2pM(:), ddx2pM(:)
  REAL   , ALLOCATABLE   ::  dx3DM(:), dx3pM(:), ddx3pM(:)
  
  REAL   , ALLOCATABLE   ::  dx1GM(:), dx1uM(:), ddx1uM(:)
  REAL   , ALLOCATABLE   ::  dx2GM(:), dx2vM(:), ddx2vM(:)
  REAL   , ALLOCATABLE   ::  dx3GM(:), dx3wM(:), ddx3wM(:)  

#ifndef FTOPY
  !--- C-FORTRAN interface --- 
  TYPE(C_PTR), bind(C, name='_dx1p_c') :: dx1pptr
  TYPE(C_PTR), bind(C, name='_dx2p_c') :: dx2pptr
  TYPE(C_PTR), bind(C, name='_dx3p_c') :: dx3pptr

  TYPE(C_PTR), bind(C, name='_dx1u_c') :: dx1uptr
  TYPE(C_PTR), bind(C, name='_dx2v_c') :: dx2vptr
  TYPE(C_PTR), bind(C, name='_dx3w_c') :: dx3wptr
#endif

  !===========================================================================================================
  !=== Arbeitsfelder =========================================================================================
  !===========================================================================================================
  !--- Geschwindigkeiten -------------------------------------------------------------------------------------
  REAL   , ALLOCATABLE, TARGET   ::  vel(:,:,:,:)
  !--- nicht-linearer Term -----------------------------------------------------------------------------------
  REAL   , ALLOCATABLE, TARGET   ::  nl (:,:,:,:)
  !--- Recht-Hand-Seite --------------------------------------------------------------------------------------
  REAL   , ALLOCATABLE, TARGET   ::  rhs(:,:,:,:)
  !--- Druck -------------------------------------------------------------------------------------------------
  REAL   , ALLOCATABLE, TARGET   ::  pre(:,:,:)

#ifndef FTOPY
  !--- C-FORTRAN interface ---
  TYPE(C_PTR), bind(C, name='_vel_c') :: velptr
  TYPE(C_PTR), bind(C, name='_nl_c') :: nlptr
  TYPE(C_PTR), bind(C, name='_rhs_c') :: rhsptr
  TYPE(C_PTR), bind(C, name='_pre_c') :: preptr
#endif
  
  !--- Ausfluss-RB (Geschwindigkeitsfeld) --------------------------------------------------------------------
  ! Da die RHS für die Konzentrationsfelder nicht über die Runge-Kutta-Zwischenschritte hinweg gespeichert 
  ! werden, müssen mindestens die zugehörigen Randbedingungen gespeichert werden.
  REAL   , ALLOCATABLE   ::  bc11(:,:,:), nlbc11(:,:,:)
  REAL   , ALLOCATABLE   ::  bc12(:,:,:), nlbc12(:,:,:)
  REAL   , ALLOCATABLE   ::  bc13(:,:,:), nlbc13(:,:,:)
  
  REAL   , ALLOCATABLE   ::  bc21(:,:,:), nlbc21(:,:,:)
  REAL   , ALLOCATABLE   ::  bc22(:,:,:), nlbc22(:,:,:)
  REAL   , ALLOCATABLE   ::  bc23(:,:,:), nlbc23(:,:,:)
  
  REAL   , ALLOCATABLE   ::  bc31(:,:,:), nlbc31(:,:,:)
  REAL   , ALLOCATABLE   ::  bc32(:,:,:), nlbc32(:,:,:)
  REAL   , ALLOCATABLE   ::  bc33(:,:,:), nlbc33(:,:,:)
  
  REAL   , ALLOCATABLE   ::  drift1(:,:,:)
  REAL   , ALLOCATABLE   ::  drift2(:,:,:)
  REAL   , ALLOCATABLE   ::  drift3(:,:,:)
  !--- Residuum ----------------------------------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  res (:,:,:)
  
  !--- Druckgradient (eine Komponente) -----------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  gpre(:,:,:)
  
  !--- Gewichte für Divergenzfreiheit ------------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  weight(:,:,:)
  
   !--- Null-Raum-Vektor --------------------------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  psi    (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_vel(:,:,:,:)
  
  REAL   , ALLOCATABLE   ::  psi_rel1 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel2 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel3 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel4 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel5 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel6 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel7 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel8 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel9 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel10(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel11(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel12(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel13(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel14(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel15(:,:,:)
  
  REAL   , ALLOCATABLE   ::  th11(:,:,:)
  REAL   , ALLOCATABLE   ::  th12(:,:,:)
  REAL   , ALLOCATABLE   ::  th13(:,:,:)
  
  REAL   , ALLOCATABLE   ::  th21(:,:,:)
  REAL   , ALLOCATABLE   ::  th22(:,:,:)
  REAL   , ALLOCATABLE   ::  th23(:,:,:)
  
  REAL   , ALLOCATABLE   ::  th31(:,:,:)
  REAL   , ALLOCATABLE   ::  th32(:,:,:)
  REAL   , ALLOCATABLE   ::  th33(:,:,:)
  
  !--- Multigrid ---------------------------------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  vec1C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec2A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec2B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec2C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec3A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec3B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec3C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec4A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec4B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec4C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec5A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec5B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec5C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec6A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec6B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec6C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec7A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec7B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec7C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec8A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec8B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec8C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec9A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec9B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec9C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec10A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec10B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec10C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec11A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec11B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec11C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec12A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec12B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec12C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec13A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec13B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec13C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec14A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec14B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec14C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec15A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec15B(:,:,:)
  
  
  !--- BiCGstab / Richardson ---------------------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  pp(:,:,:)
  REAL   , ALLOCATABLE   ::  Ap(:,:,:)
  REAL   , ALLOCATABLE   ::  rr(:,:,:)
  REAL   , ALLOCATABLE   ::  rh(:,:,:)
  REAL   , ALLOCATABLE   ::  Ar(:,:,:)
  REAL   , ALLOCATABLE   ::  z1(:,:,:)
  REAL   , ALLOCATABLE   ::  z2(:,:,:)
  
  !--- product_div_grad --------------------------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  dig(:,:,:)
  
  !--- Hilfsfelder (Druckiteration) --------------------------------------------------------------------------
  ! - rhs wird auch in test_moment nochmals verwendet und kann daher in outer_iteration nicht belegt werden!
  ! - wird auch fuer interpolierte Geschwindigkeiten in rhs_NS und rhs_conc verwendet.
  REAL   , ALLOCATABLE, TARGET   ::  work1(:,:,:)
  REAL   , ALLOCATABLE, TARGET   ::  work2(:,:,:)
  REAL   , ALLOCATABLE, TARGET   ::  work3(:,:,:)

#ifndef FTOPY
  !--- C-FORTRAN interface ---
  TYPE(C_PTR), bind(C, name='_work1_c') :: work1ptr
  TYPE(C_PTR), bind(C, name='_work2_c') :: work2ptr
  TYPE(C_PTR), bind(C, name='_work3_c') :: work3ptr
#endif
  
  !--- Linienrelaxation --------------------------------------------------------------------------------------
  REAL   , ALLOCATABLE   ::  vec1(:), dia1(:), SOR1(:), band1(:,:) ! TEST!!! siehe unten ...
  REAL   , ALLOCATABLE   ::  vec2(:), dia2(:), SOR2(:), band2(:,:)
  REAL   , ALLOCATABLE   ::  vec3(:), dia3(:), SOR3(:), band3(:,:)
  
  
  !===========================================================================================================
  !=== Indizierung (Intervallgrenzen, Verschiebungen) ========================================================
  !===========================================================================================================
  !--- Block-Index -------------------------------------------------------------------------------------------
  !INTEGER               ::  iB, jB, kB ! TEST!!! iB(1:3,1:n_grids_max) hierher verschieben ...
  
  !--- Indexverschiebung (Block --> global) ------------------------------------------------------------------
  INTEGER                ::  iShift, jShift, kShift
  
  !--- Domaingrösse (Periodizität-bereinigt) -----------------------------------------------------------------
  INTEGER                ::  dim1, dim2, dim3
  
  !--- Druck / Konzentrationen (inklusive Rand) --------------------------------------------------------------
  INTEGER                ::  S1p, S2p, S3p
  INTEGER                ::  N1p, N2p, N3p
  
  !--- Geschwindigkeiten (inklusive Rand) --------------------------------------------------------------------
  INTEGER                ::  S11B, S21B, S31B
  INTEGER                ::  S12B, S22B, S32B
  INTEGER                ::  S13B, S23B, S33B
  
  INTEGER                ::  N11B, N21B, N31B
  INTEGER                ::  N12B, N22B, N32B
  INTEGER                ::  N13B, N23B, N33B
  
  !--- Geschwindigkeiten (exklusive Rand) --------------------------------------------------------------------
  INTEGER                ::  S11, S21, S31
  INTEGER                ::  S12, S22, S32
  INTEGER                ::  S13, S23, S33
  
  INTEGER                ::  N11, N21, N31
  INTEGER                ::  N12, N22, N32
  INTEGER                ::  N13, N23, N33
  
  !--- grobe Gitter (Multigrid, INklusive Rand) --------------------------------------------------------------
  INTEGER                ::  S1R, S2R, S3R
  INTEGER                ::  d1R, d2R, d3R
  
  !--- grobe Gitter (Multigrid, EXklusive Rand) --------------------------------------------------------------
  INTEGER                ::  S11R, S22R, S33R
  INTEGER                ::  d11R, d22R, d33R
  
  !--- Ueberlappungskonvention der Blöcke (Multigrid, siehe mod_setup) ---------------------------------------
  INTEGER, PARAMETER     ::  ls1 = -1
  INTEGER, PARAMETER     ::  ls2 = -1
  INTEGER, PARAMETER     ::  ls3 = -1

#ifndef FTOPY
  !--- C-FORTRAN interface ---
  INTEGER, PROTECTED, bind(C, name='_ls1_c') :: ls1_c = ls1
  INTEGER, PROTECTED, bind(C, name='_ls2_c') :: ls2_c = ls2
  INTEGER, PROTECTED, bind(C, name='_ls3_c') :: ls3_c = ls3
#endif

  !--- Austauschrichtung (Multigrid) -------------------------------------------------------------------------
  ! ex = -1: unten <--  oben
  ! ex =  0: unten <--> oben
  ! ex =  1: unten  --> oben
  INTEGER                ::  ex1, ex2, ex3
  
  
  !===========================================================================================================
  !=== Randbedingungen =======================================================================================
  !===========================================================================================================
  !                              _
  !    Symmetrie-RB:   BC = -2    |
  !    periodische RB: BC = -1    |- symmetrische, zentrale Stencils
  !    Nachbar-Block:  BC =  0   _|
  !    Dirichlet-RB:   BC =  1    |
  !    Neumann-RB:     BC =  2    |- schiefe, nicht-zentrale Stencils
  !    Robin-RB:       BC =  3   _|
  !
  !--- global ------------------------------------------------------------------------------------------------
  LOGICAL, TARGET        ::  outlet(1:3,1:2,1:3)
  
  INTEGER                ::  BC_1L_global, BC_1U_global
  INTEGER                ::  BC_2L_global, BC_2U_global
  INTEGER                ::  BC_3L_global, BC_3U_global
  
  !--- lokal (Block) -----------------------------------------------------------------------------------------
  INTEGER                ::  BC_1L, BC_1U
  INTEGER                ::  BC_2L, BC_2U
  INTEGER                ::  BC_3L, BC_3U

  !--- field properties --------------------------------------------------------------------------------------
  INTEGER                ::  n_gather(1:3,1:n_grids_max)
  INTEGER                ::  NN (1:3,1:n_grids_max)
  INTEGER                ::  NB (1:3,1:n_grids_max)
  INTEGER, TARGET        ::  iB (1:3,1:n_grids_max)
  INTEGER                ::  SNF(1:2,1:3,1:n_grids_max)
  INTEGER                ::  SNB(1:2,1:3,1:n_grids_max)
  INTEGER                ::  BC (1:2,1:3,1:n_grids_max)
  INTEGER                ::  ngb(1:2,1:3,1:n_grids_max)
  INTEGER                ::  comm1(1:n_grids_max), comm2(1:n_grids_max)
  INTEGER                ::  rankc2(1:n_grids_max)
  LOGICAL                ::  participate_yes(1:n_grids_max)
  INTEGER, ALLOCATABLE   ::  recvR(  :,:), recvI(  :,:)
  INTEGER, ALLOCATABLE   ::  dispR(  :,:), dispI(  :,:)
  INTEGER, ALLOCATABLE   ::  offsR(:,:,:), offsI(:,:,:)
  INTEGER, ALLOCATABLE   ::  sizsR(:,:,:), sizsI(:,:,:)

#ifndef FTOPY
  !--- C-FORTRAN interface ---
  TYPE(C_PTR), bind(C, name='_ib_c') :: ibptr
  TYPE(C_PTR), bind(C, name='_outlet_c') :: outletptr
#endif
  
  !===========================================================================================================
  !=== physikalische Parameter ===============================================================================
  !===========================================================================================================
  REAL                   ::  L1, L2, L3
  REAL                   ::  L1_half, L2_half, L3_half
  REAL                   ::  L1_amp, L2_amp, L3_amp
  REAL                   ::  Re
  LOGICAL                ::  acos_yes
  
  
  !===========================================================================================================
  !=== numerische Parameter ==================================================================================
  !===========================================================================================================
  !--- allgemein ---------------------------------------------------------------------------------------------
  REAL                   ::  CFL
  REAL                   ::  time, dtime, subtime, time_start, time_end, dtime_max, dtime0, dtime_old
  INTEGER                ::  timestep, timestep_old, substep, n_timesteps
  LOGICAL                ::  mapping_yes, upwind_yes
  LOGICAL                ::  Euler_yes, Stokes_yes, twostep_yes
  LOGICAL, PARAMETER     ::  filter_BC_yes = .TRUE. ! TEST!!!
  INTEGER                ::  timeint_mode, forcing_mode
  INTEGER                ::  bulkflow_dir
  INTEGER                ::  n_lp_vel    , n_hp_vel
  REAL                   ::  chi_vel
  REAL                   ::  vel_bulk ! TEST!!!
  INTEGER                ::  timestep_mode !kschlegel
  REAL                   ::  dtime_const !kschlegel
  REAL                   ::  tol_dtime_pid_controller !kschlegel
  LOGICAL                ::  refine_2dir_yes !kschlegel
  REAL                   ::  refine_2dir_xi !kschlegel

#ifndef FTOPY
  !--- C-FORTRAN interface ---
  LOGICAL, PROTECTED, bind(C, name='_filter_BC_yes_c') :: filter_bc_yes_c = filter_BC_yes
#endif
  
  !--- Runge-Kutta-Koeffizienten -----------------------------------------------------------------------------
  REAL   , PARAMETER     ::  aRK(1:3) = (/8./15.,  5./12., 3./ 4./)
  REAL   , PARAMETER     ::  bRK(1:3) = (/  0.  ,-17./60.,-5./12./)
  INTEGER, PARAMETER     ::  RK_steps = 3

#ifndef FTOPY
  !--- C-FORTRAN interface ---
  REAL, PROTECTED, TARGET  ::  ark_c(1:3) = ark(1:3)
  REAL, PROTECTED, TARGET  ::  brk_c(1:3) = brk(1:3)
  TYPE(C_PTR), bind(C, name='_ark_c') :: arkptr
  TYPE(C_PTR), bind(C, name='_brk_c') :: brkptr
  INTEGER, PROTECTED, bind(C, name='_rk_steps_c') :: rk_steps_c = RK_steps
#endif

  !--- look-up table fuer Stabilitaetsgebiet der Zeitintegration (angle = pi/2,pi) ---------------------------
  REAL   , PARAMETER     ::  stabilitylimit(0:40) = (/1.732050813, 1.943689093, 2.089210537, 2.201001743,  &
                                      &               2.290031261, 2.361554127, 2.418567407, 2.462989697,  &
                                      &               2.496169963, 2.519146008, 2.532795254, 2.537935070,  &
                                      &               2.535397854, 2.526091466, 2.511046932, 2.491448818,  &
                                      &               2.468639045, 2.444084180, 2.419302172, 2.395757241,  &
                                      &               2.374745783, 2.357302135, 2.344145473, 2.335672458,  &
                                      &               2.331985072, 2.332936948, 2.338183901, 2.347230689,  &
                                      &               2.359471631, 2.374225928, 2.390769340, 2.408363261,  &
                                      &               2.426281290, 2.443832601, 2.460381269, 2.475360992,  &
                                      &               2.488285197, 2.498753090, 2.506452564, 2.511161051,  &
                                      &               2.512745327 /)


  INTEGER, PARAMETER     ::  n_stab = 41

#ifndef FTOPY
  !--- C-FORTRAN interface ---
  INTEGER, PROTECTED, bind(C, name='_n_stab_c') :: n_stab_c = n_stab
#endif

  !--- Helmholtz-Vorfaktoren ---------------------------------------------------------------------------------
  REAL                   ::  thetaL, multL
  
  !--- zeitliche Steuerung -----------------------------------------------------------------------------------
  INTEGER                ::  Int_dtime, Int_lev_pre
  
  INTEGER , TARGET       ::  stride_large(1:3), stride_med(1:3), stride_small(1:3)
  LOGICAL                ::  write_large, write_med, write_small
	INTEGER                ::  interval
	REAL                   ::  time_out_vect, dtime_out_vect
  REAL                   ::  time_out_scal, dtime_out_scal
  REAL                   ::  time_out_kalm, dtime_out_kalm
	REAL   , ALLOCATABLE   ::  dtime_kalm_phases(:)
  
  LOGICAL                ::  write_out_kalm, write_out_scal, write_out_vect
  LOGICAL                ::  new_dtime, finish_yes
    
  INTEGER                ::  write_count,write_stats_count,write_kalm_count
  INTEGER                ::  restart
  CHARACTER(LEN=3)       ::  restart_char
  
  INTEGER                ::  n_conc_old

#ifndef FTOPY
  !--- C-FORTRAN interface ---
  TYPE(C_PTR), bind(C, name='_stride_large_c') :: stride_largeptr
  TYPE(C_PTR), bind(C, name='_stride_med_c') :: stride_medptr
  TYPE(C_PTR), bind(C, name='_stride_small_c') :: stride_smallptr
#endif


  !===========================================================================================================
  !=== weitere Steuerungsoptionen ============================================================================
  !===========================================================================================================
  INTEGER                ::  task
  LOGICAL                ::  read_nullspace_yes
  LOGICAL                ::  nullspace_yes, nullspace_coarse_yes
  LOGICAL                ::  nullspace_ortho_yes
  
  LOGICAL                ::  write_stout_yes
  LOGICAL                ::  log_iteration_yes
  LOGICAL                ::  write_restart_yes
  LOGICAL                ::  write_lambda2_yes
  LOGICAL                ::  write_test_yes
  LOGICAL                ::  write_covariance_yes !for writing covariance into xdmf file    defined in config.txt
  INTEGER                ::  num_windows !define in usr_stats and for write covariance into xdmf
  INTEGER                ::  intervals   !define number of intervals in the output of a periodic flow
  INTEGER                ::  phase       !define number of interval  in the output of a periodic flow
  INTEGER, ALLOCATABLE   ::  repetition(:) !define array of #phases for repetitions counter

  !--- globale Laufindizes -----------------------------------------------------------------------------------
  INTEGER                ::  direction
  
  !--- explizite Behandlung von Ecken bei Dirichlet-Randbedingungen ------------------------------------------
  ! (Hintergrund: Der Druck ist an diesen Orten unbestimmt, so dass er dort künstlich zu Null gesetzt wird.)
  LOGICAL, PARAMETER     ::  corner_yes = .TRUE.

#ifndef FTOPY
  !--- C-FORTRAN interface ---
  LOGICAL, PROTECTED, bind(C, name='_corner_yes_c') :: corner_yes_c = corner_yes
#endif

  !--- Systemzeit --------------------------------------------------------------------------------------------
  INTEGER                ::  elatime, ctime(1:8)
  INTEGER                ::  day, hour, minu, sec, msec
  
  
  !===========================================================================================================
  !=== Iterationsparameter ===================================================================================
  !===========================================================================================================
  !--- Abbruchkriterium / Absolute Genauigkeit der Geschwindigkeiten -----------------------------------------
  REAL                   ::  epsU, epsU0
  
  !--- Glaetter ----------------------------------------------------------------------------------------------
  LOGICAL                ::  Jacobi_yes
  
  !--- max. Anzahl Iterationen -------------------------------------------------------------------------------
  INTEGER                ::  n_it_outer
  INTEGER                ::  n_it_Poisson
  INTEGER                ::  n_it_Helmh_vel
  
  !--- erwartete Konvergenzrate (äussere Iteration) ----------------------------------------------------------
  REAL   , TARGET        ::  precRatio0 (1:RK_steps)
  REAL   , TARGET        ::  precOffset0(1:RK_steps)
  REAL   , ALLOCATABLE   ::  precRatio  (:,:)
  REAL   , ALLOCATABLE   ::  precOffset (:,:)
  
  !--- Null-Initialisierung (äussere Iteration) --------------------------------------------------------------
  LOGICAL, TARGET                ::  init_pre(1:RK_steps), init_vel(1:RK_steps)!, init_conc(1:RK_steps)

#ifndef FTOPY
  !--- C-Fortran interface ---
  TYPE(C_PTR), bind(C, name='_init_pre_c') :: init_preptr
  TYPE(C_PTR), bind(C, name='_init_vel_c') :: init_velptr
  TYPE(C_PTR), bind(C, name='_precoffset0_c') :: precoffset0ptr
  TYPE(C_PTR), bind(C, name='_precratio0_c') :: precratio0ptr
#endif
  
  !--- Vorkonditionierung (Multigrid) ------------------------------------------------------------------------
  !--- Vorkonditionierung (Multigrid) ------------------------------------------------------------------------
  INTEGER                ::  precond_outer
  INTEGER                ::  precond_Poisson
  INTEGER                ::  precond_Helmh_vel
  !INTEGER                ::  precond_Helmh_conc
  
  !--- Anzahl Glättungen pro Gitterlevel (Multigrid) ---------------------------------------------------------
  INTEGER                ::  n_relax_down, n_relax_up, n_relax_bottom
  
  !--- implizite Richtungen bei Linienrelaxation (Multigrid) -------------------------------------------------
  INTEGER, TARGET        ::  impl_dir(1:3)

#ifndef FTOPY
  !--- C-Fortran interface ---
  TYPE(C_PTR), bind(C, name='_impl_dir_c') :: impl_dirptr
#endif

  !--- Anzahl Glättungen pro Gitterlevel (Multigrid) ---------------------------------------------------------
  LOGICAL                ::  weighting_yes
  
  
  !===========================================================================================================
  !=== Iterationsstatistik ===================================================================================
  !===========================================================================================================
  REAL                   ::  dtime_average
  REAL                   ::  max_div_init(1:2)
  INTEGER                ::  number_poisson
  
  !--- Zähler ------------------------------------------------------------------------------------------------
  INTEGER                ::  countO(1:RK_steps)
  INTEGER                ::  countP(1:RK_steps,1:2)
  INTEGER                ::  countH(1:RK_steps,1:3)
  
  !--- Konvergenzrate ----------------------------------------------------------------------------------------
  REAL                   ::  ratioO(1:RK_steps)
  REAL                   ::  ratioH(1:RK_steps,1:3)
  REAL                   ::  ratioP(1:RK_steps,1:2)
  
  
  !===========================================================================================================
  !=== MPI ===================================================================================================
  !===========================================================================================================
  !--- Kommunikatoren ----------------------------------------------------------------------------------------
  INTEGER                ::  COMM_CART
  
  INTEGER                ::  COMM_SLICE1, COMM_BAR1
  INTEGER                ::  COMM_SLICE2, COMM_BAR2
  INTEGER                ::  COMM_SLICE3, COMM_BAR3
  
  !--- Dimension und Position der Blöcke innerhalb der Kommunikatoren (Gitter-Indizes) -----------------------
  ! (für MPI_GATHERv, MPI_ALLGATHERv, vgl. iShift, jShift, kShift)
  INTEGER, ALLOCATABLE   ::  bar1_size(:), bar1_offset(:)
  INTEGER, ALLOCATABLE   ::  bar2_size(:), bar2_offset(:)
  INTEGER, ALLOCATABLE   ::  bar3_size(:), bar3_offset(:)
  
  !--- Ränge der Prozesse ------------------------------------------------------------------------------------
  INTEGER                ::  rank
  INTEGER                ::  rank_bar1, rank_slice1
  INTEGER                ::  rank_bar2, rank_slice2
  INTEGER                ::  rank_bar3, rank_slice3
  
  !--- Ränge der Nachbarprozesse (in kartesischem Gitter) ----------------------------------------------------
  INTEGER                ::  rank1L, rank1U
  INTEGER                ::  rank2L, rank2U
  INTEGER                ::  rank3L, rank3U
  
  !--- Error-Handle ------------------------------------------------------------------------------------------
  INTEGER                ::  merror
  
  !--- Request-Handles ---------------------------------------------------------------------------------------
  ! (müssen offenbar global angelegt werden) 
  INTEGER                ::  req1L, req1U
  INTEGER                ::  req2L, req2U
  INTEGER                ::  req3L, req3U
  
  !===========================================================================================================
  !=== HDF5 ==================================================================================================
  !===========================================================================================================
  INTEGER                ::  herror
  
  !*** bbecsek 2015: IMMERSED BOUNDARY PARAMETERS ************************************************************
  REAL   , ALLOCATABLE, TARGET   ::  fd(:,:,:,:)             !< force density
  REAL                   ::  mu_fluid                !< dynamic viscosity of blood
  REAL                   ::  rho_fluid               !< bulk density of blood

  REAL                   ::  L_ref                   !< reference length scale
  REAL                   ::  U_ref                   !< reference velocity

  LOGICAL                ::  write_force_yes         !< whether to output force density or not
  LOGICAL                ::  write_xdmf_yes          !< whether to output .xmf files
  LOGICAL                ::  scale_output_yes        !< whether to scale all values with U_ref, L_ref (dimensionalize)

  INTEGER, ALLOCATABLE   ::  block_sizes(:)          !< sizes (# of nodes in each block ) 0 -> NB1*NB2*NB3
  INTEGER, ALLOCATABLE   ::  block_indices(:)        !< indices (start & end) 0 -> 6*NB1*NB2*NB3
  INTEGER, ALLOCATABLE   ::  block_elem_sizes(:)     !< sizes (# of elements in each block ) 0 -> NB1*NB2*NB3

  REAL, ALLOCATABLE, TARGET      ::  vel_old(:,:,:,:)        !< velocity of previous timestep (used in Picard iterations)

  INTEGER :: n_data,n_data_tot,data_shift
  REAL, ALLOCATABLE :: mean_gbl(:,:,:,:,:)
  REAL, ALLOCATABLE :: covar_gbl(:,:,:,:,:)
  REAL, ALLOCATABLE :: write_mean (:,:,:,:)
  REAL, ALLOCATABLE :: write_fluct (:,:,:,:)
  REAL, ALLOCATABLE :: write_covar(:,:,:,:)
  REAL, ALLOCATABLE :: write_gain (:,:,:,:)

  TYPE stats_t
     INTEGER :: i_data,m
     INTEGER, POINTER :: x(:),y(:),z(:)  !i,j,k of the grid node for each m component
     REAL, POINTER :: mean_xyz(:),covar_xyz(:),mean_xyzt(:,:),covar_xyzt(:,:) !stats
     REAL, POINTER :: wgt(:) !stats
     TYPE(stats_t), POINTER :: next
  END TYPE stats_t

  TYPE stats_group_t
     INTEGER :: group_id,phase
     INTEGER :: n_data,n_data_tot,data_shift
     TYPE(stats_t), pointer :: stats_first
     TYPE(stats_group_t), POINTER :: next
  END TYPE stats_group_t

  TYPE(stats_group_t), pointer :: stats_group_first

  TYPE kalman_t
     INTEGER :: i_data,flg,m,phase
     INTEGER, POINTER :: x(:),y(:),z(:)  !i,j,k of the grid node for each m component
     REAL, POINTER :: mean(:,:,:),covar(:,:,:) ! phases, m, 3 or 6
     REAL, POINTER :: muf(:),pf(:,:),obs_data(:),obs_covar(:,:),obs_oper(:,:),K(:,:) !kalman
     TYPE(kalman_t), POINTER :: next
  END TYPE kalman_t

  TYPE(kalman_t), pointer :: kalman_first


#ifndef FTOPY
  !--- C-FORTRAN interface ---
  TYPE(C_PTR), bind(C, name='_fd_c') :: fdptr
  TYPE(C_PTR), bind(C, name='_vel_old_c') :: vel_oldptr
#endif

  !--- Load velocity initalcond from hdf5-file ---
  LOGICAL                ::  vel_initcond_file_yes




#else
  
  
  !===========================================================================================================
  !=== Raemliche Dimensionen =================================================================================
  !===========================================================================================================
  ! 2D wird ueber M3==2 eingeschaltet
  INTEGER                ::  dimens
  
  
  !===========================================================================================================
  !=== Domain- und Blockspezifikationen ======================================================================
  !===========================================================================================================
  !--- zulaessige Blockgroessen ------------------------------------------------------------------------------
  ! n  x 1   2   4   8   16   32   64   128   256  512   ...
  !--------------------------------------------------------------
  ! n2 = 2,  4,  8, 16,  32,  64, 128,  256,  512, 1024, ...
  ! n3 = 3,  6, 12, 24,  48,  96, 192,  384,  768, 1536, ...
  ! n5 = 5, 10, 20, 40,  80, 160, 320,  640, 1280, 2560, ...
  ! n7 = 7, 14, 28, 56, 112, 224, 448,  896, 1792, 3584, ...
  ! n9 = 9, 18, 36, 72, 144, 288, 576, 1152, 2304, 4608, ...
  ! ...
  INTEGER, PARAMETER     ::  N1 = 1+(M1-1)/NB1
  INTEGER, PARAMETER     ::  N2 = 1+(M2-1)/NB2
  INTEGER, PARAMETER     ::  N3 = 1+(M3-1)/NB3
  
  !--- Anzahl grobe Gitter (Multigrid) -----------------------------------------------------------------------
  INTEGER, PARAMETER     ::  n_grids_max = 15
  INTEGER                ::  n_grids, n_grids_limit
  
  !--- Dimensionen -------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  dim_ncb1c = SIZE(ncb1c)
  INTEGER, PARAMETER     ::  dim_ncb1g = SIZE(ncb1g)
  INTEGER, PARAMETER     ::  dim_ncb1d = SIZE(ncb1d)
  
  INTEGER, PARAMETER     ::  dim_ncb2c = SIZE(ncb2c)
  INTEGER, PARAMETER     ::  dim_ncb2g = SIZE(ncb2g)
  INTEGER, PARAMETER     ::  dim_ncb2d = SIZE(ncb2d)
                         
  INTEGER, PARAMETER     ::  dim_ncb3c = SIZE(ncb3c)
  INTEGER, PARAMETER     ::  dim_ncb3g = SIZE(ncb3g)
  INTEGER, PARAMETER     ::  dim_ncb3d = SIZE(ncb3d)
  
  !--- Anzahl Stencil-Koeffizienten (Feld)------------- ------------------------------------------------------
  ! Anzahl der Koeffizienten im Feld (zentrale Differenzen angenommen):
  INTEGER, PARAMETER     ::  nc1c = ncb1c(dim_ncb1c)
  INTEGER, PARAMETER     ::  nc1s = ncb1g(dim_ncb1g)
                         
  INTEGER, PARAMETER     ::  nc2c = ncb2c(dim_ncb2c)
  INTEGER, PARAMETER     ::  nc2s = ncb2g(dim_ncb2g)
                         
  INTEGER, PARAMETER     ::  nc3c = ncb3c(dim_ncb3c)
  INTEGER, PARAMETER     ::  nc3s = ncb3g(dim_ncb3g)
  
  
  !===========================================================================================================
  !=== Intervallgrenzen der Differenzen-Koeffizienten-Arrays =================================================
  !===========================================================================================================
  !--- zentral -----------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  b1U = nc1s/2
  INTEGER, PARAMETER     ::  b2U = nc2s/2
  INTEGER, PARAMETER     ::  b3U = nc3s/2
  
  INTEGER, PARAMETER     ::  b1L = -b1U
  INTEGER, PARAMETER     ::  b2L = -b2U
  INTEGER, PARAMETER     ::  b3L = -b3U
  
  !--- upwind (nicht-linear) ---------------------------------------------------------------------------------
  ! (aktuell wird nicht zwischen auf- und abwärtsgerichteten Stencils unterschieden, um auch am Rand arbeiten
  !  zu können, wo KEINE upwind-Differenzierung verwendet wird)
  INTEGER, PARAMETER     ::  n1L = b1L
  INTEGER, PARAMETER     ::  n2L = b2L
  INTEGER, PARAMETER     ::  n3L = b3L
  
  INTEGER, PARAMETER     ::  n1U = b1U
  INTEGER, PARAMETER     ::  n2U = b2U
  INTEGER, PARAMETER     ::  n3U = b3U
  
  !--- Divergenz ---------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  d1L = b1L
  INTEGER, PARAMETER     ::  d2L = b2L
  INTEGER, PARAMETER     ::  d3L = b3L
  
  INTEGER, PARAMETER     ::  d1U = b1U-1
  INTEGER, PARAMETER     ::  d2U = b2U-1
  INTEGER, PARAMETER     ::  d3U = b3U-1
  
  !--- Gradient ----------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  g1L = b1L+1
  INTEGER, PARAMETER     ::  g2L = b2L+1
  INTEGER, PARAMETER     ::  g3L = b3L+1
  
  INTEGER, PARAMETER     ::  g1U = b1U
  INTEGER, PARAMETER     ::  g2U = b2U
  INTEGER, PARAMETER     ::  g3U = b3U
  
  
  !===========================================================================================================
  !=== Differenzen-Koeffizienten-Arrays ======================================================================
  !===========================================================================================================
  !--- 1. Ableitung (zentral) --------------------------------------------------------------------------------
  REAL                   ::  cp1  (b1L:b1U,0:N1)
  REAL                   ::  cp2  (b2L:b2U,0:N2)
  REAL                   ::  cp3  (b3L:b3U,0:N3)
  
  REAL                   ::  cu1  (b1L:b1U,0:N1)
  REAL                   ::  cv2  (b2L:b2U,0:N2)
  REAL                   ::  cw3  (b3L:b3U,0:N3)
  
  !--- 1. Ableitung (upwind) ---------------------------------------------------------------------------------
  REAL                   ::  cNp1D(n1L:n1U,0:N1)
  REAL                   ::  cNp2D(n2L:n2U,0:N2)
  REAL                   ::  cNp3D(n3L:n3U,0:N3)
  
  REAL                   ::  cNp1U(n1L:n1U,0:N1)
  REAL                   ::  cNp2U(n2L:n2U,0:N2)
  REAL                   ::  cNp3U(n3L:n3U,0:N3)
  
  REAL                   ::  cNu1D(n1L:n1U,0:N1)
  REAL                   ::  cNv2D(n2L:n2U,0:N2)
  REAL                   ::  cNw3D(n3L:n3U,0:N3)
  
  REAL                   ::  cNu1U(n1L:n1U,0:N1)
  REAL                   ::  cNv2U(n2L:n2U,0:N2)
  REAL                   ::  cNw3U(n3L:n3U,0:N3)

  !--- Divergenz ---------------------------------------------------------------------------------------------
  REAL                   ::  cDu1 (d1L:d1U,0:N1)
  REAL                   ::  cDv2 (d2L:d2U,0:N2)
  REAL                   ::  cDw3 (d3L:d3U,0:N3)
  
  !--- Divergenz (transponiert) ------------------------------------------------------------------------------
  REAL                   ::  cDu1T(g1L:g1U,0:N1)
  REAL                   ::  cDv2T(g2L:g2U,0:N2)
  REAL                   ::  cDw3T(g3L:g3U,0:N3)
  
  !--- Gradient ----------------------------------------------------------------------------------------------
  REAL                   ::  cGp1 (g1L:g1U,0:N1)
  REAL                   ::  cGp2 (g2L:g2U,0:N2)
  REAL                   ::  cGp3 (g3L:g3U,0:N3)
  
  !--- Gradient (transponiert) -------------------------------------------------------------------------------
  REAL                   ::  cGp1T(d1L:d1U,0:N1)
  REAL                   ::  cGp2T(d2L:d2U,0:N2)
  REAL                   ::  cGp3T(d3L:d3U,0:N3)
  
  !--- 2. Ableitung (zentral) --------------------------------------------------------------------------------
  REAL                   ::  cp11 (b1L:b1U,0:N1)
  REAL                   ::  cp22 (b2L:b2U,0:N2)
  REAL                   ::  cp33 (b3L:b3U,0:N3)
  
  REAL                   ::  cu11 (b1L:b1U,0:N1)
  REAL                   ::  cv22 (b2L:b2U,0:N2)
  REAL                   ::  cw33 (b3L:b3U,0:N3)

  !--- Interpolation ----------------------------------------------------------------------------------------- 
  REAL                   ::  cIpu(g1L:g1U,0:N1)
  REAL                   ::  cIpv(g2L:g2U,0:N2)
  REAL                   ::  cIpw(g3L:g3U,0:N3)
  
  REAL                   ::  cIup(d1L:d1U,0:N1)
  REAL                   ::  cIvp(d2L:d2U,0:N2)
  REAL                   ::  cIwp(d3L:d3U,0:N3)

  !--- Filter ------------------------------------------------------------------------------------------------
  REAL                   ::  cFp1(b1L:b1U,0:N1)
  REAL                   ::  cFp2(b2L:b2U,0:N2)
  REAL                   ::  cFp3(b3L:b3U,0:N3)
  
  REAL                   ::  cFu1(b1L:b1U,0:N1)
  REAL                   ::  cFv2(b2L:b2U,0:N2)
  REAL                   ::  cFw3(b3L:b3U,0:N3)

  !--- Integrator (nur für Druckgitter) ----------------------------------------------------------------------
  REAL                   ::  cInt1(b1L:b1U,0:N1)
  REAL                   ::  cInt2(b2L:b2U,0:N2)
  REAL                   ::  cInt3(b3L:b3U,0:N3)
  
  !--- 2. Ableitung (Multigrid) ------------------------------------------------------------------------------ 
  ! Anmerkung: Die Koeffizientensätze unterscheiden sich z.T. lediglich durch die Randbedingungen.
  REAL                   ::  cp11R(-1:1,0:N1,1:n_grids_max)
  REAL                   ::  cp22R(-1:1,0:N2,1:n_grids_max)
  REAL                   ::  cp33R(-1:1,0:N3,1:n_grids_max)
  
  REAL                   ::  cu11R(-1:1,0:N1,1:n_grids_max)
  REAL                   ::  cv22R(-1:1,0:N2,1:n_grids_max)
  REAL                   ::  cw33R(-1:1,0:N3,1:n_grids_max)

  REAL                   ::  cdg1 (-1:1,1:N1,1:n_grids_max)
  REAL                   ::  cdg2 (-1:1,1:N2,1:n_grids_max)
  REAL                   ::  cdg3 (-1:1,1:N3,1:n_grids_max)
  
  !--- Interpolation (Multigrid) ----------------------------------------------------------------------------- 
  REAL                   ::  cI1(1:2,1:N1,1:n_grids_max)
  REAL                   ::  cI2(1:2,1:N2,1:n_grids_max)
  REAL                   ::  cI3(1:2,1:N3,1:n_grids_max)
  
  REAL                   ::  cIH1(1:2,0:N1)
  REAL                   ::  cIH2(1:2,0:N2)
  REAL                   ::  cIH3(1:2,0:N3)
  
  !--- Restriktion (Multigrid) ------------------------------------------------------------------------------- 
  REAL                   ::  cR1 (-1:1,1:N1,2:n_grids_max)
  REAL                   ::  cR2 (-1:1,1:N2,2:n_grids_max)
  REAL                   ::  cR3 (-1:1,1:N3,2:n_grids_max)
  
  REAL                   ::  cRest1(b1L:b1U,0:N1,1:n_grids_max-1) ! TEST!!!
  REAL                   ::  cRest2(b2L:b2U,0:N2,1:n_grids_max-1)
  REAL                   ::  cRest3(b3L:b3U,0:N3,1:n_grids_max-1)
  
  REAL                   ::  cRH1(1:2,1:N1)
  REAL                   ::  cRH2(1:2,1:N2)
  REAL                   ::  cRH3(1:2,1:N3)
  
  
  !===========================================================================================================
  !=== Gitterspezifikationen =================================================================================
  !===========================================================================================================
  !--- physiklische Koordinaten (global) ---------------------------------------------------------------------
  REAL                   ::  y1p(1:M1), y1u(0:M1)
  REAL                   ::  y2p(1:M2), y2v(0:M2)
  REAL                   ::  y3p(1:M3), y3w(0:M3)
  
  !--- physiklische Koordinaten (Block) ----------------------------------------------------------------------
  REAL                   ::  x1p(b1L:(N1+b1U)), x1u(b1L:(N1+b1U))
  REAL                   ::  x2p(b2L:(N2+b2U)), x2v(b2L:(N2+b2U))
  REAL                   ::  x3p(b3L:(N3+b3U)), x3w(b3L:(N3+b3U))
  
  !--- physiklische Koordinaten (Block, Multigrid) -----------------------------------------------------------
  REAL                   ::  x1pR(b1L:(N1+b1U),1:n_grids_max), x1uR(b1L:(N1+b1U),1:n_grids_max)
  REAL                   ::  x2pR(b2L:(N2+b2U),1:n_grids_max), x2vR(b2L:(N2+b2U),1:n_grids_max)
  REAL                   ::  x3pR(b3L:(N3+b3U),1:n_grids_max), x3wR(b3L:(N3+b3U),1:n_grids_max)
  
  !--- Gitterweiten (global) ---------------------------------------------------------------------------------
  REAL                   ::  dy1p(1:M1), dy1u(0:M1)
  REAL                   ::  dy2p(1:M2), dy2v(0:M2)
  REAL                   ::  dy3p(1:M3), dy3w(0:M3)
  
  !--- Gitterweiten (Block) ----------------------------------------------------------------------------------
  REAL                   ::  dx1p(1:N1), dx1u(0:N1)
  REAL                   ::  dx2p(1:N2), dx2v(0:N2)
  REAL                   ::  dx3p(1:N3), dx3w(0:N3)
  
  REAL                   ::  dx1DM(1:N1), dx1pM(1:N1), ddx1pM(1:N1)
  REAL                   ::  dx2DM(1:N2), dx2pM(1:N2), ddx2pM(1:N2)
  REAL                   ::  dx3DM(1:N3), dx3pM(1:N3), ddx3pM(1:N3)
  
  REAL                   ::  dx1GM(0:N1), dx1uM(0:N1), ddx1uM(0:N1)
  REAL                   ::  dx2GM(0:N2), dx2vM(0:N2), ddx2vM(0:N2)
  REAL                   ::  dx3GM(0:N3), dx3wM(0:N3), ddx3wM(0:N3)

  
  !===========================================================================================================
  !=== Arbeitsfelder =========================================================================================
  !===========================================================================================================
  !--- Geschwindigkeiten -------------------------------------------------------------------------------------
  REAL                   ::  vel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  !--- nicht-linearer Term -----------------------------------------------------------------------------------
  REAL                   ::  nl (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  !--- Recht-Hand-Seite --------------------------------------------------------------------------------------
  REAL                   ::  rhs(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  !--- Druck -------------------------------------------------------------------------------------------------
  REAL                   ::  pre(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !--- Ausfluss-RB (Geschwindigkeitsfeld) --------------------------------------------------------------------
  ! Da die RHS für die Konzentrationsfelder nicht über die Runge-Kutta-Zwischenschritte hinweg gespeichert 
  ! werden, müssen mindestens die zugehörigen Randbedingungen gespeichert werden.
  REAL                   ::  bc11(1:N2,1:N3,1:2), nlbc11(1:N2,1:N3,1:2)
  REAL                   ::  bc12(0:N1,1:N3,1:2), nlbc12(0:N1,1:N3,1:2)
  REAL                   ::  bc13(0:N1,1:N2,1:2), nlbc13(0:N1,1:N2,1:2)
  
  REAL                   ::  bc21(0:N2,1:N3,1:2), nlbc21(0:N2,1:N3,1:2)
  REAL                   ::  bc22(1:N1,1:N3,1:2), nlbc22(1:N1,1:N3,1:2)
  REAL                   ::  bc23(1:N1,0:N2,1:2), nlbc23(1:N1,0:N2,1:2)
  
  REAL                   ::  bc31(1:N2,0:N3,1:2), nlbc31(1:N2,0:N3,1:2)
  REAL                   ::  bc32(1:N1,0:N3,1:2), nlbc32(1:N1,0:N3,1:2)
  REAL                   ::  bc33(1:N1,1:N2,1:2), nlbc33(1:N1,1:N2,1:2)
  
  REAL                   ::  drift1(b2L:(N2+b2U),b3L:(N3+b3U),1:2)
  REAL                   ::  drift2(b1L:(N1+b1U),b3L:(N3+b3U),1:2)
  REAL                   ::  drift3(b1L:(N1+b1U),b2L:(N2+b2U),1:2)
  
  !--- Residuum ----------------------------------------------------------------------------------------------
  REAL                   ::  res (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !--- Druckgradient (eine Komponente) -----------------------------------------------------------------------
  REAL                   ::  gpre(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !--- Gewichte für Divergenzfreiheit ------------------------------------------------------------------------
  REAL                   ::  weight(1:N1,1:N2,1:N3)
  
  !--- Null-Raum-Vektor --------------------------------------------------------------------------------------
  REAL                   ::  psi      (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  psi_vel  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  REAL                   ::  psi_rel1 (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , ALLOCATABLE   ::  psi_rel2 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel3 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel4 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel5 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel6 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel7 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel8 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel9 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel10(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel11(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel12(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel13(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel14(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel15(:,:,:)
  
  REAL                   ::  th11(1:N2,1:N3,1:2)
  REAL                   ::  th12(0:N1,1:N3,1:2)
  REAL                   ::  th13(0:N1,1:N2,1:2)
  
  REAL                   ::  th21(0:N2,1:N3,1:2)
  REAL                   ::  th22(1:N1,1:N3,1:2)
  REAL                   ::  th23(1:N1,0:N2,1:2)
  
  REAL                   ::  th31(1:N2,0:N3,1:2)
  REAL                   ::  th32(1:N1,0:N3,1:2)
  REAL                   ::  th33(1:N1,1:N2,1:2)
  
  !--- Multigrid ---------------------------------------------------------------------------------------------
  REAL                   ::  vec1C (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  REAL   , ALLOCATABLE   ::  vec2A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec2B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec2C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec3A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec3B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec3C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec4A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec4B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec4C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec5A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec5B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec5C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec6A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec6B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec6C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec7A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec7B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec7C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec8A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec8B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec8C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec9A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec9B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec9C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec10A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec10B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec10C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec11A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec11B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec11C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec12A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec12B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec12C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec13A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec13B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec13C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec14A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec14B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec14C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec15A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec15B(:,:,:)
  
  
  !--- BiCGstab / Richardson ---------------------------------------------------------------------------------
  REAL                   ::  pp(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  Ap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  rr(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  rh(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  Ar(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  z1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  z2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !--- product_div_grad --------------------------------------------------------------------------------------
  REAL                   ::  dig(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !--- Hilfsfelder (Druckiteration) --------------------------------------------------------------------------
  ! - rhs wird auch in test_moment nochmals verwendet und kann daher in outer_iteration nicht belegt werden!
  ! - wird auch fuer interpolierte Geschwindigkeiten in rhs_NS und rhs_conc verwendet.
  REAL                   ::  work1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  work2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  work3(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  !--- Linienrelaxation --------------------------------------------------------------------------------------
  REAL                   ::  vec1(1:N1), dia1(1:N1), SOR1(1:N1), band1(1:2,1:N1) ! TEST!!! vec1, dia muessen auch in mod_Helmholtz ausgetauscht werden!
  REAL                   ::  vec2(1:N2), dia2(1:N2), SOR2(1:N2), band2(1:2,1:N2) ! TEST!!! SOR muesste man idealerweise auch in band eingliedern!
  REAL                   ::  vec3(1:N3), dia3(1:N3), SOR3(1:N3), band3(1:2,1:N3) ! dia3(:) --> band3(1,:), vec3(:) --> band3(2,:), etc. ...
  
  
  !===========================================================================================================
  !=== Indizierung (Intervallgrenzen, Verschiebungen) ========================================================
  !===========================================================================================================
  !--- Block-Index -------------------------------------------------------------------------------------------
  !INTEGER                ::  iB, jB, kB ! TEST!!!
  
  !--- Indexverschiebung (Block --> global) ------------------------------------------------------------------
  INTEGER                ::  iShift, jShift, kShift
  
  !--- Domaingrösse (Periodizität-bereinigt) -----------------------------------------------------------------
  INTEGER                ::  dim1, dim2, dim3
  
  !--- Druck / Konzentrationen (inklusive Rand) --------------------------------------------------------------
  INTEGER                ::  S1p, S2p, S3p
  INTEGER                ::  N1p, N2p, N3p
  
  !--- Geschwindigkeiten (inklusive Rand) --------------------------------------------------------------------
  INTEGER                ::  S11B, S21B, S31B
  INTEGER                ::  S12B, S22B, S32B
  INTEGER                ::  S13B, S23B, S33B
  
  INTEGER                ::  N11B, N21B, N31B
  INTEGER                ::  N12B, N22B, N32B
  INTEGER                ::  N13B, N23B, N33B
  
  !--- Geschwindigkeiten (exklusive Rand) --------------------------------------------------------------------
  INTEGER                ::  S11, S21, S31
  INTEGER                ::  S12, S22, S32
  INTEGER                ::  S13, S23, S33
  
  INTEGER                ::  N11, N21, N31
  INTEGER                ::  N12, N22, N32
  INTEGER                ::  N13, N23, N33
 
  !--- grobe Gitter (Multigrid, INklusive Rand) --------------------------------------------------------------
  INTEGER                ::  S1R, S2R, S3R
  INTEGER                ::  d1R, d2R, d3R
  
  !--- grobe Gitter (Multigrid, EXklusive Rand) --------------------------------------------------------------
  INTEGER                ::  S11R, S22R, S33R
  INTEGER                ::  d11R, d22R, d33R
  
  !--- Ueberlappungskonvention der Blöcke (Multigrid, siehe mod_setup) ---------------------------------------
  INTEGER, PARAMETER     ::  ls1 = -1
  INTEGER, PARAMETER     ::  ls2 = -1
  INTEGER, PARAMETER     ::  ls3 = -1
  
  !--- Austauschrichtung (Multigrid) -------------------------------------------------------------------------
  ! ex = -1: unten <--  oben
  ! ex =  0: unten <--> oben
  ! ex =  1: unten  --> oben
  INTEGER                ::  ex1, ex2, ex3
  
  
  !===========================================================================================================
  !=== Randbedingungen =======================================================================================
  !===========================================================================================================
  !                              _
  !    Symmetrie-RB:   BC = -2    |
  !    periodische RB: BC = -1    |- symmetrische, zentrale Stencils
  !    Nachbar-Block:  BC =  0   _|
  !    Dirichlet-RB:   BC =  1    |
  !    Neumann-RB:     BC =  2    |- schiefe, nicht-zentrale Stencils
  !    Robin-RB:       BC =  3   _|
  !
  !--- global ------------------------------------------------------------------------------------------------
  LOGICAL                ::  outlet   (1:3,1:2,1:3)
  
  INTEGER                ::  BC_1L_global, BC_1U_global
  INTEGER                ::  BC_2L_global, BC_2U_global
  INTEGER                ::  BC_3L_global, BC_3U_global
  
  !--- lokal (Block) -----------------------------------------------------------------------------------------
  INTEGER                ::  BC_1L, BC_1U
  INTEGER                ::  BC_2L, BC_2U
  INTEGER                ::  BC_3L, BC_3U

  !--- field properties --------------------------------------------------------------------------------------
  INTEGER                ::  n_gather(1:3,1:n_grids_max)
  INTEGER                ::  NN (1:3,1:n_grids_max)
  INTEGER                ::  NB (1:3,1:n_grids_max)
  INTEGER                ::  iB (1:3,1:n_grids_max)
  INTEGER                ::  SNF(1:2,1:3,1:n_grids_max)
  INTEGER                ::  SNB(1:2,1:3,1:n_grids_max)
  INTEGER                ::  BC (1:2,1:3,1:n_grids_max)
  INTEGER                ::  ngb(1:2,1:3,1:n_grids_max)
  INTEGER                ::  comm1(1:n_grids_max), comm2(1:n_grids_max)
  INTEGER                ::  rankc2(1:n_grids_max)
  LOGICAL                ::  participate_yes(1:n_grids_max)
  INTEGER                ::  recvR(    1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  recvI(    1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  dispR(    1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  dispI(    1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  offsR(1:3,1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  offsI(1:3,1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  sizsR(1:3,1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  sizsI(1:3,1:NB1*NB2*NB3,1:n_grids_max)
  
  
  !===========================================================================================================
  !=== physikalische Parameter ===============================================================================
  !===========================================================================================================
  REAL                   ::  L1, L2, L3
  REAL                   ::  L1_half, L2_half, L3_half
  REAL                   ::  L1_amp, L2_amp, L3_amp
  REAL                   ::  Re
  LOGICAL                ::  acos_yes
  
  
  !===========================================================================================================
  !=== numerische Parameter ==================================================================================
  !===========================================================================================================
  !--- allgemein ---------------------------------------------------------------------------------------------
  REAL                   ::  CFL
  REAL                   ::  time, dtime, subtime, time_start, time_end, dtime_max, dtime0, dtime_old
  INTEGER                ::  timestep, timestep_old, substep, n_timesteps
  LOGICAL                ::  mapping_yes, upwind_yes
  LOGICAL                ::  Euler_yes, Stokes_yes, twostep_yes
  LOGICAL, PARAMETER     ::  filter_BC_yes = .TRUE. ! TEST!!!
  INTEGER                ::  timeint_mode, forcing_mode
  INTEGER                ::  bulkflow_dir
  INTEGER                ::  n_lp_vel           , n_hp_vel
  REAL                   ::  chi_vel
  REAL                   ::  vel_bulk ! TEST!!!
  INTEGER                ::  timestep_mode !kschlegel
  REAL                   ::  dtime_const !kschlegel
  REAL                   ::  tol_dtime_pid_controller !kschlegel
  LOGICAL                ::  refine_2dir_yes !kschlegel
  REAL                   ::  refine_2dir_xi !kschlegel

  
  
  !--- Runge-Kutta-Koeffizienten -----------------------------------------------------------------------------
  REAL   , PARAMETER     ::  aRK(1:3) = (/8./15.,  5./12., 3./ 4./)
  REAL   , PARAMETER     ::  bRK(1:3) = (/  0.  ,-17./60.,-5./12./)
  INTEGER, PARAMETER     ::  RK_steps = 3
  
  !--- look-up table fuer Stabilitaetsgebiet der Zeitintegration (angle = pi/2,pi) ---------------------------
  REAL   , PARAMETER     ::  stabilitylimit(0:40) = (/1.732050813, 1.943689093, 2.089210537, 2.201001743,  &
                                      &               2.290031261, 2.361554127, 2.418567407, 2.462989697,  &
                                      &               2.496169963, 2.519146008, 2.532795254, 2.537935070,  &
                                      &               2.535397854, 2.526091466, 2.511046932, 2.491448818,  &
                                      &               2.468639045, 2.444084180, 2.419302172, 2.395757241,  &
                                      &               2.374745783, 2.357302135, 2.344145473, 2.335672458,  &
                                      &               2.331985072, 2.332936948, 2.338183901, 2.347230689,  &
                                      &               2.359471631, 2.374225928, 2.390769340, 2.408363261,  &
                                      &               2.426281290, 2.443832601, 2.460381269, 2.475360992,  &
                                      &               2.488285197, 2.498753090, 2.506452564, 2.511161051,  &
                                      &               2.512745327 /)
  INTEGER, PARAMETER     ::  n_stab = 41
  
  !--- Helmholtz-Vorfaktoren ---------------------------------------------------------------------------------
  REAL                   ::  thetaL, multL
  
  !--- zeitliche Steuerung -----------------------------------------------------------------------------------
  INTEGER                ::  Int_dtime, Int_lev_pre
  
  INTEGER                ::  stride_large(1:3), stride_med(1:3), stride_small(1:3)
  LOGICAL                ::  write_large, write_med, write_small
	INTEGER                ::  interval
  REAL                   ::  time_out_vect, dtime_out_vect
  REAL                   ::  time_out_scal, dtime_out_scal
  REAL                   ::  time_out_kalm, dtime_out_kalm 
	REAL   , ALLOCATABLE   ::  dtime_kalm_phases(:)
  
  LOGICAL                ::  write_out_kalm, write_out_scal, write_out_vect
  LOGICAL                ::  new_dtime, finish_yes
    
  INTEGER                ::  write_count,write_stats_count,write_kalm_count
  INTEGER                ::  restart
  CHARACTER(LEN=3)       ::  restart_char
  
  INTEGER                ::  n_conc_old


  !===========================================================================================================
  !=== weitere Steuerungsoptionen ============================================================================
  !===========================================================================================================
  INTEGER                ::  task
  LOGICAL                ::  read_nullspace_yes
  LOGICAL                ::  nullspace_yes, nullspace_coarse_yes
  LOGICAL                ::  nullspace_ortho_yes
  
  LOGICAL                ::  write_stout_yes
  LOGICAL                ::  log_iteration_yes
  LOGICAL                ::  write_restart_yes
  LOGICAL                ::  write_lambda2_yes
  LOGICAL                ::  write_test_yes
  INTEGER                ::  num_windows !define in usr_stats and for write covariance into xdmf
  INTEGER                ::  intervals   !define number of intervals in the output of a periodic flow
  INTEGER                ::  phase       !define number of interval  in the output of a periodic flow
  INTEGER, ALLOCATABLE   ::  repetition(:) !define array of #phases for repetitions counter

  !--- globale Laufindizes -----------------------------------------------------------------------------------
  INTEGER                ::  direction
  
  !--- explizite Behandlung von Ecken bei Dirichlet-Randbedingungen ------------------------------------------
  ! (Hintergrund: Der Druck ist an diesen Orten unbestimmt, so dass er dort künstlich zu Null gesetzt wird.)
  LOGICAL, PARAMETER     ::  corner_yes = .TRUE.
  
  !--- Systemzeit --------------------------------------------------------------------------------------------
  INTEGER                ::  elatime, ctime(1:8)
  INTEGER                ::  day, hour, minu, sec, msec


  !===========================================================================================================
  !=== Iterationsparameter ===================================================================================
  !===========================================================================================================
  !--- Abbruchkriterium / Absolute Genauigkeit der Geschwindigkeiten -----------------------------------------
  REAL                   ::  epsU, epsU0
  
  !--- Glaetter ----------------------------------------------------------------------------------------------
  LOGICAL                ::  Jacobi_yes
  
  !--- max. Anzahl Iterationen -------------------------------------------------------------------------------
  INTEGER                ::  n_it_outer
  INTEGER                ::  n_it_Poisson
  INTEGER                ::  n_it_Helmh_vel
  
  !--- erwartete Konvergenzrate (äussere Iteration) ----------------------------------------------------------
  REAL                   ::  precRatio0 (1:RK_steps)
  REAL                   ::  precOffset0(1:RK_steps)
  REAL, ALLOCATABLE      ::  precRatio  (:,:)
  REAL, ALLOCATABLE      ::  precOffset (:,:)
  
  !--- Null-Initialisierung (äussere Iteration) --------------------------------------------------------------
  LOGICAL                ::  init_pre(1:RK_steps), init_vel(1:RK_steps)!, init_conc(1:RK_steps)
  
  !--- Vorkonditionierung (Multigrid) ------------------------------------------------------------------------
  INTEGER                ::  precond_outer
  INTEGER                ::  precond_Poisson
  INTEGER                ::  precond_Helmh_vel
  
  !--- Anzahl Glättungen pro Gitterlevel (Multigrid) ---------------------------------------------------------
  INTEGER                ::  n_relax_down, n_relax_up, n_relax_bottom
  
  !--- implizite Richtungen bei Linienrelaxation (Multigrid) -------------------------------------------------
  INTEGER                ::  impl_dir(1:3)
  
  !--- Anzahl Glättungen pro Gitterlevel (Multigrid) ---------------------------------------------------------
  LOGICAL                ::  weighting_yes
  
  !===========================================================================================================
  !=== Iterationsstatistik ===================================================================================
  !===========================================================================================================
  REAL                   ::  dtime_average
  REAL                   ::  max_div_init(1:2)
  INTEGER                ::  number_poisson
  
  !--- Zähler ------------------------------------------------------------------------------------------------
  INTEGER                ::  countO(1:RK_steps)
  INTEGER                ::  countP(1:RK_steps,1:2)
  INTEGER                ::  countH(1:RK_steps,1:3)
  
  !--- Konvergenzrate ----------------------------------------------------------------------------------------
  REAL                   ::  ratioO(1:RK_steps)
  REAL                   ::  ratioH(1:RK_steps,1:3)
  REAL                   ::  ratioP(1:RK_steps,1:2)
  
  
  !===========================================================================================================
  !=== MPI ===================================================================================================
  !===========================================================================================================
  !--- Kommunikatoren ----------------------------------------------------------------------------------------
  INTEGER                ::  COMM_CART
  
  INTEGER                ::  COMM_SLICE1, COMM_BAR1
  INTEGER                ::  COMM_SLICE2, COMM_BAR2
  INTEGER                ::  COMM_SLICE3, COMM_BAR3
  
  !--- Dimension und Position der Blöcke innerhalb der Kommunikatoren (Gitter-Indizes) -----------------------
  ! (für MPI_GATHERv, MPI_ALLGATHERv, vgl. iShift, jShift, kShift)
  INTEGER                ::  bar1_size(1:NB1), bar1_offset(1:NB1)
  INTEGER                ::  bar2_size(1:NB2), bar2_offset(1:NB2)
  INTEGER                ::  bar3_size(1:NB3), bar3_offset(1:NB3)
  
  !--- Ränge der Prozesse ------------------------------------------------------------------------------------
  INTEGER                ::  rank
  INTEGER                ::  rank_bar1, rank_slice1
  INTEGER                ::  rank_bar2, rank_slice2
  INTEGER                ::  rank_bar3, rank_slice3
  
  !--- Ränge der Nachbarprozesse (in kartesischem Gitter) ----------------------------------------------------
  INTEGER                ::  rank1L, rank1U
  INTEGER                ::  rank2L, rank2U
  INTEGER                ::  rank3L, rank3U
  
  !--- Error-Handle ------------------------------------------------------------------------------------------
  INTEGER                ::  merror
  
  !--- Request-Handles ---------------------------------------------------------------------------------------
  ! (müssen offenbar global angelegt werden) 
  INTEGER                ::  req1L, req1U
  INTEGER                ::  req2L, req2U
  INTEGER                ::  req3L, req3U
  
  !===========================================================================================================
  !=== HDF5 ==================================================================================================
  !===========================================================================================================
  INTEGER                ::  herror
  
  !*** bbecsek 2015: IMMERSED BOUNDARY PARAMETERS ************************************************************
  REAL                   ::  mu_fluid                !< dynamic viscosity of blood
  REAL                   ::  rho_fluid               !< bulk density of blood

  REAL                   ::  L_ref                   !< reference length scale
  REAL                   ::  U_ref                   !< reference velocity

  LOGICAL                ::  write_force_yes         !< whether to output force density or not
  LOGICAL                ::  write_xdmf_yes          !< whether to output .xmf files
  LOGICAL                ::  scale_output_yes        !< whether to scale all values with U_ref, L_ref (dimensionalize)

  INTEGER, ALLOCATABLE   ::  block_sizes(:)          !< sizes (# of nodes in each block ) 0 -> NB1*NB2*NB3
  INTEGER, ALLOCATABLE   ::  block_indices(:)        !< indices (start & end) 0 -> 6*NB1*NB2*NB3
  INTEGER, ALLOCATABLE   ::  block_elem_sizes(:)     !< sizes (# of elements in each block ) 0 -> NB1*NB2*NB3

  REAL :: fd(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  REAL :: vel_old(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)        !< vel of previous timestep (used in Picard iterations)

  INTEGER :: n_data,n_data_tot,data_shift
  REAL, ALLOCATABLE :: mean_gbl(:,:,:,:,:)
  REAL, ALLOCATABLE :: covar_gbl(:,:,:,:,:)
  REAL, ALLOCATABLE :: write_mean (:,:,:,:)
  REAL, ALLOCATABLE :: write_fluct (:,:,:,:)
  REAL, ALLOCATABLE :: write_covar(:,:,:,:)
  REAL, ALLOCATABLE :: write_gain (:,:,:,:)

  TYPE stats_t
     INTEGER :: i_data,m
     INTEGER, POINTER :: x(:),y(:),z(:)  !i,j,k of the grid node for each m component
     REAL, POINTER :: mean_xyz(:),covar_xyz(:),mean_xyzt(:,:),covar_xyzt(:,:) !stats
     REAL, POINTER :: wgt(:) !stats
     TYPE(stats_t), POINTER :: next
  END TYPE stats_t

  TYPE stats_group_t
     INTEGER :: group_id,phase
     INTEGER :: n_data,n_data_tot,data_shift
     TYPE(stats_t), pointer :: stats_first
     TYPE(stats_group_t), POINTER :: next
  END TYPE stats_group_t

  TYPE(stats_group_t), pointer :: stats_group_first

  TYPE kalman_t
     INTEGER :: i_data,flg,m,phase
     INTEGER, POINTER :: x(:),y(:),z(:)  !i,j,k of the grid node for each m component
     REAL, POINTER :: mean(:,:,:),covar(:,:,:) ! phases, m, 3 or 6
     REAL, POINTER :: muf(:),pf(:,:),obs_data(:),obs_covar(:,:),obs_oper(:,:),K(:,:) !kalman
     TYPE(kalman_t), POINTER :: next
  END TYPE kalman_t

  TYPE(kalman_t), pointer :: kalman_first


#endif

END MODULE mod_vars
