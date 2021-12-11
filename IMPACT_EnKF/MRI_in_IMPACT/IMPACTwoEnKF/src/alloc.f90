!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011											     *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!*************************************************************************************************************
  

!> @file alloc.f90
!! @brief File in which all the containers (scalar fields, vector fields, coordinates ...) are allocated.

  !===========================================================================================================
  !=== Differenzen-Koeffizienten-Arrays ======================================================================
  !===========================================================================================================
  !--- 1. Ableitung (zentral) --------------------------------------------------------------------------------
  ALLOCATE(cp1  (b1L:b1U,0:N1))
  ALLOCATE(cp2  (b2L:b2U,0:N2))
  ALLOCATE(cp3  (b3L:b3U,0:N3))
  
  ALLOCATE(cu1  (b1L:b1U,0:N1))
  ALLOCATE(cv2  (b2L:b2U,0:N2))
  ALLOCATE(cw3  (b3L:b3U,0:N3))
  
  !--- 1. Ableitung (upwind) ---------------------------------------------------------------------------------
  ALLOCATE(cNp1D(n1L:n1U,0:N1))
  ALLOCATE(cNp2D(n2L:n2U,0:N2))
  ALLOCATE(cNp3D(n3L:n3U,0:N3))
  
  ALLOCATE(cNp1U(n1L:n1U,0:N1))
  ALLOCATE(cNp2U(n2L:n2U,0:N2))
  ALLOCATE(cNp3U(n3L:n3U,0:N3))
  
  ALLOCATE(cNu1D(n1L:n1U,0:N1))
  ALLOCATE(cNv2D(n2L:n2U,0:N2))
  ALLOCATE(cNw3D(n3L:n3U,0:N3))
  
  ALLOCATE(cNu1U(n1L:n1U,0:N1))
  ALLOCATE(cNv2U(n2L:n2U,0:N2))
  ALLOCATE(cNw3U(n3L:n3U,0:N3))
  
  !--- Divergenz ---------------------------------------------------------------------------------------------
  ALLOCATE(cDu1 (d1L:d1U,0:N1))
  ALLOCATE(cDv2 (d2L:d2U,0:N2))
  ALLOCATE(cDw3 (d3L:d3U,0:N3))
  
  !--- Divergenz (transponiert) ------------------------------------------------------------------------------
  ALLOCATE(cDu1T(g1L:g1U,0:N1))
  ALLOCATE(cDv2T(g2L:g2U,0:N2))
  ALLOCATE(cDw3T(g3L:g3U,0:N3))
  
  !--- Gradient ----------------------------------------------------------------------------------------------
  ALLOCATE(cGp1 (g1L:g1U,0:N1))
  ALLOCATE(cGp2 (g2L:g2U,0:N2))
  ALLOCATE(cGp3 (g3L:g3U,0:N3))
  
  !--- Gradient (transponiert) -------------------------------------------------------------------------------
  ALLOCATE(cGp1T(d1L:d1U,0:N1))
  ALLOCATE(cGp2T(d2L:d2U,0:N2))
  ALLOCATE(cGp3T(d3L:d3U,0:N3))
  
  !--- 2. Ableitung (zentral) --------------------------------------------------------------------------------
  ALLOCATE(cp11 (b1L:b1U,0:N1))
  ALLOCATE(cp22 (b2L:b2U,0:N2))
  ALLOCATE(cp33 (b3L:b3U,0:N3))
  
  ALLOCATE(cu11 (b1L:b1U,0:N1))
  ALLOCATE(cv22 (b2L:b2U,0:N2))
  ALLOCATE(cw33 (b3L:b3U,0:N3))
  
  !--- Interpolation ----------------------------------------------------------------------------------------- 
  ALLOCATE(cIpu(g1L:g1U,0:N1))
  ALLOCATE(cIpv(g2L:g2U,0:N2))
  ALLOCATE(cIpw(g3L:g3U,0:N3))
  
  ALLOCATE(cIup(d1L:d1U,0:N1))
  ALLOCATE(cIvp(d2L:d2U,0:N2))
  ALLOCATE(cIwp(d3L:d3U,0:N3))
  
  !--- Filter ------------------------------------------------------------------------------------------------
  ALLOCATE(cFp1(b1L:b1U,0:N1))
  ALLOCATE(cFp2(b2L:b2U,0:N2))
  ALLOCATE(cFp3(b3L:b3U,0:N3))
  
  ALLOCATE(cFu1(b1L:b1U,0:N1))
  ALLOCATE(cFv2(b2L:b2U,0:N2))
  ALLOCATE(cFw3(b3L:b3U,0:N3))
  
  !--- Integrator (nur für Druckgitter) ----------------------------------------------------------------------
  ALLOCATE(cInt1(b1L:b1U,0:N1))
  ALLOCATE(cInt2(b2L:b2U,0:N2))
  ALLOCATE(cInt3(b3L:b3U,0:N3))
  
  !--- 2. Ableitung (Multigrid) ------------------------------------------------------------------------------ 
  ALLOCATE(cp11R(-1:1,0:N1,1:n_grids_max))
  ALLOCATE(cp22R(-1:1,0:N2,1:n_grids_max))
  ALLOCATE(cp33R(-1:1,0:N3,1:n_grids_max))
  
  ALLOCATE(cu11R(-1:1,0:N1,1:n_grids_max))
  ALLOCATE(cv22R(-1:1,0:N2,1:n_grids_max))
  ALLOCATE(cw33R(-1:1,0:N3,1:n_grids_max))
  
  ALLOCATE(cdg1 (-1:1,1:N1,1:n_grids_max))
  ALLOCATE(cdg2 (-1:1,1:N2,1:n_grids_max))
  ALLOCATE(cdg3 (-1:1,1:N3,1:n_grids_max))
  
  !--- Interpolation (Multigrid) ----------------------------------------------------------------------------- 
  ALLOCATE(cI1(1:2,1:N1,1:n_grids_max))
  ALLOCATE(cI2(1:2,1:N2,1:n_grids_max))
  ALLOCATE(cI3(1:2,1:N3,1:n_grids_max))
  
  ALLOCATE(cIH1(1:2,0:N1))
  ALLOCATE(cIH2(1:2,0:N2))
  ALLOCATE(cIH3(1:2,0:N3))
  
  !--- Restriktion (Multigrid) ------------------------------------------------------------------------------- 
  ALLOCATE(cR1 (-1:1,1:N1,2:n_grids_max))
  ALLOCATE(cR2 (-1:1,1:N2,2:n_grids_max))
  ALLOCATE(cR3 (-1:1,1:N3,2:n_grids_max))
  
  ALLOCATE(cRest1(b1L:b1U,0:N1,1:n_grids_max-1)) ! TEST!!!
  ALLOCATE(cRest2(b2L:b2U,0:N2,1:n_grids_max-1))
  ALLOCATE(cRest3(b3L:b3U,0:N3,1:n_grids_max-1))
  
  ALLOCATE(cRH1(1:2,1:N1))
  ALLOCATE(cRH2(1:2,1:N2))
  ALLOCATE(cRH3(1:2,1:N3))
  
  
  !===========================================================================================================
  !=== Gitterspezifikationen =================================================================================
  !===========================================================================================================
  !--- physiklische Koordinaten (global) ---------------------------------------------------------------------
  ALLOCATE(y1p(1:M1))
  ALLOCATE(y2p(1:M2))
  ALLOCATE(y3p(1:M3))
  
  ALLOCATE(y1u(0:M1))
  ALLOCATE(y2v(0:M2))
  ALLOCATE(y3w(0:M3))

  !--- physiklische Koordinaten (Block) ----------------------------------------------------------------------
  ALLOCATE(x1p(b1L:(N1+b1U)))
  ALLOCATE(x2p(b2L:(N2+b2U)))
  ALLOCATE(x3p(b3L:(N3+b3U)))
  
  ALLOCATE(x1u(b1L:(N1+b1U)))
  ALLOCATE(x2v(b2L:(N2+b2U)))
  ALLOCATE(x3w(b3L:(N3+b3U)))

  !--- physiklische Koordinaten (Block, Multigrid) -----------------------------------------------------------
  ALLOCATE(x1pR(b1L:(N1+b1U),1:n_grids_max))
  ALLOCATE(x2pR(b2L:(N2+b2U),1:n_grids_max))
  ALLOCATE(x3pR(b3L:(N3+b3U),1:n_grids_max))
  
  ALLOCATE(x1uR(b1L:(N1+b1U),1:n_grids_max))
  ALLOCATE(x2vR(b2L:(N2+b2U),1:n_grids_max))
  ALLOCATE(x3wR(b3L:(N3+b3U),1:n_grids_max))
  
  !--- Gitterweiten (global) ---------------------------------------------------------------------------------
  ALLOCATE(dy1p(1:M1))
  ALLOCATE(dy2p(1:M2))
  ALLOCATE(dy3p(1:M3))
  
  ALLOCATE(dy1u(0:M1))
  ALLOCATE(dy2v(0:M2))
  ALLOCATE(dy3w(0:M3))

  !--- Gitterweiten (Block) ----------------------------------------------------------------------------------
  ALLOCATE(dx1p(1:N1))
  ALLOCATE(dx2p(1:N2))
  ALLOCATE(dx3p(1:N3))
  
  ALLOCATE(dx1u(0:N1))
  ALLOCATE(dx2v(0:N2))
  ALLOCATE(dx3w(0:N3))
  
  ALLOCATE(dx1DM(1:N1))
  ALLOCATE(dx2DM(1:N2))
  ALLOCATE(dx3DM(1:N3))
  
  ALLOCATE(dx1pM(1:N1))
  ALLOCATE(dx2pM(1:N2))
  ALLOCATE(dx3pM(1:N3))
  
  ALLOCATE(ddx1pM(1:N1))
  ALLOCATE(ddx2pM(1:N2))
  ALLOCATE(ddx3pM(1:N3))
  
  ALLOCATE(dx1GM(0:N1))
  ALLOCATE(dx2GM(0:N2))
  ALLOCATE(dx3GM(0:N3))
  
  ALLOCATE(dx1uM(0:N1))
  ALLOCATE(dx2vM(0:N2))
  ALLOCATE(dx3wM(0:N3))
  
  ALLOCATE(ddx1uM(0:N1))
  ALLOCATE(ddx2vM(0:N2))
  ALLOCATE(ddx3wM(0:N3))

  !===========================================================================================================
  !=== Arbeitsfelder =========================================================================================
  !===========================================================================================================
  !--- Geschwindigkeiten -------------------------------------------------------------------------------------
  ALLOCATE(vel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3))
  !--- nicht-linearer Term -----------------------------------------------------------------------------------
  ALLOCATE(nl (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3))
  !--- Recht-Hand-Seite --------------------------------------------------------------------------------------
  ALLOCATE(rhs(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3))
  !--- Druck -------------------------------------------------------------------------------------------------
  ALLOCATE(pre(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))

  !--- Ausfluss-RB (Geschwindigkeitsfeld) --------------------------------------------------------------------
  ALLOCATE(bc11(1:N2,1:N3,1:2), nlbc11(1:N2,1:N3,1:2))
  ALLOCATE(bc12(0:N1,1:N3,1:2), nlbc12(0:N1,1:N3,1:2))
  ALLOCATE(bc13(0:N1,1:N2,1:2), nlbc13(0:N1,1:N2,1:2))
  
  ALLOCATE(bc21(0:N2,1:N3,1:2), nlbc21(0:N2,1:N3,1:2))
  ALLOCATE(bc22(1:N1,1:N3,1:2), nlbc22(1:N1,1:N3,1:2))
  ALLOCATE(bc23(1:N1,0:N2,1:2), nlbc23(1:N1,0:N2,1:2))
  
  ALLOCATE(bc31(1:N2,0:N3,1:2), nlbc31(1:N2,0:N3,1:2))
  ALLOCATE(bc32(1:N1,0:N3,1:2), nlbc32(1:N1,0:N3,1:2))
  ALLOCATE(bc33(1:N1,1:N2,1:2), nlbc33(1:N1,1:N2,1:2))
  
  ALLOCATE(drift1(b2L:(N2+b2U),b3L:(N3+b3U),1:2))
  ALLOCATE(drift2(b1L:(N1+b1U),b3L:(N3+b3U),1:2))
  ALLOCATE(drift3(b1L:(N1+b1U),b2L:(N2+b2U),1:2))

  !--- Residuum ----------------------------------------------------------------------------------------------
  ALLOCATE(res (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  
  !--- Druckgradient (eine Komponente) -----------------------------------------------------------------------
  ALLOCATE(gpre(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  
  !--- Gewichte für Divergenzfreiheit ------------------------------------------------------------------------
  ALLOCATE(weight(1:N1,1:N2,1:N3))
  
  !--- Null-Raum-Vektor --------------------------------------------------------------------------------------
  ALLOCATE(psi     (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(psi_vel (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3))
  ALLOCATE(psi_rel1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  
  ALLOCATE(th11(1:N2,1:N3,1:2))
  ALLOCATE(th12(0:N1,1:N3,1:2))
  ALLOCATE(th13(0:N1,1:N2,1:2))
  
  ALLOCATE(th21(0:N2,1:N3,1:2))
  ALLOCATE(th22(1:N1,1:N3,1:2))
  ALLOCATE(th23(1:N1,0:N2,1:2))
  
  ALLOCATE(th31(1:N2,0:N3,1:2))
  ALLOCATE(th32(1:N1,0:N3,1:2))
  ALLOCATE(th33(1:N1,1:N2,1:2))
  
  !--- Multigrid ---------------------------------------------------------------------------------------------
  ALLOCATE(vec1C(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  
  
  !--- BiCGstab / Richardson ---------------------------------------------------------------------------------
  ALLOCATE(pp(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(Ap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(rr(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(rh(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(Ar(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(z1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(z2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  
  !--- product_div_grad --------------------------------------------------------------------------------------
  ALLOCATE(dig(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  
  !--- Hilfsfelder (Druckiteration) --------------------------------------------------------------------------
  ALLOCATE(work1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(work2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(work3(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))

  !--- Linienrelaxation --------------------------------------------------------------------------------------
  ALLOCATE(vec1(1:N1))
  ALLOCATE(vec2(1:N2))
  ALLOCATE(vec3(1:N3))
  
  ALLOCATE(dia1(1:N1))
  ALLOCATE(dia2(1:N2))
  ALLOCATE(dia3(1:N3))
  
  ALLOCATE(SOR1(1:N1))
  ALLOCATE(SOR2(1:N2))
  ALLOCATE(SOR3(1:N3))
  
  ALLOCATE(band1(1:2,1:N1))
  ALLOCATE(band2(1:2,1:N2))
  ALLOCATE(band3(1:2,1:N3))
    
  !===========================================================================================================
  !=== Randbedingungen =======================================================================================
  !===========================================================================================================
    
  !--- field properties --------------------------------------------------------------------------------------
  ALLOCATE(recvR(    1:NB1*NB2*NB3,1:n_grids_max))
  ALLOCATE(recvI(    1:NB1*NB2*NB3,1:n_grids_max))
  ALLOCATE(dispR(    1:NB1*NB2*NB3,1:n_grids_max))
  ALLOCATE(dispI(    1:NB1*NB2*NB3,1:n_grids_max))
  ALLOCATE(offsR(1:3,1:NB1*NB2*NB3,1:n_grids_max))
  ALLOCATE(offsI(1:3,1:NB1*NB2*NB3,1:n_grids_max))
  ALLOCATE(sizsR(1:3,1:NB1*NB2*NB3,1:n_grids_max))
  ALLOCATE(sizsI(1:3,1:NB1*NB2*NB3,1:n_grids_max))
  
  

  !===========================================================================================================
  !=== MPI ===================================================================================================
  !===========================================================================================================
  !--- Dimension und Position der Blöcke innerhalb der Kommunikatoren (Gitter-Indizes) -----------------------
  ALLOCATE(bar1_size(1:NB1))
  ALLOCATE(bar2_size(1:NB2))
  ALLOCATE(bar3_size(1:NB3))
  
  ALLOCATE(bar1_offset(1:NB1))
  ALLOCATE(bar2_offset(1:NB2))
  ALLOCATE(bar3_offset(1:NB3))

  !============================================================================================================
  !=== IBM (bbecsek) ==========================================================================================
  !============================================================================================================
  !--- Force density ------------------------------------------------------------------------------------------
  ALLOCATE(fd(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3))
  ALLOCATE(vel_old(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3))

  !--- Block Sizes --------------------------------------------------------------------------------------------
  ALLOCATE(block_sizes(0:NB1*NB2*NB3-1))
  ALLOCATE(block_indices(0:6*NB1*NB2*NB3-1))
  ALLOCATE(block_elem_sizes(0:NB1*NB2*NB3-1))
