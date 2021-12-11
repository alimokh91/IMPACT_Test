!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

MODULE mod_dims
  
  IMPLICIT NONE
  
  !===========================================================================================================
  !=== verschiedene Dimensionen ==============================================================================
  !===========================================================================================================
  
  !===========================================================================================================
  !=== Domain- und Blockspezifikationen ======================================================================
  !===========================================================================================================
  !--- Anzahl Blöcke -----------------------------------------------------------------------------------------
#ifdef ALLOC
  INTEGER             ::  NB1
  INTEGER             ::  NB2
  INTEGER             ::  NB3
#else
  INTEGER, PARAMETER  ::  NB1 = 2
  INTEGER, PARAMETER  ::  NB2 = 2
  INTEGER, PARAMETER  ::  NB3 = 2
#endif
  
  !--- Gesamte Domain (über alle Blöcke) ---------------------------------------------------------------------
#ifdef ALLOC
  INTEGER             ::  M1
  INTEGER             ::  M2
  INTEGER             ::  M3
#else
  INTEGER, PARAMETER  ::  M1 = 3*2**4+1
  INTEGER, PARAMETER  ::  M2 = 7*2**6+1
  INTEGER, PARAMETER  ::  M3 = 5*2**5+1
#endif
  
  
  !===========================================================================================================
  !=== Konvergenzordnung der Differenzenkoeffizienten (Anzahl Koeffizienten) =================================
  !===========================================================================================================
    
  !--- Anzahl Stencil-Koeffizienten (Rand) -------------------------------------------------------------------
  !INTEGER, PARAMETER  ::  ncb1c(2) = (/2,3/)
  !INTEGER, PARAMETER  ::  ncb1f(2) = (/2,3/)
  !INTEGER, PARAMETER  ::  ncb1r(2) = (/2,3/)
  !INTEGER, PARAMETER  ::  ncb1g(2) = (/0,2/)
  !INTEGER, PARAMETER  ::  ncb1d(2) = (/2,2/)
  
  !INTEGER, PARAMETER  ::  ncb1c(3) = (/3,4,5/)
  !INTEGER, PARAMETER  ::  ncb1f(3) = (/3,4,5/)
  !INTEGER, PARAMETER  ::  ncb1r(3) = (/2,3,3/)
  !INTEGER, PARAMETER  ::  ncb1g(3) = (/2,3,4/)
  !INTEGER, PARAMETER  ::  ncb1d(3) = (/3,4,4/)
  
  ! Stabil   (xi >= 2, Re=10000, N=17)
  !INTEGER, PARAMETER  ::  ncb1c(4) = (/4,5,5,7/)
  !INTEGER, PARAMETER  ::  ncb1f(4) = (/4,5,5,5/)
  !INTEGER, PARAMETER  ::  ncb1r(4) = (/2,3,3,3/)
  !INTEGER, PARAMETER  ::  ncb1g(4) = (/3,4,4,6/)
  !INTEGER, PARAMETER  ::  ncb1d(4) = (/4,4,6,6/)
  
  INTEGER, PARAMETER  ::  ncb1c(4) = (/4,5,5,7/)
  INTEGER, PARAMETER  ::  ncb1d(4) = (/4,4,6,6/)
  INTEGER, PARAMETER  ::  ncb1f(4) = (/4,5,5,5/) ! TEST!!!dim_ncb1c wird als Array-Laenge angenommen ...
  INTEGER, PARAMETER  ::  ncb1g(4) = (/3,4,4,6/)
  INTEGER, PARAMETER  ::  ncb1r(4) = (/2,3,3,3/) ! TEST!!!dim_ncb1c wird als Array-Laenge angenommen ...
  
  INTEGER, PARAMETER  ::  ncb2c(4) = (/4,5,5,7/)
  INTEGER, PARAMETER  ::  ncb2d(4) = (/4,4,6,6/)
  INTEGER, PARAMETER  ::  ncb2f(4) = (/4,5,5,5/)
  INTEGER, PARAMETER  ::  ncb2g(4) = (/3,4,4,6/)
  INTEGER, PARAMETER  ::  ncb2r(4) = (/2,3,3,3/)
  
  INTEGER, PARAMETER  ::  ncb3c(4) = (/4,5,5,7/)
  INTEGER, PARAMETER  ::  ncb3d(4) = (/4,4,6,6/)
  INTEGER, PARAMETER  ::  ncb3f(4) = (/4,5,5,5/)
  INTEGER, PARAMETER  ::  ncb3g(4) = (/3,4,4,6/)
  INTEGER, PARAMETER  ::  ncb3r(4) = (/2,3,3,3/)
  
  ! Instabil (äquidistant, Re=10000, N=17)
  !INTEGER, PARAMETER  ::  ncb1c(5) = (/5,6,6,7,9/)
  !INTEGER, PARAMETER  ::  ncb1f(5) = (/4,5,5,5,5/)
  !INTEGER, PARAMETER  ::  ncb1r(5) = (/2,3,3,3,3/)
  !INTEGER, PARAMETER  ::  ncb1g(5) = (/0,5,4,6,8/)
  !INTEGER, PARAMETER  ::  ncb1d(5) = (/5,4,6,8,8/)
  
  ! Instabil (äquidistant, Re=10000, N=17)
  !INTEGER, PARAMETER  ::  ncb1c(6) = (/6,7,7,7,9 ,11/)
  !INTEGER, PARAMETER  ::  ncb1f(6) = (/4,5,5,5,5 ,5 /)
  !INTEGER, PARAMETER  ::  ncb1r(6) = (/2,3,3,3,3 ,3 /)
  !INTEGER, PARAMETER  ::  ncb1g(6) = (/5,6,6,6,8 ,10/)
  !INTEGER, PARAMETER  ::  ncb1d(6) = (/6,6,6,8,10,10/)
  
  ! Stabil  (Re=10000, N=65, leicht gestreckt, explizites Forcing)
  !INTEGER, PARAMETER  ::  ncb1c(5) = (/3,7,7,7,9/)
  !INTEGER, PARAMETER  ::  ncb1f(5) = (/4,5,5,5,5/)
  !INTEGER, PARAMETER  ::  ncb1r(5) = (/2,3,3,3,3/)
  !INTEGER, PARAMETER  ::  ncb1g(5) = (/0,6,6,6,8/)
  !INTEGER, PARAMETER  ::  ncb1d(5) = (/6,6,6,8,8/)
  
  
END MODULE mod_dims
