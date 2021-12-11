!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch                    *
!* October 2014												     *
!*************************************************************************************************************

!> module containing all subroutines needed to generate the FD stencils.
MODULE mod_coeffs
  
  
  USE mod_dims
  USE mod_vars
  USE mod_exchange
  !USE mpi  
  
  PRIVATE
  
  PUBLIC FD_coeffs, test_diff, diff_coeffs, diff_coeffs_exact, test_coeffs
  PUBLIC FD_coeffs_solver, FD_coeffs_solver_integral
  PUBLIC interp_coeffs, interp_coeffs_Helm, restr_coeffs, restr_coeffs_Helm
  PUBLIC get_stencil, get_stencil_Helm, get_stencil_transp
  PUBLIC get_weights
  PUBLIC Matrix_invert
  
  
  INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  !> subroutine to compute FD coefficients
  SUBROUTINE FD_coeffs
  
  IMPLICIT NONE  
  
  !===========================================================================================================
  CALL diff_coeffs(1,1,-1,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BC_1L,BC_1U,0,N1,n1L,n1U,cNp1D)
  CALL diff_coeffs(1,1, 1,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BC_1L,BC_1U,0,N1,n1L,n1U,cNp1U)
  CALL diff_coeffs(1,1,-1,mapping_yes,dim_ncb1c,ncb1c,x1u              ,x1u              ,BC_1L,BC_1U,1,N1,n1L,n1U,cNu1D)
  CALL diff_coeffs(1,1, 1,mapping_yes,dim_ncb1c,ncb1c,x1u              ,x1u              ,BC_1L,BC_1U,1,N1,n1L,n1U,cNu1U)
  
  CALL diff_coeffs(1,0, 0,.FALSE.    ,dim_ncb1c,ncb1f,x1p              ,x1p              ,BC_1L,BC_1U,0,N1,b1L,b1U,cFp1 )
  CALL diff_coeffs(1,0, 0,.FALSE.    ,dim_ncb1c,ncb1f,x1u              ,x1u              ,BC_1L,BC_1U,1,N1,b1L,b1U,cFu1 )
  
  CALL diff_coeffs(1,1, 0,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BC_1L,BC_1U,0,N1,b1L,b1U,cp1  )
  CALL diff_coeffs(1,1, 0,mapping_yes,dim_ncb1c,ncb1c,x1u              ,x1u              ,BC_1L,BC_1U,1,N1,b1L,b1U,cu1  )
  
  CALL diff_coeffs(1,2, 0,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BC_1L,BC_1U,0,N1,b1L,b1U,cp11 )
  CALL diff_coeffs(1,2, 0,mapping_yes,dim_ncb1c,ncb1c,x1u              ,x1u              ,BC_1L,BC_1U,1,N1,b1L,b1U,cu11 )
  
  CALL diff_coeffs(1,1, 0,mapping_yes,dim_ncb1g,ncb1g,x1p(g1L:(N1+g1U)),x1u(g1L:(N1+g1U)),BC_1L,BC_1U,2,N1,g1L,g1U,cGp1 )
  CALL diff_coeffs(1,1, 0,mapping_yes,dim_ncb1d,ncb1d,x1u(d1L:(N1+d1U)),x1p(d1L:(N1+d1U)),BC_1L,BC_1U,3,N1,d1L,d1U,cDu1 )
  
  CALL diff_coeffs(1,0, 0,mapping_yes,dim_ncb1g,ncb1g,x1p(g1L:(N1+g1U)),x1u(g1L:(N1+g1U)),BC_1L,BC_1U,2,N1,g1L,g1U,cIpu )
  CALL diff_coeffs(1,0, 0,mapping_yes,dim_ncb1d,ncb1d,x1u(d1L:(N1+d1U)),x1p(d1L:(N1+d1U)),BC_1L,BC_1U,3,N1,d1L,d1U,cIup )
  !-----------------------------------------------------------------------------------------------------------
  CALL diff_coeffs(2,1,-1,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BC_2L,BC_2U,0,N2,n2L,n2U,cNp2D)
  CALL diff_coeffs(2,1, 1,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BC_2L,BC_2U,0,N2,n2L,n2U,cNp2U)
  CALL diff_coeffs(2,1,-1,mapping_yes,dim_ncb2c,ncb2c,x2v              ,x2v              ,BC_2L,BC_2U,1,N2,n2L,n2U,cNv2D)
  CALL diff_coeffs(2,1, 1,mapping_yes,dim_ncb2c,ncb2c,x2v              ,x2v              ,BC_2L,BC_2U,1,N2,n2L,n2U,cNv2U)
  
  CALL diff_coeffs(2,0, 0,.FALSE.    ,dim_ncb2c,ncb2f,x2p              ,x2p              ,BC_2L,BC_2U,0,N2,b2L,b2U,cFp2 )
  CALL diff_coeffs(2,0, 0,.FALSE.    ,dim_ncb2c,ncb2f,x2v              ,x2v              ,BC_2L,BC_2U,1,N2,b2L,b2U,cFv2 )
  
  CALL diff_coeffs(2,1, 0,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BC_2L,BC_2U,0,N2,b2L,b2U,cp2  )
  CALL diff_coeffs(2,1, 0,mapping_yes,dim_ncb2c,ncb2c,x2v              ,x2v              ,BC_2L,BC_2U,1,N2,b2L,b2U,cv2  )
  
  CALL diff_coeffs(2,2, 0,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BC_2L,BC_2U,0,N2,b2L,b2U,cp22 )
  CALL diff_coeffs(2,2, 0,mapping_yes,dim_ncb2c,ncb2c,x2v              ,x2v              ,BC_2L,BC_2U,1,N2,b2L,b2U,cv22 )
  
  CALL diff_coeffs(2,1, 0,mapping_yes,dim_ncb2g,ncb2g,x2p(g2L:(N2+g2U)),x2v(g2L:(N2+g2U)),BC_2L,BC_2U,2,N2,g2L,g2U,cGp2 )
  CALL diff_coeffs(2,1, 0,mapping_yes,dim_ncb2d,ncb2d,x2v(d2L:(N2+d2U)),x2p(d2L:(N2+d2U)),BC_2L,BC_2U,3,N2,d2L,d2U,cDv2 )
  
  CALL diff_coeffs(2,0, 0,mapping_yes,dim_ncb2g,ncb2g,x2p(g2L:(N2+g2U)),x2v(g2L:(N2+g2U)),BC_2L,BC_2U,2,N2,g2L,g2U,cIpv )
  CALL diff_coeffs(2,0, 0,mapping_yes,dim_ncb2d,ncb2d,x2v(d2L:(N2+d2U)),x2p(d2L:(N2+d2U)),BC_2L,BC_2U,3,N2,d2L,d2U,cIvp )
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
  CALL diff_coeffs(3,1,-1,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BC_3L,BC_3U,0,N3,n3L,n3U,cNp3D)
  CALL diff_coeffs(3,1, 1,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BC_3L,BC_3U,0,N3,n3L,n3U,cNp3U)
  CALL diff_coeffs(3,1,-1,mapping_yes,dim_ncb3c,ncb3c,x3w              ,x3w              ,BC_3L,BC_3U,1,N3,n3L,n3U,cNw3D)
  CALL diff_coeffs(3,1, 1,mapping_yes,dim_ncb3c,ncb3c,x3w              ,x3w              ,BC_3L,BC_3U,1,N3,n3L,n3U,cNw3U)
  
  CALL diff_coeffs(3,0, 0,.FALSE.    ,dim_ncb3c,ncb3f,x3p              ,x3p              ,BC_3L,BC_3U,0,N3,b3L,b3U,cFp3 )
  CALL diff_coeffs(3,0, 0,.FALSE.    ,dim_ncb3c,ncb3f,x3w              ,x3w              ,BC_3L,BC_3U,1,N3,b3L,b3U,cFw3 )
  
  CALL diff_coeffs(3,1, 0,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BC_3L,BC_3U,0,N3,b3L,b3U,cp3  )
  CALL diff_coeffs(3,1, 0,mapping_yes,dim_ncb3c,ncb3c,x3w              ,x3w              ,BC_3L,BC_3U,1,N3,b3L,b3U,cw3  )
  
  CALL diff_coeffs(3,2, 0,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BC_3L,BC_3U,0,N3,b3L,b3U,cp33 )
  CALL diff_coeffs(3,2, 0,mapping_yes,dim_ncb3c,ncb3c,x3w              ,x3w              ,BC_3L,BC_3U,1,N3,b3L,b3U,cw33 )
  
  CALL diff_coeffs(3,1, 0,mapping_yes,dim_ncb3g,ncb3g,x3p(g3L:(N3+g3U)),x3w(g3L:(N3+g3U)),BC_3L,BC_3U,2,N3,g3L,g3U,cGp3 )
  CALL diff_coeffs(3,1, 0,mapping_yes,dim_ncb3d,ncb3d,x3w(d3L:(N3+d3U)),x3p(d3L:(N3+d3U)),BC_3L,BC_3U,3,N3,d3L,d3U,cDw3 )
  
  CALL diff_coeffs(3,0, 0,mapping_yes,dim_ncb3g,ncb3g,x3p(g3L:(N3+g3U)),x3w(g3L:(N3+g3U)),BC_3L,BC_3U,2,N3,g3L,g3U,cIpw )
  CALL diff_coeffs(3,0, 0,mapping_yes,dim_ncb3d,ncb3d,x3w(d3L:(N3+d3U)),x3p(d3L:(N3+d3U)),BC_3L,BC_3U,3,N3,d3L,d3U,cIwp )
  END IF
  !===========================================================================================================
  
  !===========================================================================================================
                   CALL diff_coeffs(1,-1,0,mapping_yes,dim_ncb1c,ncb1c,x1p,x1p,BC_1L,BC_1U,0,N1,b1L,b1U,cInt1)
                   CALL diff_coeffs(2,-1,0,mapping_yes,dim_ncb2c,ncb2c,x2p,x2p,BC_2L,BC_2U,0,N2,b2L,b2U,cInt2)
  IF (dimens == 3) CALL diff_coeffs(3,-1,0,mapping_yes,dim_ncb3c,ncb3c,x3p,x3p,BC_3L,BC_3U,0,N3,b3L,b3U,cInt3)
  !===========================================================================================================
  
  
  END SUBROUTINE FD_coeffs
  
  
  
  
  
  
  
  
  
  SUBROUTINE test_coeffs
  
  IMPLICIT NONE
  
  
  ! ACHTUNG!!! Upwind-Differenzen noch nicht getestet!!!
  !===========================================================================================================
  IF (iB(2,1) == 1 .AND. iB(3,1) == 1) THEN
     CALL test_diff(1,0,x1p              ,x1p              ,N1,b1L,b1U,cFp1 ,BC_1L,0,'cFp1' )
     CALL test_diff(1,0,x1u              ,x1u              ,N1,b1L,b1U,cFu1 ,BC_1L,1,'cFu1' )
     CALL test_diff(1,1,x1u              ,x1u              ,N1,b1L,b1U,cu1  ,BC_1L,1,'cu1'  )
     CALL test_diff(1,1,x1p              ,x1p              ,N1,b1L,b1U,cp1  ,BC_1L,0,'cp1'  )
     CALL test_diff(1,2,x1u              ,x1u              ,N1,b1L,b1U,cu11 ,BC_1L,1,'cu11' )
     CALL test_diff(1,2,x1p              ,x1p              ,N1,b1L,b1U,cp11 ,BC_1L,0,'cp11' )
     
     !CALL test_diff(1,1,x1p(g1L:(N1+g1U)),x1u(g1L:(N1+g1U)),N1,g1L,g1U,cGp1 ,BC_1L,2,'cGp1' ) ! TEST!!! Funktioniert nicht mehr ...
     !CALL test_diff(1,1,x1u(d1L:(N1+d1U)),x1p(d1L:(N1+d1U)),N1,d1L,d1U,cDu1 ,BC_1L,3,'cDu1' )
     !CALL test_diff(1,0,x1p(g1L:(N1+g1U)),x1u(g1L:(N1+g1U)),N1,g1L,g1U,cIpu ,BC_1L,2,'cIpu' )
     !CALL test_diff(1,0,x1u(d1L:(N1+d1U)),x1p(d1L:(N1+d1U)),N1,d1L,d1U,cIup ,BC_1L,3,'cIup' )
     
     CALL test_diff(1,1,x1u(n1L:(N1+n1U)),x1u(n1L:(N1+n1U)),N1,n1L,n1U,cNu1D,BC_1L,1,'cNu1D')
     CALL test_diff(1,1,x1u(n1L:(N1+n1U)),x1u(n1L:(N1+n1U)),N1,n1L,n1U,cNu1U,BC_1L,1,'cNu1U')
     CALL test_diff(1,1,x1p(n1L:(N1+n1U)),x1p(n1L:(N1+n1U)),N1,n1L,n1U,cNp1D,BC_1L,0,'cNp1D')
     CALL test_diff(1,1,x1p(n1L:(N1+n1U)),x1p(n1L:(N1+n1U)),N1,n1L,n1U,cNp1U,BC_1L,0,'cNp1U')
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (iB(1,1) == 1 .AND. iB(3,1) == 1) THEN
     CALL test_diff(2,0,x2p              ,x2p              ,N2,b2L,b2U,cFp2 ,BC_2L,0,'cFp2' )
     CALL test_diff(2,0,x2v              ,x2v              ,N2,b2L,b2U,cFv2 ,BC_2L,1,'cFv2' )
     CALL test_diff(2,1,x2p              ,x2p              ,N2,b2L,b2U,cp2  ,BC_2L,0,'cp2'  )
     CALL test_diff(2,1,x2v              ,x2v              ,N2,b2L,b2U,cv2  ,BC_2L,1,'cv2'  )
     CALL test_diff(2,2,x2p              ,x2p              ,N2,b2L,b2U,cp22 ,BC_2L,0,'cp22' )
     CALL test_diff(2,2,x2v              ,x2v              ,N2,b2L,b2U,cv22 ,BC_2L,1,'cv22' )
     
     !CALL test_diff(2,1,x2p(g2L:(N2+g2U)),x2v(g2L:(N2+g2U)),N2,g2L,g2U,cGp2 ,BC_2L,2,'cGp2' )
     !CALL test_diff(2,1,x2v(d2L:(N2+d2U)),x2p(d2L:(N2+d2U)),N2,d2L,d2U,cDv2 ,BC_2L,3,'cDv2' )
     !CALL test_diff(2,0,x2p(g2L:(N2+g2U)),x2v(g2L:(N2+g2U)),N2,g2L,g2U,cIpv ,BC_2L,2,'cIpv' )
     !CALL test_diff(2,0,x2v(d2L:(N2+d2U)),x2p(d2L:(N2+d2U)),N2,d2L,d2U,cIvp ,BC_2L,3,'cIvp' )
     
     CALL test_diff(2,1,x2v(n2L:(N2+n2U)),x2v(n2L:(N2+n2U)),N2,n2L,n2U,cNv2D,BC_2L,1,'cNv2D')
     CALL test_diff(2,1,x2v(n2L:(N2+n2U)),x2v(n2L:(N2+n2U)),N2,n2L,n2U,cNv2U,BC_2L,1,'cNv2U')
     CALL test_diff(2,1,x2p(n2L:(N2+n2U)),x2p(n2L:(N2+n2U)),N2,n2L,n2U,cNp2D,BC_2L,0,'cNp2D')
     CALL test_diff(2,1,x2p(n2L:(N2+n2U)),x2p(n2L:(N2+n2U)),N2,n2L,n2U,cNp2U,BC_2L,0,'cNp2U')
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. dimens == 3) THEN
     CALL test_diff(3,0,x3p              ,x3p              ,N3,b3L,b3U,cFp3 ,BC_3L,0,'cFp3' )
     CALL test_diff(3,0,x3w              ,x3w              ,N3,b3L,b3U,cFw3 ,BC_3L,1,'cFw3' )
     CALL test_diff(3,1,x3p              ,x3p              ,N3,b3L,b3U,cp3  ,BC_3L,0,'cp3'  )
     CALL test_diff(3,1,x3w              ,x3w              ,N3,b3L,b3U,cw3  ,BC_3L,1,'cw3'  )
     CALL test_diff(3,2,x3p              ,x3p              ,N3,b3L,b3U,cp33 ,BC_3L,0,'cp33' )
     CALL test_diff(3,2,x3w              ,x3w              ,N3,b3L,b3U,cw33 ,BC_3L,1,'cw33' )
     
     !CALL test_diff(3,1,x3p(g3L:(N3+g3U)),x3w(g3L:(N3+g3U)),N3,g3L,g3U,cGp3 ,BC_3L,2,'cGp3' )
     !CALL test_diff(3,1,x3w(d3L:(N3+d3U)),x3p(d3L:(N3+d3U)),N3,d3L,d3U,cDw3 ,BC_3L,3,'cDw3' )
     !CALL test_diff(3,0,x3p(g3L:(N3+g3U)),x3w(g3L:(N3+g3U)),N3,g3L,g3U,cIpw ,BC_3L,2,'cIpw' )
     !CALL test_diff(3,0,x3w(d3L:(N3+d3U)),x3p(d3L:(N3+d3U)),N3,d3L,d3U,cIwp ,BC_3L,3,'cIwp' )
     
     CALL test_diff(3,1,x3w(n3L:(N3+n3U)),x3w(n3L:(N3+n3U)),N3,n3L,n3U,cNw3D,BC_3L,1,'cNw3D')
     CALL test_diff(3,1,x3w(n3L:(N3+n3U)),x3w(n3L:(N3+n3U)),N3,n3L,n3U,cNw3U,BC_3L,1,'cNw3U')
     CALL test_diff(3,1,x3p(n3L:(N3+n3U)),x3p(n3L:(N3+n3U)),N3,n3L,n3U,cNp3D,BC_3L,0,'cNp3D')
     CALL test_diff(3,1,x3p(n3L:(N3+n3U)),x3p(n3L:(N3+n3U)),N3,n3L,n3U,cNp3U,BC_3L,0,'cNp3U')
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE test_coeffs
  
  
  
  
  
  
  
  
  
  
  !> subroutine for computing the difference coefficients for derivatives of order ABL in direction DIR.
  SUBROUTINE diff_coeffs(dir,abl,upwind,mapping_yes,dim_ncb,n_coeff_bound,xC,xE,BCL,BCU,grid_type,Nmax,bL,bU,cc)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)   ::  dir !< direction of the derivative
  INTEGER, INTENT(in)   ::  abl !< order of derivative, mapping only supported with abl<=2
  INTEGER, INTENT(in)   ::  upwind !< {-1,0,1}  
  LOGICAL, INTENT(in)   ::  mapping_yes !< if true, stencils are computed on equidistant grid, then mapped to physical grid
  INTEGER, INTENT(in)   ::  dim_ncb !< size of coefficient stencil ? 
  INTEGER, INTENT(in)   ::  n_coeff_bound(1:dim_ncb) !< number of coefficents in stencil ?
  INTEGER, INTENT(in)   ::  Nmax !< local max gridpoint number
  INTEGER, INTENT(in)   ::  bL !< ?
  INTEGER, INTENT(in)   ::  bU !< ?
  REAL   , INTENT(in)   ::  xC(bL:(Nmax+bU)) !< coefficient coordinates (onto where value are transfered)
  REAL   , INTENT(in)   ::  xE(bL:(Nmax+bU)) !< coordinates of transferred value
  INTEGER, INTENT(in)   ::  BCL !< BC at lower limit
  INTEGER, INTENT(in)   ::  BCU !< BC at upper limit
  INTEGER, INTENT(in)   ::  grid_type !< type of grid 0: pressure, 1: u, 2: v, 3: w ?
  
  REAL   , INTENT(out)  ::  cc(bL:bU,0:Nmax) !< coefficients
  
  INTEGER               ::  n_coeff
  INTEGER               ::  dim_n_coeff_bound
  INTEGER               ::  i, ii, iC, iStart, SShift
  INTEGER               ::  k, kk
  
  INTEGER               ::  left, right
  
  REAL                  ::  dxi(1:2)
  
  REAL                  ::  dxL, dxU ! Für Integrationskoeffizienten
  
  REAL   , ALLOCATABLE  ::  cc_xi(:,:)
  REAL   , ALLOCATABLE  ::  deltaX(:)
  
  LOGICAL               ::  filter_yes
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Upwinding wird auf dem Rand unterdrücken, da Differenzenstencils dort ohnehin schief sind.!
  !              - Generell koennte/sollte die Initialisierung der Index-Grenzen auch schon vor dieser       !
  !                Routine ausgefuehrt werden, um hier die Uebersichtlichkeit zu verbessern.                 !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (mapping_yes .AND. abl .GT. 2) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Can`t handle derivatives > 2 combined with mapping ...'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  IF (dir == 1) SShift = iShift ! TEST!!! Besser als Argument übergeben ...
  IF (dir == 2) SShift = jShift
  IF (dir == 3) SShift = kShift
  
  !===========================================================================================================
  !=== Startindex für Entwicklungspunkte =====================================================================
  !===========================================================================================================
  IF (grid_type == 0 .OR. grid_type == 3) THEN
     iStart = 1
  ELSE
     iStart = 0
  END IF
  !===========================================================================================================
  
  filter_yes = .FALSE.
  IF (abl == 0 .AND. (grid_type == 0 .OR. grid_type == 1)) filter_yes = .TRUE.
  
  
  dim_n_coeff_bound = SIZE(n_coeff_bound)
  
  cc = 0.
  
  DO i = iStart, Nmax
     
     !========================================================================================================
     !=== Stencil-Breite auslesen ============================================================================
     !========================================================================================================
     IF      (BCL .GT. 0 .AND. i .LE. (iStart - 1 + dim_n_coeff_bound)) THEN
        n_coeff = n_coeff_bound(i + 1 - iStart)
        !IF (upwind /= 0 .AND. i == iStart) n_coeff = 0
     ELSE IF (BCU .GT. 0 .AND. i .GE. (Nmax   + 1 - dim_n_coeff_bound)) THEN
        n_coeff = n_coeff_bound(Nmax + 1 - i)
        !IF (upwind /= 0 .AND. i == Nmax  ) n_coeff = 0
     ELSE
        n_coeff = n_coeff_bound(dim_n_coeff_bound)
     END IF
     
     !========================================================================================================
     
     IF (n_coeff .GT. 0 .AND. n_coeff .GT. abl) THEN
        
        !=====================================================================================================
        !=== Anzahl der Koeffizienten RECHTS vom Entwicklungspunkt bestimmen =================================
        !=====================================================================================================
        IF (grid_type == 0 .OR. grid_type == 1) THEN
           !==================================================================================================
           !=== "normales" Gitter, keine Zwischengitterpunkte ================================================
           !==================================================================================================
           IF      (BCL .GT. 0 .AND. i .LE. (iStart - 1 + n_coeff/2)) THEN
              right = n_coeff - i + iStart - 1
           ELSE IF (BCU .GT. 0 .AND. i .GE. (Nmax   + 1 - n_coeff/2)) THEN
              right = Nmax - i
           ELSE
              ! Ausserhalb des Randes werden zentrale Differenzen gefordert (= ungerade Anzahl Koeffizienten)!
              IF (MOD(n_coeff,2) == 0) THEN
                 IF (rank == 0) THEN
                    WRITE(*,'(a   )') 'ERROR! Choose odd number of coefficients!'
                    WRITE(*,'(a,i4)') '    direction =', dir
                    WRITE(*,'(a,i4)') '    grid_type =', grid_type
                    WRITE(*,'(a,i4)') '            i =', i+SShift
                 END IF
              END IF
              right = (n_coeff-1)/2
           END IF
        ELSE IF (grid_type == 2) THEN
           !==================================================================================================
           !=== Druck ==> Impuls =============================================================================
           !==================================================================================================
           IF      (BCL .GT. 0 .AND. i .LE. n_coeff/2) THEN
              right = n_coeff - i - iStart
           ELSE IF (BCU .GT. 0 .AND. i .GE. (Nmax + 0 - n_coeff/2)) THEN
              right = Nmax - i
           ELSE
              ! Ausserhalb des Randes werden zentrale Differenzen gefordert (= gerade Anzahl Koeffizienten)!
              IF (MOD(n_coeff,2) /= 0) THEN
                 IF (rank == 0) THEN
                    WRITE(*,'(a   )') 'ERROR! Choose even number of coefficients!'
                    WRITE(*,'(a,i4)') '    direction =', dir
                    WRITE(*,'(a,i4)') '    grid_type =', grid_type
                    WRITE(*,'(a,i4)') '            i =', i+SShift
                 END IF
                 CALL MPI_FINALIZE(merror)
                 STOP
              END IF
              right = n_coeff/2
           END IF
        ELSE IF (grid_type == 3) THEN
           !==================================================================================================
           !=== Geschwindigkeit ==> Konti ====================================================================
           !==================================================================================================
           IF      (BCL .GT. 0 .AND. i .LE. n_coeff/2) THEN
              right = n_coeff - i - iStart
           ELSE IF (BCU .GT. 0 .AND. i .GE. (Nmax + 1 - n_coeff/2)) THEN
              right = Nmax - i
           ELSE
              ! Ausserhalb des Randes werden zentrale Differenzen gefordert (= gerade Anzahl Koeffizienten)!
              IF (MOD(n_coeff,2) /= 0) THEN
                 IF (rank == 0) THEN
                    WRITE(*,'(a   )') 'ERROR! Choose even number of coefficients!'
                    WRITE(*,'(a,i4)') '    direction =', dir
                    WRITE(*,'(a,i4)') '    grid_type =', grid_type
                    WRITE(*,'(a,i4)') '            i =', i+SShift
                 END IF
                 CALL MPI_FINALIZE(merror)
                 STOP
              END IF
              right = n_coeff/2 - 1
           END IF
        END IF
        !=====================================================================================================
        
        
        !=====================================================================================================
        !=== Anzahl der Koeffizienten LINKS vom Entwicklungspunkt bestimmen ==================================
        !=====================================================================================================
        ! (0-ter bzw. 1-ter Koeffizient links vom Entwicklungspunkt
        !          == "zentraler" Koeffizient im Speicher)
        left = right - n_coeff + 1
        !=====================================================================================================
        
        
        !=====================================================================================================
        IF (upwind == -1) THEN
           IF (.NOT. (BCL .GT. 0 .AND. i .LT. (iStart + n_coeff/2))) THEN
              n_coeff = n_coeff - 1
              left    = left    + 1
           END IF
        END IF
        IF (upwind ==  1) THEN
           IF (.NOT. (BCU .GT. 0 .AND. i .GT. (Nmax   - n_coeff/2))) THEN
              n_coeff = n_coeff - 1
              right   = right   - 1
           END IF
        END IF
        !=====================================================================================================
        
        
        !=====================================================================================================
        !=== Stencilanordnung testen =========================================================================
        !=====================================================================================================
        IF (right .GT. bU .OR. left .LT. bL) THEN
           IF (rank == 0) THEN
              !WRITE(*,'(a   )') 'WARNING! The FD-Stencil does probably not fit into provided array!'
              WRITE(*,'(a   )') 'ERROR! Stencil doesn`t fit into provided array!'
              WRITE(*,'(a,i4)') '    direction =', dir
              WRITE(*,'(a,i4)') '    grid_type =', grid_type
              WRITE(*,'(a,i4)') '            i =', i+SShift
           END IF
           CALL MPI_FINALIZE(merror)
           STOP
        END IF
        !=====================================================================================================
        
        
        ALLOCATE(deltaX(1:n_coeff))
        
        
        !=====================================================================================================
        !=== räumliche Abstände zum Entwicklungspunkt ========================================================
        !=====================================================================================================
        DO ii = 1, n_coeff
           ! Koeffizienten-Punkte ("iC"):
           iC = i + left + (ii-1)
           
           IF (mapping_yes) THEN
              IF      (grid_type == 2) THEN
                 deltaX(ii) = REAL(iC-i)-0.5
              ELSE IF (grid_type == 3) THEN
                 deltaX(ii) = REAL(iC-i)+0.5
              ELSE
                 deltaX(ii) = REAL(iC-i)
              END IF
           ELSE
              deltaX(ii) = xC(iC) - xE(i)
           END IF
        END DO
        !=====================================================================================================
        
        
        
        !=====================================================================================================
        !=== Integrations-Intervall ==========================================================================
        !=====================================================================================================
        ! Anmerkungen: - Kein Mapping aus Genauigkeitsgründen vorgesehen.
        !              - Nur für Druckgitter vorgesehen.
        IF (abl == -1) THEN
           
           IF (BCL .GT. 0 .AND. i == iStart) THEN
              dxL = 0.
           ELSE
              dxL = (xC(i-1) - xE(i)) / 2. ! der Einfachheit halber anstelle der Geschwindigkeitspunkte ...
           END IF
           IF (BCU .GT. 0 .AND. i == Nmax  ) THEN
              dxU = 0.
           ELSE
              dxU = (xC(i+1) - xE(i)) / 2.
           END IF
           
        END IF
        !=====================================================================================================
        
        
        !=====================================================================================================
        !=== Bestimmung der Koeffizienten ====================================================================
        !=====================================================================================================
        ! Anmerkung: Filter und Integratoren sollen nicht gemappt werden (siehe z.B. Stolz).
        IF (mapping_yes .AND. abl /= 0 .AND. abl /= -1) THEN
           
           !--- derivatives (mapped) ---
           
           ALLOCATE(cc_xi(left:right,1:abl))
           
           cc_xi = 0.
           dxi   = 0.
           
           DO k = 1, abl
              IF (n_coeff .LE. 7) THEN ! TEST!!!
                 !kk = k
                 !INCLUDE 'FD_coeffs_expl_map.f90'
                 CALL diff_coeffs_exact(k,n_coeff,deltaX(1),cc_xi(left:right,k)) ! TEST!!!
              ELSE
                 CALL FD_coeffs_solver(k,filter_yes,n_coeff,deltaX,cc_xi(left:right,k))
              END IF
              
              DO ii = left, right
                 dxi(k) = dxi(k) + cc_xi(ii,k)*xC(i+ii)
              END DO
           END DO
           
           !---------------------------------------------------!
           ! Mapping:                                          !
           ! d/dx     = (xi') * d/dxi                          !
           ! d^2/dx^2 = (xi") * d/dxi + (xi')**2 * d^2/dxi^2   !
           !                                                   !
           ! Mapping-Faktoren:                                 !
           ! xi' =   1  /  x'                                  !
           ! xi" = - x" / (x')**3                              !
           !---------------------------------------------------!
           
           IF (abl == 1) cc(left:right,i) = cc_xi(left:right,1)/dxi(1)
           IF (abl == 2) cc(left:right,i) = cc_xi(left:right,2)/dxi(1)**2 - cc_xi(left:right,1)*dxi(2)/dxi(1)**3
           
           DEALLOCATE(cc_xi)
           
        ELSE IF (mapping_yes .AND. abl == 0 .AND. .NOT. filter_yes) THEN ! TEST!!!
           
           !--- interpolation (mapped) ---
           
           !k  = 1
           !kk = 0
           !ALLOCATE(cc_xi(left:right,k:k))
           !INCLUDE 'FD_coeffs_expl_map.f90'
           !cc(left:right,i) = cc_xi(left:right,k)
           !DEALLOCATE(cc_xi)
           
           IF (n_coeff .LE. 7) THEN ! TEST!!!
              CALL diff_coeffs_exact(abl,n_coeff,deltaX(1),cc(left:right,i))
           ELSE
              CALL FD_coeffs_solver(abl,filter_yes,n_coeff,deltaX,cc(left:right,i))
           END IF
           
        ELSE
           
           IF (abl == -1) THEN
              !--- integration coefficients ---
              CALL FD_coeffs_solver_integral(n_coeff,deltaX,dxL,dxU,cc(left:right,i))
           ELSE
              !--- interpolation and derivatives ---
              CALL FD_coeffs_solver(abl,filter_yes,n_coeff,deltaX,cc(left:right,i))
           END IF
           
        END IF
        !=====================================================================================================
        
        
        DEALLOCATE(deltaX)
        
        
        !=====================================================================================================
        !=== Symmetrie =======================================================================================
        !=====================================================================================================
        IF (BCL == -2 .AND. abl == -1 .AND. i == iStart) THEN
           cc(0,i) = 0.5*cc(0,i)
           DO ii = bL, -1
              cc(ii,i) = 0.
           END DO
        END IF
        IF (BCU == -2 .AND. abl == -1 .AND. i == Nmax  ) THEN
           cc(0,i) = 0.5*cc(0,i)
           DO ii = 1, bU
              cc(ii,i) = 0.
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        ! TEST!!! abl == -1 ???
        IF (BCL == -2 .AND. (i+bL) .LT. 1) THEN
           IF (grid_type == 0 .OR. (grid_type == 2 .AND. i .GE. 1)) THEN
              DO ii = 1, 1-(i+bL)
                 cc(1-i+ii,i) = cc(1-i+ii,i) + cc(1-i-ii,i)
                 cc(1-i-ii,i) = 0.
              END DO
           END IF
           IF ((grid_type == 1 .AND. i .GE. 1) .OR. grid_type == 3) THEN
              DO ii = 1, 1-(i+bL)
                 cc(0-i+ii,i) = cc(0-i+ii,i) - cc(1-i-ii,i)
                 cc(1-i-ii,i) = 0.
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BCU == -2 .AND. (i+bU) .GT. (Nmax+0)) THEN
           IF (grid_type == 0 .OR. (grid_type == 2 .AND. i .LE. Nmax-1)) THEN
              DO ii = 1, i+bU-(Nmax-0)
                 cc(Nmax-i  -ii,i) = cc(Nmax-i  -ii,i) + cc(Nmax-i  +ii,i)
                 cc(Nmax-i  +ii,i) = 0.
              END DO
           END IF
        END IF
        IF (BCU == -2 .AND. (i+bU) .GT. (Nmax-1)) THEN
           IF ((grid_type == 1 .AND. i .LE. Nmax-1) .OR. grid_type == 3) THEN
              DO ii = 1, i+bU-(Nmax-1)
                 cc(Nmax-i  -ii,i) = cc(Nmax-i  -ii,i) - cc(Nmax-i-1+ii,i)
                 cc(Nmax-i-1+ii,i) = 0.
              END DO
           END IF
        END IF
        !=====================================================================================================
        
     END IF
     
  END DO
  
  
  END SUBROUTINE diff_coeffs
  
  
  
  
  
  
  
  
  
  
  !> subroutine to compute FD coefficients exactly, the routine is only called if there is
  !! <= 7 coefficients.
  SUBROUTINE diff_coeffs_exact(abl,n_coeff,deltaX,cc) ! TEST!!!
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  abl           !< order of derivative
  INTEGER, INTENT(in   ) ::  n_coeff       !< number of coefficients
  REAL   , INTENT(in   ) ::  deltaX        !< discrete step
  REAL   , INTENT(out  ) ::  cc(1:n_coeff) !< coefficients
  
  ! TEST!!!
  ! - korrekt ausgerichtet?
  ! - Vorzeichen ok?
  
  
  IF (n_coeff == 2) THEN
     IF (deltaX ==  0.5 .AND. abl == 0) cc(1:n_coeff) = (/ 3.,-1./)/2.
     IF (deltaX == -0.5 .AND. abl == 0) cc(1:n_coeff) = (/ 1., 1./)/2.
     IF (deltaX == -1.5 .AND. abl == 0) cc(1:n_coeff) = (/-1., 3./)/2.
     
     IF (deltaX ==  0.5 .AND. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
     IF (deltaX == -0.5 .AND. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
     IF (deltaX == -1.5 .AND. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
     
     IF (deltaX ==  0.0 .AND. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
     IF (deltaX == -1.0 .AND. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
  END IF
  IF (n_coeff == 3) THEN
     IF (deltaX ==  0.5 .AND. abl == 0) cc(1:n_coeff) = (/15.,-10., 3./)/8.
     IF (deltaX == -0.5 .AND. abl == 0) cc(1:n_coeff) = (/ 3.,  6.,-1./)/8.
     IF (deltaX == -1.5 .AND. abl == 0) cc(1:n_coeff) = (/-1.,  6., 3./)/8.
     IF (deltaX == -2.5 .AND. abl == 0) cc(1:n_coeff) = (/ 3.,-10.,15./)/8.
     
     IF (deltaX ==  0.5 .AND. abl == 1) cc(1:n_coeff) = (/-2.,  3.,-1./)/1.
     IF (deltaX == -0.5 .AND. abl == 1) cc(1:n_coeff) = (/-1.,  1., 0./)/1.
     IF (deltaX == -1.5 .AND. abl == 1) cc(1:n_coeff) = (/ 0., -1., 1./)/1.
     IF (deltaX == -2.5 .AND. abl == 1) cc(1:n_coeff) = (/ 1., -3., 2./)/1.
     
     IF (deltaX ==  0.0 .AND. abl == 1) cc(1:n_coeff) = (/-3.,  4.,-1./)/2.
     IF (deltaX == -1.0 .AND. abl == 1) cc(1:n_coeff) = (/-1.,  0., 1./)/2.
     IF (deltaX == -2.0 .AND. abl == 1) cc(1:n_coeff) = (/ 1., -4., 3./)/2.
     
     IF (deltaX ==  0.0 .AND. abl == 2) cc(1:n_coeff) = (/ 1., -2., 1./)/1.
     IF (deltaX == -1.0 .AND. abl == 2) cc(1:n_coeff) = (/ 1., -2., 1./)/1.
     IF (deltaX == -2.0 .AND. abl == 2) cc(1:n_coeff) = (/ 1., -2., 1./)/1.
  END IF
  IF (n_coeff == 4) THEN
     IF (deltaX ==  0.5 .AND. abl == 0) cc(1:n_coeff) = (/ 35.,-35.,  21., -5./)/16.
     IF (deltaX == -0.5 .AND. abl == 0) cc(1:n_coeff) = (/  5., 15.,  -5.,  1./)/16.
     IF (deltaX == -1.5 .AND. abl == 0) cc(1:n_coeff) = (/ -1.,  9.,   9., -1./)/16.
     IF (deltaX == -2.5 .AND. abl == 0) cc(1:n_coeff) = (/  1., -5.,  15.,  5./)/16.
     IF (deltaX == -3.5 .AND. abl == 0) cc(1:n_coeff) = (/ -5., 21., -35., 35./)/16.
     
     IF (deltaX ==  0.5 .AND. abl == 1) cc(1:n_coeff) = (/-71.,141., -93., 23./)/24.
     IF (deltaX == -0.5 .AND. abl == 1) cc(1:n_coeff) = (/-23., 21.,   3., -1./)/24.
     IF (deltaX == -1.5 .AND. abl == 1) cc(1:n_coeff) = (/  1.,-27.,  27., -1./)/24.
     IF (deltaX == -2.5 .AND. abl == 1) cc(1:n_coeff) = (/  1., -3., -21., 23./)/24.
     IF (deltaX == -3.5 .AND. abl == 1) cc(1:n_coeff) = (/-23., 93.,-141., 71./)/24.
     
     IF (deltaX ==  0.0 .AND. abl == 1) cc(1:n_coeff) = (/-11., 18.,  -9.,  2./)/6.
     IF (deltaX == -1.0 .AND. abl == 1) cc(1:n_coeff) = (/ -2., -3.,   6., -1./)/6.
     IF (deltaX == -2.0 .AND. abl == 1) cc(1:n_coeff) = (/  1., -6.,   3.,  2./)/6.
     IF (deltaX == -3.0 .AND. abl == 1) cc(1:n_coeff) = (/ -2.,  9., -18., 11./)/6.
     
     IF (deltaX ==  0.0 .AND. abl == 2) cc(1:n_coeff) = (/  2., -5.,   4., -1./)/1.
     IF (deltaX == -1.0 .AND. abl == 2) cc(1:n_coeff) = (/  1., -2.,   1.,  0./)/1.
     IF (deltaX == -2.0 .AND. abl == 2) cc(1:n_coeff) = (/  0.,  1.,  -2.,  1./)/1.
     IF (deltaX == -3.0 .AND. abl == 2) cc(1:n_coeff) = (/ -1.,  4.,  -5.,  2./)/1.
  END IF
  IF (n_coeff == 5) THEN
     IF (deltaX ==  0.5 .AND. abl == 0) cc(1:n_coeff) = (/315.,-420., 378.,-180., 35./)/128.
     IF (deltaX == -0.5 .AND. abl == 0) cc(1:n_coeff) = (/ 35., 140., -70.,  28., -5./)/128.
     IF (deltaX == -1.5 .AND. abl == 0) cc(1:n_coeff) = (/ -5.,  60.,  90., -20.,  3./)/128.
     IF (deltaX == -2.5 .AND. abl == 0) cc(1:n_coeff) = (/  3., -20.,  90.,  60., -5./)/128.
     IF (deltaX == -3.5 .AND. abl == 0) cc(1:n_coeff) = (/ -5.,  28., -70., 140., 35./)/128.
     IF (deltaX == -4.5 .AND. abl == 0) cc(1:n_coeff) = (/ 35.,-180., 378.,-420.,315./)/128.
     
     IF (deltaX ==  0.5 .AND. abl == 1) cc(1:n_coeff) = (/-93., 229.,-225., 111.,-22./)/24.
     IF (deltaX == -0.5 .AND. abl == 1) cc(1:n_coeff) = (/-22.,  17.,   9.,  -5.,  1./)/24.
     IF (deltaX == -1.5 .AND. abl == 1) cc(1:n_coeff) = (/  1., -27.,  27.,  -1.,  0./)/24.
     IF (deltaX == -2.5 .AND. abl == 1) cc(1:n_coeff) = (/  0.,   1., -27.,  27., -1./)/24.
     IF (deltaX == -3.5 .AND. abl == 1) cc(1:n_coeff) = (/ -1.,   5.,  -9., -17., 22./)/24.
     IF (deltaX == -4.5 .AND. abl == 1) cc(1:n_coeff) = (/ 22.,-111., 225.,-229., 93./)/24.
     
     IF (deltaX ==  0.0 .AND. abl == 1) cc(1:n_coeff) = (/-25.,  48., -36.,  16., -3./)/12.
     IF (deltaX == -1.0 .AND. abl == 1) cc(1:n_coeff) = (/ -3., -10.,  18.,  -6.,  1./)/12.
     IF (deltaX == -2.0 .AND. abl == 1) cc(1:n_coeff) = (/  1.,  -8.,   0.,   8., -1./)/12.
     IF (deltaX == -3.0 .AND. abl == 1) cc(1:n_coeff) = (/ -1.,   6., -18.,  10.,  3./)/12.
     IF (deltaX == -4.0 .AND. abl == 1) cc(1:n_coeff) = (/  3., -16.,  36., -48., 25./)/12.
     
     IF (deltaX ==  0.0 .AND. abl == 2) cc(1:n_coeff) = (/ 35.,-104., 114., -56., 11./)/12.
     IF (deltaX == -1.0 .AND. abl == 2) cc(1:n_coeff) = (/ 11., -20.,   6.,   4., -1./)/12.
     IF (deltaX == -2.0 .AND. abl == 2) cc(1:n_coeff) = (/ -1.,  16., -30.,  16., -1./)/12.
     IF (deltaX == -3.0 .AND. abl == 2) cc(1:n_coeff) = (/ -1.,   4.,   6., -20., 11./)/12.
     IF (deltaX == -4.0 .AND. abl == 2) cc(1:n_coeff) = (/ 11., -56., 114.,-104., 35./)/12.
  END IF
  IF (n_coeff == 6) THEN
     IF (deltaX ==  0.5 .AND. abl == 0) cc(1:n_coeff) = (/  693.,-1155.,  1386., -990.,   385., -63./)/256.
     IF (deltaX == -0.5 .AND. abl == 0) cc(1:n_coeff) = (/   63.,  315.,  -210.,  126.,   -45.,   7./)/256.
     IF (deltaX == -1.5 .AND. abl == 0) cc(1:n_coeff) = (/   -7.,  105.,   210.,  -70.,    21.,  -3./)/256.
     IF (deltaX == -2.5 .AND. abl == 0) cc(1:n_coeff) = (/    3.,  -25.,   150.,  150.,   -25.,   3./)/256.
     IF (deltaX == -3.5 .AND. abl == 0) cc(1:n_coeff) = (/   -3.,   21.,   -70.,  210.,   105.,  -7./)/256.
     IF (deltaX == -4.5 .AND. abl == 0) cc(1:n_coeff) = (/    7.,  -45.,   126., -210.,   315.,  63./)/256.
     IF (deltaX == -5.5 .AND. abl == 0) cc(1:n_coeff) = (/  -63.,  385.,  -990., 1386., -1155., 693./)/256.
     
     IF (deltaX ==  0.5 .AND. abl == 1) cc(1:n_coeff) = (/-9129.,26765.,-34890.,25770.,-10205.,1689./)/1920.
     IF (deltaX == -0.5 .AND. abl == 1) cc(1:n_coeff) = (/-1689., 1005.,  1430.,-1110.,   435., -71./)/1920.
     IF (deltaX == -1.5 .AND. abl == 1) cc(1:n_coeff) = (/   71.,-2115.,  2070.,   10.,   -45.,   9./)/1920.
     IF (deltaX == -2.5 .AND. abl == 1) cc(1:n_coeff) = (/   -9.,  125., -2250., 2250.,  -125.,   9./)/1920.
     IF (deltaX == -3.5 .AND. abl == 1) cc(1:n_coeff) = (/   -9.,   45.,   -10.,-2070.,  2115., -71./)/1920.
     IF (deltaX == -4.5 .AND. abl == 1) cc(1:n_coeff) = (/   71., -435.,  1110.,-1430., -1005.,1689./)/1920.
     IF (deltaX == -5.5 .AND. abl == 1) cc(1:n_coeff) = (/-1689.,10205.,-25770.,34890.,-26765.,9129./)/1920.
     
     IF (deltaX ==  0.0 .AND. abl == 1) cc(1:n_coeff) = (/ -137.,  300.,  -300.,  200.,   -75.,  12./)/60.
     IF (deltaX == -1.0 .AND. abl == 1) cc(1:n_coeff) = (/  -12.,  -65.,   120.,  -60.,    20.,  -3./)/60.
     IF (deltaX == -2.0 .AND. abl == 1) cc(1:n_coeff) = (/    3.,  -30.,   -20.,   60.,   -15.,   2./)/60.
     IF (deltaX == -3.0 .AND. abl == 1) cc(1:n_coeff) = (/   -2.,   15.,   -60.,   20.,    30.,  -3./)/60.
     IF (deltaX == -4.0 .AND. abl == 1) cc(1:n_coeff) = (/    3.,  -20.,    60., -120.,    65.,  12./)/60.
     IF (deltaX == -5.0 .AND. abl == 1) cc(1:n_coeff) = (/  -12.,   75.,  -200.,  300.,  -300., 137./)/60.
     
     IF (deltaX ==  0.0 .AND. abl == 2) cc(1:n_coeff) = (/   45., -154.,   214., -156.,    61., -10./)/12.
     IF (deltaX == -1.0 .AND. abl == 2) cc(1:n_coeff) = (/   10.,  -15.,    -4.,   14.,    -6.,   1./)/12.
     IF (deltaX == -2.0 .AND. abl == 2) cc(1:n_coeff) = (/   -1.,   16.,   -30.,   16.,    -1.,   0./)/12.
     IF (deltaX == -3.0 .AND. abl == 2) cc(1:n_coeff) = (/    0.,   -1.,    16.,  -30.,    16.,  -1./)/12.
     IF (deltaX == -4.0 .AND. abl == 2) cc(1:n_coeff) = (/    1.,   -6.,    14.,   -4.,   -15.,  10./)/12.
     IF (deltaX == -5.0 .AND. abl == 2) cc(1:n_coeff) = (/  -10.,   61.,  -156.,  214.,  -154.,  45./)/12.
  END IF
  IF (n_coeff == 7) THEN
     IF (deltaX ==  0.5 .AND. abl == 0) cc(1:n_coeff) = (/  3003., -6006.,  9009.,-8580.,   5005., -1638.,  231./)/1024.
     IF (deltaX == -0.5 .AND. abl == 0) cc(1:n_coeff) = (/   231.,  1386., -1155.,  924.,   -495.,   154.,  -21./)/1024.
     IF (deltaX == -1.5 .AND. abl == 0) cc(1:n_coeff) = (/   -21.,   378.,   945., -420.,    189.,   -54.,    7./)/1024.
     IF (deltaX == -2.5 .AND. abl == 0) cc(1:n_coeff) = (/     7.,   -70.,   525.,  700.,   -175.,    42.,   -5./)/1024.
     IF (deltaX == -3.5 .AND. abl == 0) cc(1:n_coeff) = (/    -5.,    42.,  -175.,  700.,    525.,   -70.,    7./)/1024.
     IF (deltaX == -4.5 .AND. abl == 0) cc(1:n_coeff) = (/     7.,   -54.,   189., -420.,    945.,   378.,  -21./)/1024.
     IF (deltaX == -5.5 .AND. abl == 0) cc(1:n_coeff) = (/   -21.,   154.,  -495.,  924.,  -1155.,  1386.,  231./)/1024.
     IF (deltaX == -6.5 .AND. abl == 0) cc(1:n_coeff) = (/   231., -1638.,  5005.,-8580.,   9009., -6006., 3003./)/1024.
     
     IF (deltaX ==  0.5 .AND. abl == 1) cc(1:n_coeff) = (/-10756., 36527.,-59295., 58310.,-34610., 11451.,-1627./)/1920.
     IF (deltaX == -0.5 .AND. abl == 1) cc(1:n_coeff) = (/ -1627.,   633.,  2360., -2350.,  1365.,  -443.,   62./)/1920.
     IF (deltaX == -1.5 .AND. abl == 1) cc(1:n_coeff) = (/    62., -2061.,  1935.,   190.,  -180.,    63.,   -9./)/1920.
     IF (deltaX == -2.5 .AND. abl == 1) cc(1:n_coeff) = (/    -9.,   125., -2250.,  2250.,  -125.,     9.,    0./)/1920.
     IF (deltaX == -3.5 .AND. abl == 1) cc(1:n_coeff) = (/     0.,    -9.,   125., -2250.,  2250.,  -125.,    9./)/1920.
     IF (deltaX == -4.5 .AND. abl == 1) cc(1:n_coeff) = (/     9.,   -63.,   180.,  -190., -1935.,  2061.,  -62./)/1920.
     IF (deltaX == -5.5 .AND. abl == 1) cc(1:n_coeff) = (/   -62.,   443., -1365.,  2350., -2360.,  -633., 1627./)/1920.
     IF (deltaX == -6.5 .AND. abl == 1) cc(1:n_coeff) = (/  1627.,-11451., 34610.,-58310., 59295.,-36527.,10756./)/1920.
     
     IF (deltaX == -0.0 .AND. abl == 1) cc(1:n_coeff) = (/  -147.,   360.,  -450.,   400.,  -225.,    72.,  -10./)/60.
     IF (deltaX == -1.0 .AND. abl == 1) cc(1:n_coeff) = (/   -10.,   -77.,   150.,  -100.,    50.,   -15.,    2./)/60.
     IF (deltaX == -2.0 .AND. abl == 1) cc(1:n_coeff) = (/     2.,   -24.,   -35.,    80.,   -30.,     8.,   -1./)/60.
     IF (deltaX == -3.0 .AND. abl == 1) cc(1:n_coeff) = (/    -1.,     9.,   -45.,     0.,    45.,    -9.,    1./)/60.
     IF (deltaX == -4.0 .AND. abl == 1) cc(1:n_coeff) = (/     1.,    -8.,    30.,   -80.,    35.,    24.,   -2./)/60.
     IF (deltaX == -5.0 .AND. abl == 1) cc(1:n_coeff) = (/    -2.,    15.,   -50.,   100.,  -150.,    77.,   10./)/60.
     IF (deltaX == -6.0 .AND. abl == 1) cc(1:n_coeff) = (/    10.,   -72.,   225.,  -400.,   450.,  -360.,  147./)/60.
     
     IF (deltaX == -0.0 .AND. abl == 2) cc(1:n_coeff) = (/   812., -3132.,  5265., -5080.,  2970.,  -972.,  137./)/180.
     IF (deltaX == -1.0 .AND. abl == 2) cc(1:n_coeff) = (/   137.,  -147.,  -225.,   470.,  -285.,    93.,  -13./)/180.
     IF (deltaX == -2.0 .AND. abl == 2) cc(1:n_coeff) = (/   -13.,   228.,  -420.,   200.,    15.,   -12.,    2./)/180.
     IF (deltaX == -3.0 .AND. abl == 2) cc(1:n_coeff) = (/     2.,   -27.,   270.,  -490.,   270.,   -27.,    2./)/180.
     IF (deltaX == -4.0 .AND. abl == 2) cc(1:n_coeff) = (/     2.,   -12.,    15.,   200.,  -420.,   228.,  -13./)/180.
     IF (deltaX == -5.0 .AND. abl == 2) cc(1:n_coeff) = (/   -13.,    93.,  -285.,   470.,  -225.,  -147.,  137./)/180.
     IF (deltaX == -6.0 .AND. abl == 2) cc(1:n_coeff) = (/   137.,  -972.,  2970., -5080.,  5265., -3132.,  812./)/180.
  END IF
  
  
  END SUBROUTINE diff_coeffs_exact
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE test_diff(dir,abl,xC,xE,Nmax,bL,bU,cc,BCL,grid_type,name)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)   ::  dir
  INTEGER, INTENT(in)   ::  abl
  INTEGER, INTENT(in)   ::  Nmax
  INTEGER, INTENT(in)   ::  bL
  INTEGER, INTENT(in)   ::  bU
  REAL   , INTENT(in)   ::  xC(bL:(Nmax+bU))
  REAL   , INTENT(in)   ::  xE(bL:(Nmax+bU))
  REAL   , INTENT(in)   ::  cc(bL:bU,0:Nmax)
  INTEGER, INTENT(in)   ::  BCL
  INTEGER, INTENT(in)   ::  grid_type
  
  CHARACTER(*), INTENT(in) ::  name
  
  INTEGER               ::  i, ii, iStartC, iStartE, SShift
  INTEGER               ::  nw, nw_max, n_dx
  REAL                  ::  wave, wave_mod(1:2), dx
  
  CHARACTER(LEN=1)      ::  part
  REAL                  ::  pi
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Hier sollen alle Differenzen-Stencils innerhalb eines Blockes getestet werden, analog zu  !
  !                ihrer Berechnung. Daher wird bewusst nicht mittels MPI in eine Datei geschrieben, sondern !
  !                in verschiedene Blöcke, um Fehler bei der MPI-Programmierung auszuschliessen.             !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !--- pi ---
  pi = 2.*ABS(ACOS(0.))
  
  
  IF (dir == 1) THEN
     WRITE(part,'(i1.1)') iB(1,1)
     SShift = iShift
  ELSE IF (dir == 2) THEN
     WRITE(part,'(i1.1)') iB(2,1)
     SShift = jShift
  ELSE IF (dir == 3) THEN
     WRITE(part,'(i1.1)') iB(3,1)
     SShift = kShift
  ELSE
     IF (rank == 0) WRITE(*,*) 'ERROR! Wrong input at subroutine `test_diff`!'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  
  OPEN(10,FILE='test_'//name//'_transfer_block'//part//'.txt',STATUS='UNKNOWN')
  OPEN(11,FILE='test_'//name//'_coeffs_block'  //part//'.txt',STATUS='UNKNOWN')
  
  
  !===========================================================================================================
  !=== Startindex für Entwicklungspunkte =====================================================================
  !===========================================================================================================
  IF (BCL .LE. 0 .OR. grid_type == 0 .OR. grid_type == 3) THEN
     iStartE = 1
  ELSE
     iStartE = 0
  END IF
  
  
  !===========================================================================================================
  !=== Startindex für Koeffizienten ==========================================================================
  !===========================================================================================================
  IF (BCL .GT. 0 .AND. (grid_type == 1 .OR. grid_type == 3)) THEN
     iStartC = 0
  ELSE
     iStartC = 1
  END IF
  
  
  !===========================================================================================================
  !=== Test der Koeffizienten ================================================================================
  !===========================================================================================================
  nw_max = 100
  
  DO i = iStartE, Nmax
     
     DO nw = 1, nw_max
        
        wave     = pi*REAL(nw)/REAL(nw_max)
        wave_mod = 0.
        
        !=== Referenz-Gitterweite bestimmen ==================================================================
        dx   = 0.
        n_dx = 0
        
        DO ii = bL, bU-1
           !IF (i+ii .GE. iStartC .AND. i+ii .LT. Nmax .AND. ABS(xC(i+ii+1)-xC(i+ii)) .GT. dx) THEN
           IF (i+ii .GE. iStartC .AND. i+ii .LT. Nmax) THEN
              dx   = dx + ABS(xC(i+ii+1)-xC(i+ii))
              n_dx = n_dx + 1
           END IF
        END DO
        
        dx = dx / REAL(n_dx)
        
        
        !=== Transferfunktion plotten ========================================================================
        IF (abl == 0) THEN
           DO ii = bL, bU
              wave_mod(1) = wave_mod(1) + cc(ii,i)*cos(wave/dx*(xC(i+ii)-xE(i)))
              wave_mod(2) = wave_mod(2) + cc(ii,i)*sin(wave/dx*(xC(i+ii)-xE(i)))
           END DO
        ELSE IF (abl == 1) THEN
           DO ii = bL, bU
              wave_mod(1) = wave_mod(1) + cc(ii,i)*sin(wave/dx*(xC(i+ii)-xE(i)))*dx
              wave_mod(2) = wave_mod(2) + cc(ii,i)*cos(wave/dx*(xC(i+ii)-xE(i)))*dx/wave**2
           END DO
        ELSE IF (abl == 2) THEN
           DO ii = bL, bU
              wave_mod(1) = wave_mod(1) - cc(ii,i)*cos(wave/dx*(xC(i+ii)-xE(i)))*dx**2
              wave_mod(2) = wave_mod(2) - cc(ii,i)*sin(wave/dx*(xC(i+ii)-xE(i)))*dx**2/wave
           END DO
        END IF
        
        WRITE(10,'(i7,3E26.17)') i+SShift, wave, wave_mod(1:2)
        
     END DO
     
     WRITE(10,*)
     WRITE(11,'(i7,100E26.17)') i+SShift, cc(:,i)
     
  END DO
  
  CLOSE(10)
  CLOSE(11)
  
  
  END SUBROUTINE test_diff
  
  
  
  
  
  
  
  
  
  
  !> routine that solves for the FD differentiation and transfer coefficients by matrix inversion. 
  SUBROUTINE FD_coeffs_solver(abl,filter_yes,n_coeff,deltaX,cc)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)   ::  abl               !< order of derivative
  LOGICAL, INTENT(in)   ::  filter_yes        !< ? 
  INTEGER, INTENT(in)   ::  n_coeff           !< number of coefficients in FD stencil
  REAL   , INTENT(in)   ::  deltaX(1:n_coeff) !< local dx
  REAL   , INTENT(out)  ::  cc    (1:n_coeff) !< coefficients
  
  INTEGER               ::  i, j
  
  REAL                  ::  polyn_vals    (1:n_coeff,1:n_coeff)
  REAL                  ::  polyn_vals_inv(1:n_coeff,1:n_coeff)
  
  REAL                  ::  const
  
  
  !===========================================================================================================
  !=== Aufstellen des Gleichungssystems ======================================================================
  !===========================================================================================================
  DO i = 1, n_coeff
     DO j = 1, n_coeff
        polyn_vals(i,j) = deltaX(i)**(j-1)
     END DO
  END DO
  
  
  !===========================================================================================================
  !=== Zusatzbedingungen =====================================================================================
  !===========================================================================================================
  IF (filter_yes) THEN
     ! G(pi) = 0.
     DO i = 1, n_coeff
        polyn_vals(i,n_coeff) = (-1.)**i
     END DO
  END IF
  
  
  !===========================================================================================================
  !=== Lösen des Gleichungssystems ===========================================================================
  !===========================================================================================================
  CALL Matrix_invert(n_coeff,polyn_vals,polyn_vals_inv)
  
  
  !===========================================================================================================
  !=== Funktionswert am Entwicklungspunkt (deltaX = 0) =======================================================
  !===========================================================================================================
  const = 1.
  DO i = 1, abl-1
     const = const*REAL(1+i)
  END DO
  
  
  !===========================================================================================================
  !=== Koeffizienten bestimmen (explizite Differenzen) =======================================================
  !===========================================================================================================
  cc(1:n_coeff) = const*polyn_vals_inv(1+abl,1:n_coeff)
  
  
  END SUBROUTINE FD_coeffs_solver
  
  
  
  
  
  
  
  
  
  
  ! subroutine for computing FD integration coefficients. 
  SUBROUTINE FD_coeffs_solver_integral(n_coeff,deltaX,dxL,dxU,cc)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)   ::  n_coeff           !< number of coefficient in stencil
  REAL   , INTENT(in)   ::  deltaX(1:n_coeff) !< local dx
  REAL   , INTENT(out)  ::  cc    (1:n_coeff) !< coefficients
  REAL   , INTENT(in)   ::  dxL, dxU          !< ?
  
  INTEGER               ::  i, j, k
  
  REAL                  ::  polyn_vals    (1:n_coeff,1:n_coeff)
  REAL                  ::  polyn_vals_inv(1:n_coeff,1:n_coeff)
  
  REAL                  ::  const
  
  
  !===========================================================================================================
  !=== Aufstellen des Gleichungssystems ======================================================================
  !===========================================================================================================
  DO i = 1, n_coeff
     DO j = 1, n_coeff
        polyn_vals(i,j) = deltaX(i)**(j-1)
     END DO
  END DO
  
  
  !===========================================================================================================
  !=== Lösen des Gleichungssystems ===========================================================================
  !===========================================================================================================
  CALL Matrix_invert(n_coeff,polyn_vals,polyn_vals_inv)
  
  
  !===========================================================================================================
  !=== Koeffizienten bestimmen (explizite Differenzen) =======================================================
  !===========================================================================================================
  !  (Matrix-Vektor-Multiplikation)
  cc = 0.
  DO j = 1, n_coeff
     DO i = 1, n_coeff
        cc(j) = cc(j) + (dxU**i - dxL**i)/REAL(i)*polyn_vals_inv(i,j)
     END DO
  END DO
  
  
  END SUBROUTINE FD_coeffs_solver_integral
  
  
  
  
  
  
  
  
  
  
  !> routine that creates the div(grad()) operator. For the pressure equation?
  SUBROUTINE get_stencil
  
  IMPLICIT NONE
  
  INTEGER                ::  i, imax
  INTEGER                ::  j, jmax
  INTEGER                ::  k, kmax
  
  INTEGER                ::  g
  
  REAL                   ::  cDu1R(-1:0,1:N1)
  REAL                   ::  cDv2R(-1:0,1:N2)
  REAL                   ::  cDw3R(-1:0,1:N3)
  
  REAL                   ::  cGp1R( 0:1,0:N1)
  REAL                   ::  cGp2R( 0:1,0:N2)
  REAL                   ::  cGp3R( 0:1,0:N3)
  
  
  cdg1 = 0.
  cdg2 = 0.
  cdg3 = 0.
  
  DO g = 1, n_grids
     
     !========================================================================================================
     imax = NN(1,g)
     
     cDu1R = 0.
     cGp1R = 0.
     
     !---------------------------------------------------------------------------!
     ! Achtung: Falls nicht periodisch, wird hier am Rand ein Fehler gemacht,    !
     !          der nur über die Sonderbehandlung weiter unten korrigiert wird!  !
     !---------------------------------------------------------------------------!
     
     DO i = 1, imax
        cDu1R(-1,i) = -1./(x1uR(i  ,g) - x1uR(i-1,g))
        cDu1R( 0,i) =  1./(x1uR(i  ,g) - x1uR(i-1,g))
     END DO
     
     DO i = 1, imax-1
        cGp1R( 0,i) = -1./(x1pR(i+1,g) - x1pR(i  ,g))
        cGp1R( 1,i) =  1./(x1pR(i+1,g) - x1pR(i  ,g))
     END DO
     
     IF (BC(1,1,g) .GT. 0) THEN
        cGp1R( :,0   ) =  0.
        cDu1R( :,1   ) =  0.
        !cDu1R( 0,1   ) =  1./(x1u (1  ) - x1u (0  ))
        cDu1R( 0,1   ) =  1./(y1u (1  ) - y1u (0  )) ! TEST!!! Der Schoenheit halber, s.u.
     ELSE
        cGp1R( 0,0   ) = -1./(x1pR(1,g) - x1pR(0,g))
        cGp1R( 1,0   ) =  1./(x1pR(1,g) - x1pR(0,g))
     END IF
     
     IF (BC(2,1,g) .GT. 0) THEN
        cGp1R( :,imax) =  0.
        cDu1R( :,imax) =  0.
        !cDu1R(-1,imax) = -1./(x1u (N1      ) - x1u (N1-1  )) ! TEST!!! Das geht in die Hose ...
        cDu1R(-1,imax) = -1./(y1u (M1      ) - y1u (M1-1  ))
     ELSE
        cGp1R( 0,imax) = -1./(x1pR(imax+1,g) - x1pR(imax,g))
        cGp1R( 1,imax) =  1./(x1pR(imax+1,g) - x1pR(imax,g))
     END IF
     !--------------------------------------------------------------------------------------------------------
     DO i = 1, imax
        cdg1(-1,i,g) = cDu1R(-1,i)*cGp1R(0,i-1)
        cdg1( 0,i,g) = cDu1R(-1,i)*cGp1R(1,i-1) + cDu1R(0,i)*cGp1R(0,i)
        cdg1( 1,i,g) =                            cDu1R(0,i)*cGp1R(1,i)
     END DO
     
     ! Faktor 2 kommt von der Approximation DH¯¹D ~= DI¯¹G, wobei I hier das Interpolationspolynom am Rand darstellt 
     IF (BC(1,1,g) .GT. 0) cdg1( :,1   ,g) = 2.*cdg1(:,1   ,g) ! TEST!!! Ist das wirklich so optimal?
     IF (BC(2,1,g) .GT. 0) cdg1( :,imax,g) = 2.*cdg1(:,imax,g)
     
     IF (BC(1,1,g) == -2) THEN
        cdg1( 1,1   ,g) = cdg1( 1,1   ,g) + cdg1(-1,1   ,g)
        cdg1(-1,1   ,g) = 0.
     END IF
     IF (BC(2,1,g) == -2) THEN
        cdg1(-1,imax,g) = cdg1(-1,imax,g) + cdg1( 1,imax,g)
        cdg1( 1,imax,g) = 0.
     END IF
     !========================================================================================================
     jmax = NN(2,g)
     
     cDv2R = 0.
     cGp2R = 0.
     
     DO j = 1, jmax
        cDv2R(-1,j) = -1./(x2vR(j  ,g) - x2vR(j-1,g))
        cDv2R( 0,j) =  1./(x2vR(j  ,g) - x2vR(j-1,g))
     END DO
     
     DO j = 1, jmax-1
        cGp2R( 0,j) = -1./(x2pR(j+1,g) - x2pR(j  ,g))
        cGp2R( 1,j) =  1./(x2pR(j+1,g) - x2pR(j  ,g))
     END DO
     
     IF (BC(1,2,g) .GT. 0) THEN
        cGp2R( :,0   ) =  0.
        cDv2R( :,1   ) =  0.
        !cDv2R( 0,1   ) =  1./(x2v (1  ) - x2v (0    ))
        cDv2R( 0,1   ) =  1./(y2v (1  ) - y2v (0    )) ! TEST!!! Der Schoenheit halber, s.u. ! TEST!!! Das Ergebnis ist nicht 100% korrekt. Warum???
     ELSE
        cGp2R( 0,0   ) = -1./(x2pR(1,g) - x2pR(0  ,g))
        cGp2R( 1,0   ) =  1./(x2pR(1,g) - x2pR(0  ,g))
     END IF
     
     IF (BC(2,2,g) .GT. 0) THEN
        cGp2R( :,jmax) =  0.
        cDv2R( :,jmax) =  0.
        !cDv2R(-1,jmax) = -1./(x2v (N2      ) - x2v (N2-1  )) ! TEST!!! Das geht in die Hose ...
        cDv2R(-1,jmax) = -1./(y2v (M2      ) - y2v (M2-1  ))
     ELSE
        cGp2R( 0,jmax) = -1./(x2pR(jmax+1,g) - x2pR(jmax,g))
        cGp2R( 1,jmax) =  1./(x2pR(jmax+1,g) - x2pR(jmax,g))
     END IF
     !--------------------------------------------------------------------------------------------------------
     DO j = 1, jmax
        cdg2(-1,j,g) = cDv2R(-1,j)*cGp2R(0,j-1)
        cdg2( 0,j,g) = cDv2R(-1,j)*cGp2R(1,j-1) + cDv2R(0,j)*cGp2R(0,j)
        cdg2( 1,j,g) =                            cDv2R(0,j)*cGp2R(1,j)
     END DO
     
     ! Faktor 2 kommt von der Approximation DH¯¹D ~= DI¯¹G, wobei I hier das Interpolationspolynom am Rand darstellt 
     IF (BC(1,2,g) .GT. 0) cdg2( :,1   ,g) = 2.*cdg2(:,1   ,g)
     IF (BC(2,2,g) .GT. 0) cdg2( :,jmax,g) = 2.*cdg2(:,jmax,g)
     
     IF (BC(1,2,g) == -2) THEN
        cdg2( 1,1   ,g) = cdg2( 1,1   ,g) + cdg2(-1,1   ,g)
        cdg2(-1,1   ,g) = 0.
     END IF
     IF (BC(2,2,g) == -2) THEN
        cdg2(-1,jmax,g) = cdg2(-1,jmax,g) + cdg2( 1,jmax,g)
        cdg2( 1,jmax,g) = 0.
     END IF
     !========================================================================================================
     IF (dimens == 3) THEN
     
     kmax = NN(3,g)
     
     cDw3R = 0.
     cGp3R = 0.
     
     DO k = 1, kmax
        cDw3R(-1,k) = -1./(x3wR(k  ,g) - x3wR(k-1,g))
        cDw3R( 0,k) =  1./(x3wR(k  ,g) - x3wR(k-1,g))
     END DO
     
     DO k = 1, kmax-1
        cGp3R( 0,k) = -1./(x3pR(k+1,g) - x3pR(k  ,g))
        cGp3R( 1,k) =  1./(x3pR(k+1,g) - x3pR(k  ,g))
     END DO
     
     IF (BC(1,3,g) .GT. 0) THEN
        cGp3R( :,0   ) =  0.
        cDw3R( :,1   ) =  0.
        !cDw3R( 0,1   ) =  1./(x3w (1  ) - x3w (0  ))
        cDw3R( 0,1   ) =  1./(y3w (1  ) - y3w (0  )) ! TEST!!! Der Schoenheit halber, s.u.
     ELSE
        cGp3R( 0,0   ) = -1./(x3pR(1,g) - x3pR(0,g))
        cGp3R( 1,0   ) =  1./(x3pR(1,g) - x3pR(0,g))
     END IF
     
     IF (BC(2,3,g) .GT. 0) THEN
        cGp3R( :,kmax) =  0.
        cDw3R( :,kmax) =  0.
        !cDw3R(-1,kmax) = -1./(x3w (N3      ) - x3w (N3-1  )) ! TEST!!! Das geht in die Hose ...
        cDw3R(-1,kmax) = -1./(y3w (M3      ) - y3w (M3-1  ))
     ELSE
        cGp3R( 0,kmax) = -1./(x3pR(kmax+1,g) - x3pR(kmax,g))
        cGp3R( 1,kmax) =  1./(x3pR(kmax+1,g) - x3pR(kmax,g))
     END IF
     !--------------------------------------------------------------------------------------------------------
     DO k = 1, kmax
        cdg3(-1,k,g) = cDw3R(-1,k)*cGp3R(0,k-1)
        cdg3( 0,k,g) = cDw3R(-1,k)*cGp3R(1,k-1) + cDw3R(0,k)*cGp3R(0,k)
        cdg3( 1,k,g) =                            cDw3R(0,k)*cGp3R(1,k)
     END DO
     
     ! Faktor 2 kommt von der Approximation DH¯¹D ~= DI¯¹G, wobei I hier das Interpolationspolynom am Rand darstellt 
     IF (BC(1,3,g) .GT. 0) cdg3( :,1   ,g) = 2.*cdg3(:,1   ,g)
     IF (BC(2,3,g) .GT. 0) cdg3( :,kmax,g) = 2.*cdg3(:,kmax,g)
     
     IF (BC(1,3,g) == -2) THEN
        cdg3( 1,1   ,g) = cdg3( 1,1   ,g) + cdg3(-1,1   ,g)
        cdg3(-1,1   ,g) = 0.
     END IF
     IF (BC(2,3,g) == -2) THEN
        cdg3(-1,kmax,g) = cdg3(-1,kmax,g) + cdg3( 1,kmax,g)
        cdg3( 1,kmax,g) = 0.
     END IF
     
     END IF
     !========================================================================================================
     
  END DO
  
  
  END SUBROUTINE get_stencil
  
  
  
  
  
  
  
  
  
  
  !> subroutine to compute the transpose of the div(grad()) operator (pressure equation). 
  SUBROUTINE get_stencil_transp
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k, g
  REAL   , ALLOCATABLE   ::  work(:,:)
  
  
  DO g = 1, n_grids
     
     !========================================================================================================
     ALLOCATE(work(-1:1,0:(NN(1,g)+1)))
     
     work = 0.
     DO i = 1, NN(1,g)
        work(-1:1,i) = cdg1(-1:1,i,g)
     END DO
     
     cdg1(:,:,g) = 0.
     DO i = 1, NN(1,g)
        cdg1(-1,i,g) = work( 1,i-1)
        cdg1( 0,i,g) = work( 0,i  )
        cdg1( 1,i,g) = work(-1,i+1)
     END DO
     
     DEALLOCATE(work)
     
     IF (BC(1,1,g) == 0 .OR. BC(1,1,g) == -1) THEN
        i = 1
        cdg1(-1,i,g) = 1./(x1uR(i-1,g) - x1uR(i-2,g))/(x1pR(i  ,g) - x1pR(i-1,g))
        cdg1( 1,i,g) = 1./(x1uR(i+1,g) - x1uR(i  ,g))/(x1pR(i+1,g) - x1pR(i  ,g))
     END IF
     IF (BC(2,1,g) == 0 .OR. BC(2,1,g) == -1) THEN
        i = NN(1,g)
        cdg1(-1,i,g) = 1./(x1uR(i-1,g) - x1uR(i-2,g))/(x1pR(i  ,g) - x1pR(i-1,g))
        cdg1( 1,i,g) = 1./(x1uR(i+1,g) - x1uR(i  ,g))/(x1pR(i+1,g) - x1pR(i  ,g))
     END IF
     !========================================================================================================
     ALLOCATE(work(-1:1,0:(NN(2,g)+1)))
     
     work = 0.
     DO j = 1, NN(2,g)
        work(-1:1,j) = cdg2(-1:1,j,g)
     END DO
     
     cdg2(:,:,g) = 0.
     DO j = 1, NN(2,g)
        cdg2(-1,j,g) = work( 1,j-1)
        cdg2( 0,j,g) = work( 0,j  )
        cdg2( 1,j,g) = work(-1,j+1)
     END DO
     
     DEALLOCATE(work)
     
     IF (BC(1,2,g) == 0 .OR. BC(1,2,g) == -1) THEN
        j = 1
        cdg2(-1,j,g) = 1./(x2vR(j-1,g) - x2vR(j-2,g))/(x2pR(j  ,g) - x2pR(j-1,g))
        cdg2( 1,j,g) = 1./(x2vR(j+1,g) - x2vR(j  ,g))/(x2pR(j+1,g) - x2pR(j  ,g))
     END IF
     IF (BC(2,2,g) == 0 .OR. BC(2,2,g) == -1) THEN
        j = NN(2,g)
        cdg2(-1,j,g) = 1./(x2vR(j-1,g) - x2vR(j-2,g))/(x2pR(j  ,g) - x2pR(j-1,g))
        cdg2( 1,j,g) = 1./(x2vR(j+1,g) - x2vR(j  ,g))/(x2pR(j+1,g) - x2pR(j  ,g))
     END IF
     !========================================================================================================
     IF (dimens == 3) THEN
     
     ALLOCATE(work(-1:1,0:(NN(3,g)+1)))
     
     work = 0.
     DO k = 1, NN(3,g)
        work(-1:1,k) = cdg3(-1:1,k,g)
     END DO
     
     cdg3(:,:,g) = 0.
     DO k = 1, NN(3,g)
        cdg3(-1,k,g) = work( 1,k-1)
        cdg3( 0,k,g) = work( 0,k  )
        cdg3( 1,k,g) = work(-1,k+1)
     END DO
     
     DEALLOCATE(work)
     
     IF (BC(1,3,g) == 0 .OR. BC(1,3,g) == -1) THEN
        k = 1
        cdg3(-1,k,g) = 1./(x3wR(k-1,g) - x3wR(k-2,g))/(x3pR(k  ,g) - x3pR(k-1,g))
        cdg3( 1,k,g) = 1./(x3wR(k+1,g) - x3wR(k  ,g))/(x3pR(k+1,g) - x3pR(k  ,g))
     END IF
     IF (BC(2,3,g) == 0 .OR. BC(2,3,g) == -1) THEN
        k = NN(3,g)
        cdg3(-1,k,g) = 1./(x3wR(k-1,g) - x3wR(k-2,g))/(x3pR(k  ,g) - x3pR(k-1,g))
        cdg3( 1,k,g) = 1./(x3wR(k+1,g) - x3wR(k  ,g))/(x3pR(k+1,g) - x3pR(k  ,g))
     END IF
     
     END IF
     !========================================================================================================
     
  END DO
  
  
  ! sicherheitshalber (cdg3 = 0. sollte eigentlich schon erfüllt sein):
  IF (dimens == 2) cdg3 = 0.
  
  
  END SUBROUTINE get_stencil_transp
  
  
  
  
  
  
  
  
  
  
  !> subroutine that creates the Helmholtz operator for the velocity equaiton?
  SUBROUTINE get_stencil_Helm
  
  IMPLICIT NONE
  
  INTEGER                ::  imax, jmax, kmax
  INTEGER                ::  g
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Bei Symmetrie ergeben sich automatisch auf der Symmetrie-Ebene Dirichlet-RB für die       !
  !                Normalkomponente der Geschwindigkeit (betrifft cu11R, cv22R, cw33R und tritt aufgrund des !
  !                Gittertyps nur auf den gröberen Gittern auf). Die zugehörigen Restriktionskoeffizienten   !
  !                liefern automatisch entsprechend die "Randbedingung" vel_n = 0.                           !
  !              - Bei Symmetrie-RB werden in "CALL diff_coeffs" die Koeffizienten am Rand gemäss Gittertyp  !
  !                manipuliert, d.h. bei cp11R, cp22R, cp33R werden Geschwindigkeitskomponenten parallel     !
  !                zur Symmetrieebene (bzw. Druck) angenommen, so dass cu11R, cv22R, cw33R nach dem Kopieren !
  !                von cp11R, cp22R, cp33R entsprechend für Geschwindigkeitskomponenten normal zur           !
  !                Symmetrie-Ebene geändert werden müssen.                                                   !
  !----------------------------------------------------------------------------------------------------------!
  
  ! sicherheitshalber (dimens == 2, ...):
  cp11R = 0.
  cp22R = 0.
  cp33R = 0.
  
  cu11R = 0.
  cv22R = 0.
  cw33R = 0.
  
  
  !===========================================================================================================
  !=== grobe & feinstes Gitter ===============================================================================
  !===========================================================================================================
  DO g = 1, n_grids
     
     !========================================================================================================
     imax = NN(1,g)
     
     CALL diff_coeffs(1,2,0,mapping_yes,2,(/2,3/),x1pR(-1:imax+1,g),x1pR(-1:imax+1,g),BC_1L,BC_1U,0,imax,-1,1,cp11R(-1,0,g))
     
     IF (BC_1L == 1 .OR. BC_1L == 3) THEN
        cp11R( :,1   ,g) =  0.
        cp11R( 0,1   ,g) =  1.
     END IF
     IF (BC_1U == 1 .OR. BC_1U == 3) THEN
        cp11R( :,imax,g) =  0.
        cp11R( 0,imax,g) =  1.
     END IF
     
     IF (BC_1L == 2 .OR. BC_1L == 4) THEN
        cp11R(-1,1   ,g) =  0.
        cp11R( 0,1   ,g) = -1./(x1pR(2   ,g) - x1pR(1     ,g))
        cp11R( 1,1   ,g) =  1./(x1pR(2   ,g) - x1pR(1     ,g))
     END IF
     IF (BC_1U == 2 .OR. BC_1U == 4) THEN
        cp11R(-1,imax,g) = -1./(x1pR(imax,g) - x1pR(imax-1,g))
        cp11R( 0,imax,g) =  1./(x1pR(imax,g) - x1pR(imax-1,g))
        cp11R( 1,imax,g) =  0.
     END IF
     !--------------------------------------------------------------------------------------------------------
     cu11R(:,:,g) = cp11R(:,:,g)
     
     IF (BC_1L ==  2) THEN
        cu11R( :,1   ,g) =  0.
        cu11R( 0,1   ,g) =  1.
     END IF
     IF (BC_1U ==  2) THEN
        cu11R( :,imax,g) =  0.
        cu11R( 0,imax,g) =  1.
     END IF
     
     IF (BC_1L ==  3) THEN
        cu11R(-1,1   ,g) =  0.
        cu11R( 0,1   ,g) = -1./(x1pR(2   ,g) - x1pR(1     ,g))
        cu11R( 1,1   ,g) =  1./(x1pR(2   ,g) - x1pR(1     ,g))
     END IF
     IF (BC_1U ==  3) THEN
        cu11R(-1,imax,g) = -1./(x1pR(imax,g) - x1pR(imax-1,g))
        cu11R( 0,imax,g) =  1./(x1pR(imax,g) - x1pR(imax-1,g))
        cu11R( 1,imax,g) =  0.
     END IF
     
     IF (BC_1L == -2) THEN
        cu11R(-1,1   ,g) =  0.
        cu11R( 1,1   ,g) =  0.
     END IF
     IF (BC_1U == -2) THEN
        cu11R(-1,imax,g) =  0.
        cu11R( 1,imax,g) =  0.
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     jmax = NN(2,g)
     
     CALL diff_coeffs(2,2,0,mapping_yes,2,(/2,3/),x2pR(-1:jmax+1,g),x2pR(-1:jmax+1,g),BC_2L,BC_2U,0,jmax,-1,1,cp22R(-1,0,g))
     
     IF (BC_2L == 1 .OR. BC_2L == 3) THEN
        cp22R(:,1    ,g) =  0.
        cp22R(0,1    ,g) =  1.
     END IF
     IF (BC_2U == 1 .OR. BC_2U == 3) THEN
        cp22R( :,jmax,g) =  0.
        cp22R( 0,jmax,g) =  1.
     END IF
     
     IF (BC_2L == 2 .OR. BC_2L == 4) THEN
        cp22R(-1,1   ,g) =  0.
        cp22R( 0,1   ,g) = -1./(x2pR(2   ,g) - x2pR(1     ,g))
        cp22R( 1,1   ,g) =  1./(x2pR(2   ,g) - x2pR(1     ,g))
     END IF
     IF (BC_2U == 2 .OR. BC_2U == 4) THEN
        cp22R(-1,jmax,g) = -1./(x2pR(jmax,g) - x2pR(jmax-1,g))
        cp22R( 0,jmax,g) =  1./(x2pR(jmax,g) - x2pR(jmax-1,g))
        cp22R( 1,jmax,g) =  0.
     END IF
     !--------------------------------------------------------------------------------------------------------
     cv22R(:,:,g) = cp22R(:,:,g)
     
     IF (BC_2L ==  2) THEN
        cv22R( :,1   ,g) =  0.
        cv22R( 0,1   ,g) =  1.
     END IF
     IF (BC_2U ==  2) THEN
        cv22R( :,jmax,g) =  0.
        cv22R( 0,jmax,g) =  1.
     END IF
     
     IF (BC_2L ==  3) THEN
        cv22R(-1,1   ,g) =  0.
        cv22R( 0,1   ,g) = -1./(x2pR(2   ,g) - x2pR(1     ,g))
        cv22R( 1,1   ,g) =  1./(x2pR(2   ,g) - x2pR(1     ,g))
     END IF
     IF (BC_2U ==  3) THEN
        cv22R(-1,jmax,g) = -1./(x2pR(jmax,g) - x2pR(jmax-1,g))
        cv22R( 0,jmax,g) =  1./(x2pR(jmax,g) - x2pR(jmax-1,g))
        cv22R( 1,jmax,g) =  0.
     END IF
     
     IF (BC_2L == -2) THEN
        cv22R(-1,1   ,g) =  0.
        cv22R( 1,1   ,g) =  0.
     END IF
     IF (BC_2U == -2) THEN
        cv22R(-1,jmax,g) =  0.
        cv22R( 1,jmax,g) =  0.
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
     
     kmax = NN(3,g)
     
     CALL diff_coeffs(3,2,0,mapping_yes,2,(/2,3/),x3pR(-1:kmax+1,g),x3pR(-1:kmax+1,g),BC_3L,BC_3U,0,kmax,-1,1,cp33R(-1,0,g))
     
     IF (BC_3L == 1 .OR. BC_3L == 3) THEN
        cp33R( :,1   ,g) =  0.
        cp33R( 0,1   ,g) =  1.
     END IF
     IF (BC_3U == 1 .OR. BC_3U == 3) THEN
        cp33R( :,kmax,g) =  0.
        cp33R( 0,kmax,g) =  1.
     END IF
     
     IF (BC_3L == 2 .OR. BC_3L == 4) THEN
        cp33R(-1,1   ,g) =  0.
        cp33R( 0,1   ,g) = -1./(x3pR(2   ,g) - x3pR(1     ,g))
        cp33R( 1,1   ,g) =  1./(x3pR(2   ,g) - x3pR(1     ,g))
     END IF
     IF (BC_3U == 2 .OR. BC_3U == 4) THEN
        cp33R(-1,kmax,g) = -1./(x3pR(kmax,g) - x3pR(kmax-1,g))
        cp33R( 0,kmax,g) =  1./(x3pR(kmax,g) - x3pR(kmax-1,g))
        cp33R( 1,kmax,g) =  0.
     END IF
     !--------------------------------------------------------------------------------------------------------
     cw33R(:,:,g) = cp33R(:,:,g)
     
     IF (BC_3L ==  2) THEN
        cw33R( :,1   ,g) =  0.
        cw33R( 0,1   ,g) =  1.
     END IF
     IF (BC_3U ==  2) THEN
        cw33R( :,kmax,g) =  0.
        cw33R( 0,kmax,g) =  1.
     END IF
     
     IF (BC_3L ==  3) THEN
        cw33R(-1,1   ,g) =  0.
        cw33R( 0,1   ,g) = -1./(x3pR(2   ,g) - x3pR(1     ,g))
        cw33R( 1,1   ,g) =  1./(x3pR(2   ,g) - x3pR(1     ,g))
     END IF
     IF (BC_3U ==  3) THEN
        cw33R(-1,kmax,g) = -1./(x3pR(kmax,g) - x3pR(kmax-1,g))
        cw33R( 0,kmax,g) =  1./(x3pR(kmax,g) - x3pR(kmax-1,g))
        cw33R( 1,kmax,g) =  0.
     END IF
     
     IF (BC_3L == -2) THEN
        cw33R(-1,1   ,g) =  0.
        cw33R( 1,1   ,g) =  0.
     END IF
     IF (BC_3U == -2) THEN
        cw33R(-1,kmax,g) =  0.
        cw33R( 1,kmax,g) =  0.
     END IF
     !--------------------------------------------------------------------------------------------------------
     END IF
     !========================================================================================================
     
  END DO
  
  g = 1
  
  !===========================================================================================================
  !=== feinstes Gitter (Tangential-Koeffizienten überschreiben) ==============================================
  !===========================================================================================================
  CALL diff_coeffs(1,2,0,mapping_yes,2,(/2,3/),x1u(-1:N1+1),x1u(-1:N1+1),BC_1L,BC_1U,1,N1,-1,1,cu11R(-1,0,g))
  
  IF (BC_1L == 1 .OR. BC_1L == 2) THEN
     cu11R(-1,0 ,g) = 0.
     cu11R( 0,0 ,g) = 1.- (x1p(1 )-x1u(0   )) / (x1u(1 )-x1u(0   ))
     cu11R( 1,0 ,g) =     (x1p(1 )-x1u(0   )) / (x1u(1 )-x1u(0   ))
  END IF
  IF (BC_1U == 1 .OR. BC_1U == 2) THEN
     cu11R(-1,N1,g) = 1.- (x1p(N1)-x1u(N1-1)) / (x1u(N1)-x1u(N1-1))
     cu11R( 0,N1,g) =     (x1p(N1)-x1u(N1-1)) / (x1u(N1)-x1u(N1-1))
     cu11R( 1,N1,g) = 0.
  END IF
  
  IF (BC_1L == 3 .OR. BC_1L == 4) THEN
     cu11R(-1,0 ,g) =  0.
     cu11R( 0,0 ,g) = -1./(x1u(1 )-x1u(0   ))
     cu11R( 1,0 ,g) =  1./(x1u(1 )-x1u(0   ))
  END IF
  IF (BC_1U == 3 .OR. BC_1U == 4) THEN
     cu11R(-1,N1,g) = -1./(x1u(N1)-x1u(N1-1))
     cu11R( 0,N1,g) =  1./(x1u(N1)-x1u(N1-1))
     cu11R( 1,N1,g) =  0.
  END IF
  !-----------------------------------------------------------------------------------------------------------
  CALL diff_coeffs(2,2,0,mapping_yes,2,(/2,3/),x2v(-1:N2+1),x2v(-1:N2+1),BC_2L,BC_2U,1,N2,-1,1,cv22R(-1,0,g))
  
  IF (BC_2L == 1 .OR. BC_2L == 2) THEN
     cv22R(-1,0 ,g) = 0.
     cv22R( 0,0 ,g) = 1.- (x2p(1 )-x2v(0   )) / (x2v(1 )-x2v(0   ))
     cv22R( 1,0 ,g) =     (x2p(1 )-x2v(0   )) / (x2v(1 )-x2v(0   ))
  END IF
  IF (BC_2U == 1 .OR. BC_2U == 2) THEN
     cv22R(-1,N2,g) = 1.- (x2p(N2)-x2v(N2-1)) / (x2v(N2)-x2v(N2-1))
     cv22R( 0,N2,g) =     (x2p(N2)-x2v(N2-1)) / (x2v(N2)-x2v(N2-1))
     cv22R( 1,N2,g) = 0.
  END IF
  
  IF (BC_2L == 3 .OR. BC_2L == 4) THEN
     cv22R(-1,0 ,g) =  0.
     cv22R( 0,0 ,g) = -1./(x2v(1 )-x2v(0   ))
     cv22R( 1,0 ,g) =  1./(x2v(1 )-x2v(0   ))
  END IF
  IF (BC_2U == 3 .OR. BC_2U == 4) THEN
     cv22R(-1,N2,g) = -1./(x2v(N2)-x2v(N2-1))
     cv22R( 0,N2,g) =  1./(x2v(N2)-x2v(N2-1))
     cv22R( 1,N2,g) =  0.
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
  
  CALL diff_coeffs(3,2,0,mapping_yes,2,(/2,3/),x3w(-1:N3+1),x3w(-1:N3+1),BC_3L,BC_3U,1,N3,-1,1,cw33R(-1,0,g))
  
  IF (BC_3L == 1 .OR. BC_3L == 2) THEN
     cw33R(-1,0 ,g) = 0.
     cw33R( 0,0 ,g) = 1.- (x3p(1 )-x3w(0   )) / (x3w(1 )-x3w(0   ))
     cw33R( 1,0 ,g) =     (x3p(1 )-x3w(0   )) / (x3w(1 )-x3w(0   ))
  END IF
  IF (BC_3U == 1 .OR. BC_3U == 2) THEN
     cw33R(-1,N3,g) = 1.- (x3p(N3)-x3w(N3-1)) / (x3w(N3)-x3w(N3-1))
     cw33R( 0,N3,g) =     (x3p(N3)-x3w(N3-1)) / (x3w(N3)-x3w(N3-1))
     cw33R( 1,N3,g) = 0.
  END IF
  
  IF (BC_3L == 3 .OR. BC_3L == 4) THEN
     cw33R(-1,0 ,g) =  0.
     cw33R( 0,0 ,g) = -1./(x3w(1 )-x3w(0   ))
     cw33R( 1,0 ,g) =  1./(x3w(1 )-x3w(0   ))
  END IF
  IF (BC_3U == 3 .OR. BC_3U == 4) THEN
     cw33R(-1,N3,g) = -1./(x3w(N3)-x3w(N3-1))
     cw33R( 0,N3,g) =  1./(x3w(N3)-x3w(N3-1))
     cw33R( 1,N3,g) =  0.
  END IF
  
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE get_stencil_Helm
  
  
  
  
  
  
  
  
  
  
  !> subroutine that interpolates the coefficient from a coarser to a finer grid. (pressure equation)
  SUBROUTINE interp_coeffs
  
  IMPLICIT NONE
  
  INTEGER               ::  i, i0, di
  INTEGER               ::  j, j0, dj
  INTEGER               ::  k, k0, dk
  INTEGER               ::  g
  
  REAL                  ::  Dx12, Dx10
  
  
  cI1 = 0.
  cI2 = 0.
  cI3 = 0.
  
  
  !===========================================================================================================
  !=== Interpolation, linienweise, 1d ========================================================================
  !===========================================================================================================
  DO g = 1, n_grids-1
     
     di = (M1-1)/((NN(1,g)-1)*NB(1,g))
     dj = (M2-1)/((NN(2,g)-1)*NB(2,g))
     dk = (M3-1)/((NN(3,g)-1)*NB(3,g))
     
     !--------------------------------------------------------------------------------------------------------
     DO i = 2, NN(1,g)-1, 2
        i0 = 1 + (i-1)*di + (M1-1)*(iB(1,g)-1)/NB(1,g)
        
        Dx10 = y1p(i0   )-y1p(i0-di)
        Dx12 = y1p(i0+di)-y1p(i0-di)
        
        cI1(1,i,g) = 1.- Dx10/Dx12
        cI1(2,i,g) =     Dx10/Dx12
     END DO
     !--------------------------------------------------------------------------------------------------------
     DO j = 2, NN(2,g)-1, 2
        j0 = 1 + (j-1)*dj + (M2-1)*(iB(2,g)-1)/NB(2,g)
        
        Dx10 = y2p(j0   )-y2p(j0-dj)
        Dx12 = y2p(j0+dj)-y2p(j0-dj)
        
        cI2(1,j,g) = 1.- Dx10/Dx12
        cI2(2,j,g) =     Dx10/Dx12
     END DO
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
     DO k = 2, NN(3,g)-1, 2
        k0 = 1 + (k-1)*dk + (M3-1)*(iB(3,g)-1)/NB(3,g)
        
        Dx10 = y3p(k0   )-y3p(k0-dk)
        Dx12 = y3p(k0+dk)-y3p(k0-dk)
        
        cI3(1,k,g) = 1.- Dx10/Dx12
        cI3(2,k,g) =     Dx10/Dx12
     END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END DO
  !===========================================================================================================
  
  
  END SUBROUTINE interp_coeffs
  
  
  
  
  
  
  
  
  
  
  !> subroutine that interpolates coefficients from a coarser to a finer grid. (velocity equation)
  !! Boudary conditions are extrapolated.
  SUBROUTINE interp_coeffs_Helm
  
  IMPLICIT NONE
  
  INTEGER               ::  i, j, k
  REAL                  ::  Dx12, Dx1a, Dx1b
  
  
  cIH1 = 0.
  cIH2 = 0.
  cIH3 = 0.
  
  !===========================================================================================================
  !=== Interpolation, linienweise, 1d ========================================================================
  !===========================================================================================================
  DO i = 1, N1-2, 2
     Dx1a = x1u(i  )-x1p(i)
     Dx1b = x1u(i+1)-x1p(i)
     Dx12 = x1p(i+2)-x1p(i)
     
     cIH1(1,i  ) = 1.- Dx1a/Dx12
     cIH1(2,i  ) =     Dx1a/Dx12
     
     cIH1(1,i+1) = 1.- Dx1b/Dx12
     cIH1(2,i+1) =     Dx1b/Dx12
  END DO
  
  !--- Randbedingungen (Extrapolation) ---
  IF (BC_1L .GT. 0) THEN
     Dx1a = x1u(0)-x1p(1)
     Dx12 = x1p(3)-x1p(1)
     
     cIH1(1,0) = 1.- Dx1a/Dx12
     cIH1(2,0) =     Dx1a/Dx12
  END IF
  IF (BC_1U .GT. 0) THEN
     Dx1a = x1u(N1)-x1p(N1-2)
     Dx12 = x1p(N1)-x1p(N1-2)
     
     cIH1(1,N1) = 1.- Dx1a/Dx12
     cIH1(2,N1) =     Dx1a/Dx12
  END IF
  !-----------------------------------------------------------------------------------------------------------
  DO j = 1, N2-2, 2
     Dx1a = x2v(j  )-x2p(j)
     Dx1b = x2v(j+1)-x2p(j)
     Dx12 = x2p(j+2)-x2p(j)
     
     cIH2(1,j  ) = 1.- Dx1a/Dx12
     cIH2(2,j  ) =     Dx1a/Dx12
     
     cIH2(1,j+1) = 1.- Dx1b/Dx12
     cIH2(2,j+1) =     Dx1b/Dx12
  END DO
  
  !--- Randbedingungen (Extrapolation) ---
  IF (BC_2L .GT. 0) THEN
     Dx1a = x2v(0)-x2p(1)
     Dx12 = x2p(3)-x2p(1)
     
     cIH2(1,0) = 1.- Dx1a/Dx12
     cIH2(2,0) =     Dx1a/Dx12
  END IF
  IF (BC_2U .GT. 0) THEN
     Dx1a = x2v(N2)-x2p(N2-2)
     Dx12 = x2p(N2)-x2p(N2-2)
     
     cIH2(1,N2) = 1.- Dx1a/Dx12
     cIH2(2,N2) =     Dx1a/Dx12
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
  
  DO k = 1, N3-2, 2
     Dx1a = x3w(k  )-x3p(k)
     Dx1b = x3w(k+1)-x3p(k)
     Dx12 = x3p(k+2)-x3p(k)
     
     cIH3(1,k  ) = 1.- Dx1a/Dx12
     cIH3(2,k  ) =     Dx1a/Dx12
     
     cIH3(1,k+1) = 1.- Dx1b/Dx12
     cIH3(2,k+1) =     Dx1b/Dx12
  END DO
  
  !--- Randbedingungen (Extrapolation) ---
  IF (BC_3L .GT. 0) THEN
     Dx1a = x3w(0)-x3p(1)
     Dx12 = x3p(3)-x3p(1)
     
     cIH3(1,0) = 1.- Dx1a/Dx12
     cIH3(2,0) =     Dx1a/Dx12
  END IF
  IF (BC_3U .GT. 0) THEN
     Dx1a = x3w(N3)-x3p(N3-2)
     Dx12 = x3p(N3)-x3p(N3-2)
     
     cIH3(1,N3) = 1.- Dx1a/Dx12
     cIH3(2,N3) =     Dx1a/Dx12
  END IF
  
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE interp_coeffs_Helm
  
  
  
  
  
  
  
  
  
  
  !> subroutine to restrict coefficients from a finer to a coarser grid (pressure equation).
  !! @todo: some assignments need testing.
  !! @todo: clean up and substitute variables.
  !! @todo: there is a flag (1==2) --> clarify and clean up.
  SUBROUTINE restr_coeffs ! TEST!!! aufraeumen und Variablen substituieren ...
  
  IMPLICIT NONE
  
  INTEGER               ::  g
  INTEGER               ::  i, ii, iimax
  INTEGER               ::  j, jj, jjmax
  INTEGER               ::  k, kk, kkmax
  
  INTEGER               ::  BC_1L, BC_1U, BC_2L, BC_2U, BC_3L, BC_3U ! TEST!!!
  
  
  cR1 = 0.
  cR2 = 0.
  cR3 = 0.
  
  !===========================================================================================================
  !=== Restriktion, linienweise, 1d ==========================================================================
  !===========================================================================================================
  DO g = 2, n_grids
     
     iimax = (NN(1,g)-1)/n_gather(1,g)+1
     jjmax = (NN(2,g)-1)/n_gather(2,g)+1
     kkmax = (NN(3,g)-1)/n_gather(3,g)+1
     
     !iimax = (NN(1,g)-1)*NB(1,g)/NB(1,g-1)+1 ! TEST!!! alt ...
     !jjmax = (NN(2,g)-1)*NB(2,g)/NB(2,g-1)+1
     !kkmax = (NN(3,g)-1)*NB(3,g)/NB(3,g-1)+1
     !--------------------------------------------------------------------------------------------------------
     IF (n_gather(1,g) .GT. 1) THEN
        BC_1L = BC(1,1,g-1) ! TEST!!! evtl. wieder auf dem feinen Gitter speichern ...
        BC_1U = BC(2,1,g-1)
     ELSE
        BC_1L = BC(1,1,g)
        BC_1U = BC(2,1,g)
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (n_gather(2,g) .GT. 1) THEN
        BC_2L = BC(1,2,g-1)
        BC_2U = BC(2,2,g-1)
     ELSE
        BC_2L = BC(1,2,g)
        BC_2U = BC(2,2,g)
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (n_gather(3,g) .GT. 1) THEN
        BC_3L = BC(1,3,g-1)
        BC_3U = BC(2,3,g-1)
     ELSE
        BC_3L = BC(1,3,g)
        BC_3U = BC(2,3,g)
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     DO ii = 1, iimax
        cR1(-1,ii,g) = 1./4.
        cR1( 0,ii,g) = 2./4.
        cR1( 1,ii,g) = 1./4.
     END DO
     
     IF (BC_1L .GT. 0) THEN
        cR1(-1,1,g) = 0.
        cR1( 0,1,g) = 1.
        cR1( 1,1,g) = 0.
        
        cR1(-1,2,g) = 0. ! TEST!!! Sollte evtl. noch ergaenzt werden ...
        cR1( 0,2,g) = 1.
        cR1( 1,2,g) = 0.
     END IF
     IF (BC_1L == -2) THEN
        cR1( 1,1,g) = cR1( 1,1,g) + cR1(-1,1,g)
        cR1(-1,1,g) = 0.
     END IF
     
     IF (BC_1U .GT. 0) THEN
        cR1(-1,iimax  ,g) = 0.
        cR1( 0,iimax  ,g) = 1.
        cR1( 1,iimax  ,g) = 0.
        
        cR1(-1,iimax-1,g) = 0. ! TEST!!!
        cR1( 0,iimax-1,g) = 1.
        cR1( 1,iimax-1,g) = 0.
     END IF
     IF (BC_1U == -2) THEN
        cR1(-1,iimax,g) = cR1( 1,iimax,g) + cR1(-1,iimax,g)
        cR1( 1,iimax,g) = 0.
     END IF
     !--------------------------------------------------------------------------------------------------------
     DO jj = 1, jjmax
        cR2(-1,jj,g) = 1./4.
        cR2( 0,jj,g) = 2./4.
        cR2( 1,jj,g) = 1./4.
     END DO
     
     IF (BC_2L .GT. 0) THEN
        cR2(-1,1,g) = 0.
        cR2( 0,1,g) = 1.
        cR2( 1,1,g) = 0.
        
        cR2(-1,2,g) = 0. ! TEST!!!
        cR2( 0,2,g) = 1.
        cR2( 1,2,g) = 0.
     END IF
     IF (BC_2L == -2) THEN
        cR2( 1,1,g) = cR2( 1,1,g) + cR2(-1,1,g)
        cR2(-1,1,g) = 0.
     END IF
     
     IF (BC_2U .GT. 0) THEN
        cR2(-1,jjmax  ,g) = 0.
        cR2( 0,jjmax  ,g) = 1.
        cR2( 1,jjmax  ,g) = 0.
        
        cR2(-1,jjmax-1,g) = 0. ! TEST!!!
        cR2( 0,jjmax-1,g) = 1.
        cR2( 1,jjmax-1,g) = 0.
     END IF
     IF (BC_2U == -2) THEN
        cR2(-1,jjmax,g) = cR2( 1,jjmax,g) + cR2(-1,jjmax,g)
        cR2( 1,jjmax,g) = 0.
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
     DO kk = 1, kkmax
        cR3(-1,kk,g) = 1./4.
        cR3( 0,kk,g) = 2./4.
        cR3( 1,kk,g) = 1./4.
     END DO
     
     IF (BC_3L .GT. 0) THEN
        cR3(-1,1,g) = 0.
        cR3( 0,1,g) = 1.
        cR3( 1,1,g) = 0.
        
        cR3(-1,2,g) = 0. ! TEST!!!
        cR3( 0,2,g) = 1.
        cR3( 1,2,g) = 0.
     END IF
     IF (BC_3L == -2) THEN
        cR3( 1,1,g) = cR3( 1,1,g) + cR3(-1,1,g)
        cR3(-1,1,g) = 0.
     END IF
     
     IF (BC_3U .GT. 0) THEN
        cR3(-1,kkmax  ,g) = 0.
        cR3( 0,kkmax  ,g) = 1.
        cR3( 1,kkmax  ,g) = 0.
        
        cR3(-1,kkmax-1,g) = 0. ! TEST!!!
        cR3( 0,kkmax-1,g) = 1.
        cR3( 1,kkmax-1,g) = 0.
     END IF
     IF (BC_3U == -2) THEN
        cR3(-1,kkmax,g) = cR3( 1,kkmax,g) + cR3(-1,kkmax,g)
        cR3( 1,kkmax,g) = 0.
     END IF
     
     END IF
     !--------------------------------------------------------------------------------------------------------
  END DO
  !===========================================================================================================
  
  
  IF (1 == 2) THEN ! TEST!!!
  !===========================================================================================================
  !=== Restriktion, linienweise, 1d ==========================================================================
  !===========================================================================================================
  cRest1 = 0.
  cRest2 = 0.
  cRest3 = 0.
  
  DO g = 1, n_grids-1
     iimax = NN(1,g)
     jjmax = NN(2,g)
     kkmax = NN(3,g)
     
     IF (ncb1f(dim_ncb1c) .LT. iimax) THEN ! TEST!!! Quick and dirty ...
        CALL diff_coeffs(1,0,0,.FALSE.,dim_ncb1c,ncb1r,x1pR(b1L,g),x1pR(b1L,g),BC_1L,BC_1U,0,iimax,b1L,b1U,cRest1(b1L,0,g))
     ELSE
        DO i = 1, iimax
           cRest1(-1:1,i,g) = (/0.25,0.5,0.25/) ! TEST!!! angepasste Koeffizienten? (mapping_yes = .FALSE. waere besser ...)
        END DO
        IF (BC_1L .GT. 0) cRest1(-1:1,1    ,g) = (/0. ,1. ,0. /)
        IF (BC_1L ==  -2) cRest1(-1:1,1    ,g) = (/0. ,0.5,0.5/)
        IF (BC_1U .GT. 0) cRest1(-1:1,iimax,g) = (/0. ,1. ,0. /)
        IF (BC_1U ==  -2) cRest1(-1:1,iimax,g) = (/0.5,0.5,0. /)
     END IF
     
     
     IF (ncb2f(dim_ncb2c) .LT. jjmax) THEN
        CALL diff_coeffs(2,0,0,.FALSE.,dim_ncb2c,ncb2r,x2pR(b2L,g),x2pR(b2L,g),BC_2L,BC_2U,0,jjmax,b2L,b2U,cRest2(b2L,0,g))
     ELSE
        DO j = 1, jjmax
           cRest2(-1:1,j,g) = (/0.25,0.5,0.25/)
        END DO
        IF (BC_2L .GT. 0) cRest2(-1:1,1    ,g) = (/0. ,1. ,0. /)
        IF (BC_2L ==  -2) cRest2(-1:1,1    ,g) = (/0. ,0.5,0.5/)
        IF (BC_2U .GT. 0) cRest2(-1:1,jjmax,g) = (/0. ,1. ,0. /)
        IF (BC_2U ==  -2) cRest2(-1:1,jjmax,g) = (/0.5,0.5,0. /)
     END IF
     
     
     IF (dimens == 3) THEN
     IF (ncb3f(dim_ncb3c) .LT. kkmax) THEN
        CALL diff_coeffs(3,0,0,.FALSE.,dim_ncb3c,ncb3r,x3pR(b3L,g),x3pR(b3L,g),BC_3L,BC_3U,0,kkmax,b3L,b3U,cRest3(b3L,0,g))
     ELSE
        DO k = 1, kkmax
           cRest3(-1:1,k,g) = (/0.25,0.5,0.25/)
        END DO
        IF (BC_3L .GT. 0) cRest3(-1:1,1    ,g) = (/0. ,1. ,0. /)
        IF (BC_3L ==  -2) cRest3(-1:1,1    ,g) = (/0. ,0.5,0.5/)
        IF (BC_3U .GT. 0) cRest3(-1:1,kkmax,g) = (/0. ,1. ,0. /)
        IF (BC_3U ==  -2) cRest3(-1:1,kkmax,g) = (/0.5,0.5,0. /)
     END IF
     END IF
     
  END DO
  !===========================================================================================================
  END IF
  
  
  END SUBROUTINE restr_coeffs
  
  
  
  
  
  
  
  
  
  
  !> subroutine to restrict the coefficients from a finer to a coarser grid (velocity equation).
  SUBROUTINE restr_coeffs_Helm
  
  IMPLICIT NONE
  
  INTEGER               ::  i, j, k
  REAL                  ::  Dx12, Dx1a
  
  
  cRH1 = 0.
  cRH2 = 0.
  cRH3 = 0.
  
  !===========================================================================================================
  !=== Restriktion, linienweise, 1d ==========================================================================
  !===========================================================================================================
  DO i = 1, N1, 2
     Dx1a = x1p(i)-x1u(i-1)
     Dx12 = x1u(i)-x1u(i-1)
     
     cRH1(1,i) = 1.- Dx1a/Dx12
     cRH1(2,i) =     Dx1a/Dx12
  END DO
  
  IF (BC_1L .GT. 0) THEN
     cRH1(1,1 ) = 1.
     cRH1(2,1 ) = 0.
  END IF
  IF (BC_1L == -2) THEN
     cRH1(:,1 ) = 0.
  END IF
  
  IF (BC_1U .GT. 0) THEN
     cRH1(1,N1) = 0.
     cRH1(2,N1) = 1.
  END IF
  IF (BC_1U == -2) THEN
     cRH1(:,N1) = 0.
  END IF
  !-----------------------------------------------------------------------------------------------------------
  DO j = 1, N2, 2
     Dx1a = x2p(j)-x2v(j-1)
     Dx12 = x2v(j)-x2v(j-1)
     
     cRH2(1,j) = 1.- Dx1a/Dx12
     cRH2(2,j) =     Dx1a/Dx12
  END DO
  
  IF (BC_2L .GT. 0) THEN
     cRH2(1,1 ) = 1.
     cRH2(2,1 ) = 0.
  END IF
  IF (BC_2L == -2) THEN
     cRH2(:,1 ) = 1.
  END IF
  
  IF (BC_2U .GT. 0) THEN
     cRH2(1,N2) = 0.
     cRH2(2,N2) = 1.
  END IF
  IF (BC_2U == -2) THEN
     cRH2(:,N2) = 0.
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
  
  DO k = 1, N3, 2
     Dx1a = x3p(k)-x3w(k-1)
     Dx12 = x3w(k)-x3w(k-1)
     
     cRH3(1,k) = 1.- Dx1a/Dx12
     cRH3(2,k) =     Dx1a/Dx12
  END DO
  
  IF (BC_3L .GT. 0) THEN
     cRH3(1,1 ) = 1.
     cRH3(2,1 ) = 0.
  END IF
  IF (BC_3L == -2) THEN
     cRH3(:,1 ) = 1.
  END IF
  
  IF (BC_3U .GT. 0) THEN
     cRH3(1,N3) = 0.
     cRH3(2,N3) = 1.
  END IF
  IF (BC_3U == -2) THEN
     cRH3(:,N3) = 0.
  END IF
  
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE restr_coeffs_Helm
  
  
  
  
  
  
  
  
  
  
  !> subroutine for computing the weights for checking of divergence-freeness. 
  !! @todo: only concerning explicit differences so far.
  !! @todo: borders need separate treatment.
  SUBROUTINE get_weights ! TEST!!! bezieht sich bislang nur auf explizite Differenzen ... ODER: rauswerfen!
  
  ! revised: 24.10.07
  
  IMPLICIT NONE
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Diese Abschaetzung ist im Feld etwas konservativer (Faktor 1.5-2.5) im Vergleich zu       !
  !                Differenzen 2. Konvergenzordnung. Am Rand dagegen sehr viel strikter!                     !
  !              - dx_iL, dx_iU dürfen nicht der Ableitung der Mapping-Funktion an der Wand gleichgesetzt    !
  !                werden, da diese im Extremfall auch verschwinden kann.                                    !
  !----------------------------------------------------------------------------------------------------------!
  
  
  
  IF (1 == 1) THEN ! TEST!!! Rand muesste noch extra behandelt werden!
  !===========================================================================================================
  !=== Gewichte für Überprüfung der Divergenzfreiheit ========================================================
  !===========================================================================================================
  weight = 0.
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
     
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              DO ii = d1L, d1U
                 weight(i,j,k) = weight(i,j,k) + ABS(cDu1(ii,i))
              END DO
              DO jj = d2L, d2U
                 weight(i,j,k) = weight(i,j,k) + ABS(cDv2(jj,j))
              END DO
              DO kk = d3L, d3U
                 weight(i,j,k) = weight(i,j,k) + ABS(cDw3(kk,k))
              END DO
              weight(i,j,k) = 1./weight(i,j,k)
           END DO
        END DO
     END DO
     
  !-----------------------------------------------------------------------------------------------------------
  ELSE
     
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              DO ii = d1L, d1U
                 weight(i,j,k) = weight(i,j,k) + ABS(cDu1(ii,i))
              END DO
              DO jj = d2L, d2U
                 weight(i,j,k) = weight(i,j,k) + ABS(cDv2(jj,j))
              END DO
              weight(i,j,k) = 1./weight(i,j,k)
           END DO
        END DO
     END DO
     
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  ELSE
     weight = 1.
  END IF
  
  IF (1 == 1) THEN ! TEST!!!
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_1L .GT. 0) THEN
     i = 1
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO ii = d1L, d1U
              weight(i,j,k) = weight(i,j,k) + ABS(cDu1(ii,i))
           END DO
           weight(i,j,k) = 1./weight(i,j,k)
        END DO
     END DO
  END IF
  IF (BC_1U .GT. 0) THEN
     i = N1
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO ii = d1L, d1U
              weight(i,j,k) = weight(i,j,k) + ABS(cDu1(ii,i))
           END DO
           weight(i,j,k) = 1./weight(i,j,k)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_2L .GT. 0) THEN
     j = 1
     DO k = S3p, N3p
        DO i = S1p, N1p
           DO jj = d2L, d2U
              weight(i,j,k) = weight(i,j,k) + ABS(cDv2(jj,j))
           END DO
           weight(i,j,k) = 1./weight(i,j,k)
        END DO
     END DO
  END IF
  IF (BC_2U .GT. 0) THEN
     j = N2
     DO k = S3p, N3p
        DO i = S1p, N1p
           DO jj = d2L, d2U
              weight(i,j,k) = weight(i,j,k) + ABS(cDv2(jj,j))
           END DO
           weight(i,j,k) = 1./weight(i,j,k)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_3L .GT. 0) THEN
     k = 1
     DO j = S2p, N2p
        DO i = S1p, N1p
           DO kk = d3L, d3U
              weight(i,j,k) = weight(i,j,k) + ABS(cDw3(kk,k))
           END DO
           weight(i,j,k) = 1./weight(i,j,k)
        END DO
     END DO
  END IF
  IF (BC_3U .GT. 0) THEN
     k = N3
     DO j = S2p, N2p
        DO i = S1p, N1p
           DO kk = d3L, d3U
              weight(i,j,k) = weight(i,j,k) + ABS(cDw3(kk,k))
           END DO
           weight(i,j,k) = 1./weight(i,j,k)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  END IF
  
  
  ! Kanten/Ecken sind durch die Randbedingungen immer Divergenz-frei:
  IF (BC_1L .GT. 0 .AND. BC_2L .GT. 0) weight(1 ,1 ,1:N3) = 0.
  IF (BC_1L .GT. 0 .AND. BC_2U .GT. 0) weight(1 ,N2,1:N3) = 0.
  IF (BC_1U .GT. 0 .AND. BC_2L .GT. 0) weight(N1,1 ,1:N3) = 0.
  IF (BC_1U .GT. 0 .AND. BC_2U .GT. 0) weight(N1,N2,1:N3) = 0.
  
  IF (BC_1L .GT. 0 .AND. BC_3L .GT. 0) weight(1 ,1:N2,1 ) = 0.
  IF (BC_1L .GT. 0 .AND. BC_3U .GT. 0) weight(1 ,1:N2,N3) = 0.
  IF (BC_1U .GT. 0 .AND. BC_3L .GT. 0) weight(N1,1:N2,1 ) = 0.
  IF (BC_1U .GT. 0 .AND. BC_3U .GT. 0) weight(N1,1:N2,N3) = 0.
  
  IF (BC_2L .GT. 0 .AND. BC_3L .GT. 0) weight(1:N1,1 ,1 ) = 0.
  IF (BC_2L .GT. 0 .AND. BC_3U .GT. 0) weight(1:N1,1 ,N3) = 0.
  IF (BC_2U .GT. 0 .AND. BC_3L .GT. 0) weight(1:N1,N2,1 ) = 0.
  IF (BC_2U .GT. 0 .AND. BC_3U .GT. 0) weight(1:N1,N2,N3) = 0.
  !===========================================================================================================
    
  END SUBROUTINE get_weights
  
  
  
  
  
  
  
  
  
  
  !< subroutine for inverting a matrix. Only used for FD coefficients.
  SUBROUTINE Matrix_invert(N,matrix,matrix_inv)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)  ::  N                    !< number of coefficients in stencil
  REAL   , INTENT(in)  ::  matrix     (1:N,1:N) !< original matrix
  REAL   , INTENT(out) ::  matrix_inv (1:N,1:N) !< inverted matrix
  
  REAL                 ::  matrix_left(1:N,1:N)
  REAL                 ::  mult1, mult2
  REAL                 ::  eps
  INTEGER              ::  i, j, k
  
  REAL                 ::  store
  INTEGER              ::  maxValue, maxValuePos
  
  
  
  eps = 10.**(-20) ! double precision
  !eps = 10.**(-??) ! single precision
  
  matrix_left = matrix
  
  
  ! Rechte Seite initialisieren (= Einheitsmatrix):
  matrix_inv = 0.
  
  DO i = 1, N
     matrix_inv(i,i) = 1.
  END DO  
  
  
  !--- Vorwaertsschritt -------------------------------------
  ! linke Seite umformen in obere Dreiecksmatrix
  ! rechte Seite umformen in untere Dreiecksmatrix
  
  DO j = 1, N-1
     
     ! Pivoting == Umsortieren der aktuellen Untermatrix (j:N,j:N)
     ! (Diagonalelemente der linken Dreiecksmatrix sind betragsmaessig zu maximieren)
     ! Groesster Wert in Spalte j = zukuenftiges Diagonalelement j,j
     maxValue    = ABS(matrix_left(j,j))
     maxValuePos = j
     DO i = j+1, N
        IF (ABS(matrix_left(i,j)) .GT. maxValue) THEN
           maxValue    = ABS(matrix_left(i,j))
           maxValuePos = i
        END IF
     END DO
     
     ! Zeilen vertauschen:
     IF (maxValuePos /= j) THEN
        DO i = 1, N
           store                      = matrix_left(maxValuePos,i)
           matrix_left(maxValuePos,i) = matrix_left(j,i)
           matrix_left(j,i)           = store
           
           store                      = matrix_inv(maxValuePos,i)
           matrix_inv(maxValuePos,i)  = matrix_inv(j,i)
           matrix_inv(j,i)            = store
        END DO
     END IF
     
     mult1 = matrix_left(j,j)
     IF (ABS(mult1) .LT. eps .AND. rank == 0) THEN
        WRITE(*,'(a)'      ) 'WARNING! Matrix is probably singular (1)!'
        WRITE(*,'(a,E13.5)') '       ... dividing by', mult1
        WRITE(*,'(a,i2)'   ) '                   i =', j
     END IF
     
     DO i = j+1, N
        mult2 = matrix_left(i,j)/mult1
                
        ! Aufaddieren von Zeile j auf aktuelle Zeile i (linke Seite):
        DO k = j, N
           ! To avoid underflow:
           IF (.NOT. (ABS(mult2) .LT. eps .AND. ABS(matrix_left(j,k)) .LT. eps)) THEN
              matrix_left(i,k) = matrix_left(i,k) - matrix_left(j,k)*mult2
           END IF
        END DO
        
        ! Aufaddieren von Zeile j auf aktuelle Zeile i (rechte Seite):
        DO k = 1, N
           ! To avoid underflow:
           IF (.NOT. (ABS(mult2) .LT. eps .AND. ABS(matrix_inv(j,k)) .LT. eps)) THEN
              matrix_inv(i,k) = matrix_inv(i,k) - matrix_inv(j,k)*mult2
           END IF
        END DO
        
        ! Komponente i,j explizit zu Null setzen:
        matrix_left(i,j) = 0.
     END DO
     
  END DO
  
  
  !--- Rueckwaertsschritt -----------------------------------
  ! linke Seite umformen in Einheitsmatrix
  ! rechte Seite umformen in gesuchte inverse Matrix
  
  DO j = N, 2, -1
     
     ! Multiplikator:
     mult1 = matrix_left(j,j)
     IF (ABS(mult1) .LT. eps .AND. rank == 0) THEN
        WRITE(*,'(a)'      ) 'WARNING! Matrix is probably singular (2)!'
        WRITE(*,'(a,E13.5)') '       ... dividing by', mult1
        WRITE(*,'(a,i2)'   ) '                   i =', j
     END IF
     
     DO i = 1, j-1
        mult2 = matrix_left(i,j)/mult1
        
        ! Aufaddieren von Zeile j auf aktuelle Zeile i (rechte Seite):
        DO k = 1, N
           ! To avoid underflow:
           IF (.NOT. (ABS(mult2) .LT. eps .AND. ABS(matrix_inv(j,k)) .LT. eps)) THEN
              matrix_inv(i,k) = matrix_inv(i,k) - matrix_inv(j,k)*mult2
           END IF
        END DO
        
        ! nicht-Diagonalelement explizit zu 0 setzen:
        matrix_left(i,j) = 0.
     END DO
     
  END DO
  
  
  ! linke Matrix auf Einheitsmatrix bringen:
  DO i = 1, N
     mult1 = matrix_left(i,i)
     IF (ABS(mult1) .LT. eps .AND. rank == 0) THEN
        WRITE(*,'(a)'      ) 'WARNING! Matrix is probably singular (3)!'
        WRITE(*,'(a,E13.5)') '       ... dividing by', mult1
        WRITE(*,'(a,i2)'   ) '                   i =', i
     END IF
     matrix_inv(i,:) = matrix_inv(i,:)/mult1
  END DO
  
  
  END SUBROUTINE Matrix_invert
  
  
  
  
END MODULE mod_coeffs
