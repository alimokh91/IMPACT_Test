!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!*************************************************************************************************************
 
!> @file usr_boundcond.f90
!! file containing subroutines for velocity boundary conditions


!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  SUBROUTINE boundary_vel_stat
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE usr_vars
  USE usr_func
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k
  REAL                   :: parab
  
  
  !--- boundary conditions for velocity ---
  ! note: - advective boundary conditions are automatically initialized from the initial condition
  !       - if constant, bc11 etc. can also be specified in "initial_conditions_vel" in file
  !         "usr_initcond.f90"
  !       - cf. sketch in file "usr_geometry.f90"
  !
  !   velocity component
  !   |orientation of the boundary normal
  !   ||     grid points on the boundary
  !   ||     |         |      lower/upper boundary
  !   ||     |         |      |
  ! bc11(S21B:N21B,S31B:N31B,1:2)
  ! bc12(S11B:N11B,S31B:N31B,1:2)
  ! bc13(S11B:N11B,S21B:N21B,1:2)
  ! bc21(S22B:N22B,S32B:N32B,1:2)
  ! bc22(S12B:N12B,S32B:N32B,1:2)
  ! bc23(S12B:N12B,S22B:N22B,1:2)
  ! bc31(S23B:N23B,S33B:N33B,1:2)
  ! bc32(S13B:N13B,S33B:N33B,1:2)
  ! bc33(S13B:N13B,S23B:N23B,1:2)
  !
  !
  ! processor-block boundary condition types:
  !                           _
  ! symmetry BC:    BC = -2    |
  ! periodicity BC: BC = -1    |- symmetric, central FD stencils
  ! neighbor block: BC =  0   _|
  ! Dirichlet BC:   BC =  1    |
  ! Neumann BC*:    BC =  2    |- non-symmetric, skew FD stencils
  ! Robin BC*:      BC =  3   _|
  ! *not yet implemented for velocity
  !
  !    orientation of boundary normal
  !    |lower/upper boundary
  !    ||
  ! BC_1L
  !
  !IF (1==2) THEN
  !!--- set inlet to poiseuille profile ---
  !DO k = S31B, N31B
  !  DO j = S21B, N21B
  !    CALL poiseuille_parabola( x2p(S2p) , x2p(N2p) , x2p(j) , parab)

  !    bc11(j,k,1) = parab*smooth_step(.TRUE.,time/10.)
  !  END DO
  !END DO

  !END IF
  
  
  !bc21(S22B:N22B,S32B:N32B,1) = 0.
  !bc31(S23B:N23B,S33B:N33B,1) = 0.


  
  !--- set upper and lower wall BCs last to ensure they are set in corners! ---
  !--- 2 direction
  bc12(S11B:N11B,S31B:N31B,1:2) = 0.
  bc22(S12B:N12B,S32B:N32B,1:2) = 0.
  bc32(S13B:N13B,S33B:N33B,1:2) = 0.

  !bc11(S21B:N21B,S31B:N31B,1) = 1.

  
  !IF ((.NOT. outlet(2,1,2)) .AND. BC_2L == 1) THEN ! Note: this case differentiation may save computational time but is not stringently necessary
  !   DO k = S32B, N32B
  !      DO i = S12B, N12B
  !         bc22(i,k,1) = interface((x1p(i)-0.75*L1)/(0.01*L1))
  !      END DO
  !   END DO
  !END IF
  
  
  END SUBROUTINE boundary_vel_stat
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE boundary_vel_tint
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE usr_vars
  USE usr_func
  
  IMPLICIT NONE
  
  
  !--- additional terms on RHS of time-integrated velocity boundary conditions ---
  ! note: - du/dt = RHS-nlbc
  !       - if not specified, advective boundary conditions are used
  !       - cf. sketch in file "usr_geometry.f90"
  !
  !     velocity component
  !     |orientation of the boundary normal
  !     ||     grid points on the boundary
  !     ||     |         |      lower/upper boundary
  !     ||     |         |      |
  ! nlbc11(S21B:N21B,S31B:N31B,1:2)
  ! nlbc12(S11B:N11B,S31B:N31B,1:2)
  ! nlbc13(S11B:N11B,S21B:N21B,1:2)
  ! nlbc21(S22B:N22B,S32B:N32B,1:2)
  ! nlbc22(S12B:N12B,S32B:N32B,1:2)
  ! nlbc23(S12B:N12B,S22B:N22B,1:2)
  ! nlbc31(S23B:N23B,S33B:N33B,1:2)
  ! nlbc32(S13B:N13B,S33B:N33B,1:2)
  ! nlbc33(S13B:N13B,S23B:N23B,1:2)
  !
  !
  ! processor-block boundary condition types:
  !                           _
  ! symmetry BC:    BC = -2    |
  ! periodicity BC: BC = -1    |- symmetric, central FD stencils
  ! neighbor block: BC =  0   _|
  ! Dirichlet BC:   BC =  1    |
  ! Neumann BC*:    BC =  2    |- non-symmetric, skew FD stencils
  ! Robin BC*:      BC =  3   _|
  ! *not yet implemented for velocity
  !
  !    orientation of boundary normal
  !    |lower/upper boundary
  !    ||
  ! BC_1L
  !
  
  
  END SUBROUTINE boundary_vel_tint
