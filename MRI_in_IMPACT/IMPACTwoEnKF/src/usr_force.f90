!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@arotrg.unibe.ch)                   *
!* October 2014                                                                                              *
!*************************************************************************************************************
 
!> @file usr_force.f90
!! file containing additional forcing terms for governing equaitons

!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> volume forcing of the momentum equation
  !! @note This could serve as the (most direct/straight-forward) interface between structural and fluid
  !!       solvers.
  SUBROUTINE forcing_vel
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE usr_vars
  USE usr_func
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k, m
  INTEGER                ::  ii, jj, kk
  REAL                   ::  lamb_fringe,parab,time_fact
  
  
  !--- additional volume forces in the momentum equation ---
  ! note: - du/dt = RHS-nl
  !       - cf. sketch in file "usr_geometry.f90"
  !
  !       grid points in the domain
  !       |       |       |    velocity component
  !       |       |       |    |
  ! nl(S11:N11,S21:N21,S31:N31,1)
  ! nl(S12:N12,S22:N22,S32:N32,2)
  ! nl(S13:N13,S23:N23,S33:N33,3)
  !

  ! fd is obtained from MOOSE
  !fd = 0.
  !nl = nl + fd * L_ref/(rho_fluid*U_ref**2.)

  !=== FRINGE FORCING =======================================================================================
  IF (fringe_yes) THEN
    CALL apply_fringe_forcing
  END IF
  !==========================================================================================================

  !=== WINDKESSEL LOADING ===================================================================================
  !IF (WK_yes) THEN
  !  CALL windkessel_integration_step
  !  CALL apply_windkessel_loading
  !END IF
  !==========================================================================================================

  !*** debugging
  !write(*,*)'x1p(65)=',x1p(65),' x1u(64)=',x1u(64)

  !DO k = S31, N31
  !   DO j = S21, N21
  !      DO i = S11, N11
  !         nl(i,j,k,1) = nl(i,j,k,1) - 1.*(0. - vel(i,j,k,1))*interface((x2p(j)-0.8*L2)/(0.5*L1))
  !      END DO
  !   END DO
  !END DO
  
  
  END SUBROUTINE forcing_vel
