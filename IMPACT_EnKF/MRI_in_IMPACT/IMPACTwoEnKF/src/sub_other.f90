!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!*************************************************************************************************************
  
!> @file sub_other.f90
!! file containing pseudocalls that should help the compiler to not mistake calls.

!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> pseudo call for mod_exchange::exchange. 
  SUBROUTINE pseudocall(phi)
  
  IMPLICIT NONE
  
  REAL   , INTENT(in)    ::  phi
  
  !---------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Pseudo-CALL fuer mod_exchange. Darf auf keinen Fall vom Compiler geinlined werden!                   !
  !---------------------------------------------------------------------------------------------------------------------!
  
  
  END SUBROUTINE pseudocall
  
  
  
  
  
  
  
  
  
  
  !> pseudo_call for mod_diff:apply_compact.
  SUBROUTINE pseudocall_int(phi)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)    ::  phi
  
  !---------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Pseudo-CALL fuer mod_exchange. Darf auf keinen Fall vom Compiler geinlined werden!                   !
  ! Update: bbecsek 031214: This is actually only called by mod_diff::apply_compact and not by mod_exchange::exchange!  !
  !                         If compact differences are removed this routine will become obsolete.                        !
  !---------------------------------------------------------------------------------------------------------------------!
  
  
  END SUBROUTINE pseudocall_int
  
  
  
  
