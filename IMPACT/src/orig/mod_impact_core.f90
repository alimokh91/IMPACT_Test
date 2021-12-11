MODULE mod_impact_core

IMPLICIT NONE

!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* May 2005 - Dec 2011                                                                                       *
!* by Michael John, Institute of Fluid Dynamics, ETH Zurich (john@ifd.mavt.ethz.ch)                          *
!* May 2013 - May 2013                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)           *
!*************************************************************************************************************

!> @file mod_impact_core.f90
!! Main file.

!> Main routine of IMPACT. All of the initialization is done here. Depending on the parameter task 
!! either the DNS is 
!! started, the postprocessing is called, an eigenvalue problem is solved (this feature corresponding to
!! task=3 has been migrated to a usr file by mjohn) or matrices are analyzed. The parameter is specified in
!! usr_config.f90

PRIVATE

PUBLIC :: impact_core_init
PUBLIC :: impact_core_finalize


CONTAINS


subroutine impact_core_init
  USE mod_vars
  USE mod_coeffs
  USE mod_geometry
  USE mod_setup
  USE mod_test
  USE mod_lib
  
  USE HDF5
  USE MPI

  IMPLICIT NONE

  !----------------------------------------------------------------------------------------------------------!
  ! Annotations / TO-DO list:                                                                                !
  !  - henniger 2011: unite the four init_* (general, parallel, boundaries, limits), already intermingled    !
  !  - mjohn 060513: moved task 3 (eigensolver) to a usr file because mod functions depended themselves on   !
  !                  usr functions, i.e. variables / fields of task 3 were not independent! May be recovered !
  !                  as part of a mod file if new fields are provided and dependence on usr functions lifted !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== Initialization ========================================================================================
  !===========================================================================================================
  !--- Initialize MPI ----------------------------------------------------------------------------------------
  CALL MPI_INIT(merror)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,merror)
  
  !--- Initialize HDF5 ---------------------------------------------------------------------------------------
  CALL h5open_f(herror)
  
  !--- Set alarm if queue is being used ----------------------------------------------------------------------
  CALL init_alarm

  !--- Test of input parameters ------------------------------------------------------------------------------
  CALL test_parameter
  
   !--- General -------------------------------------------------------------------------------------------
   CALL init_general
   !--- MPI -----------------------------------------------------------------------------------------------
   CALL init_parallel
   !--- Type of boundary conditions -----------------------------------------------------------------------
   CALL init_boundaries
   !--- Limits of indices ---------------------------------------------------------------------------------
   CALL init_limits
  
  !--- Physical coordinates ----------------------------------------------------------------------------------
  CALL coordinates

  IF (task .EQ. 5) THEN
  write(*,*) 'x1p(s1p) =', x1p(s1p), ', x1p(n1p) =', x1p(n1p), ', s1p = ', s1p,', n1p = ', n1p
  write(*,*) 'y1p(1) =', y1p(1), ', y1p(m1) =', y1p(m1), ', m1 = ', m1
  write(*,*) 's3p = ', s3p,', n3p = ', n3p
  write(*,*) 's33 = ', s33,', n33 = ', n33
  write(*,*) 's33b = ', s33b,', n33b = ', n33b
  write(*,*) 'n3 = ', n3
  write(*,*) 'x3w(s33) = ', x3w(s33),', x3w(n33) = ', x3w(n33)
  write(*,*) 'x3w(s33b) = ', x3w(s33b+1),', x3w(n33b) = ', x3w(n33b)
  write(*,*) 'm1 = ', m1,'m2 = ', m2 ,'m3 = ', m3
  write(*,*) 'ib(1,1) = ', iB(1,1), ' n1 = ', n1 ,';'
  write(*,*) 'ib(2,1) = ', iB(2,1), ' n2 = ', n2 ,';'
  write(*,*) 'ib(3,1) = ', iB(3,1), ' n3 = ', n3 ,';'
  END IF

  !--- reference velocity ------------------------------------------------------------------------------------
  Re = rho_fluid*U_ref*L_ref/mu_fluid
  
  !--- Determine differential coefficients -------------------------------------------------------------------
  CALL FD_coeffs
  
  !--- Get stencil of operator div(grad( )), order of convergence is 2 ---------------------------------------
  CALL get_stencil
  CALL get_stencil_Helm
  
  !--- Get interpolation coefficients (multigrid, order of convergence is 2) ---------------------------------
  CALL interp_coeffs
  CALL interp_coeffs_Helm
  CALL restr_coeffs
  CALL restr_coeffs_Helm
  
  !--- Differenzenkoeffizienten testen -----------------------------------------------------------------------
  ! bbecsek 031214: TODO: according to usr_config these routines are deprecated.
  IF (write_test_yes) CALL test_coeffs
  
  !--- Weights for stop criterion ----------------------------------------------------------------------------
  CALL get_weights
  
  !--- Beta (for analysis) -----------------------------------------------------------------------------------
  CALL get_beta
  !===========================================================================================================
end subroutine impact_core_init


subroutine impact_core_finalize
  USE mod_vars
  USE HDF5
  USE MPI

  IMPLICIT NONE

  !===========================================================================================================
  IF (rank == 0) THEN
     IF (write_stout_yes) WRITE(*,*)
     IF (write_stout_yes) WRITE(*,*) 'DONE ...'
     IF (write_stout_yes) WRITE(*,*)
  END IF
  
  !---- close HDF5 -------------------------------------------------------------------------------------------
  CALL h5close_f(herror)
  
  !---- close MPI --------------------------------------------------------------------------------------------
  CALL MPI_FINALIZE(merror)
  !===========================================================================================================
  
end subroutine impact_core_finalize

END MODULE mod_impact_core

