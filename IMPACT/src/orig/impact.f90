!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* May 2005 - Dec 2011                                                                                       *
!* by Michael John, Institute of Fluid Dynamics, ETH Zurich (john@ifd.mavt.ethz.ch)                          *
!* May 2013 - May 2013                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)		     *
!*************************************************************************************************************

!> @file impact.f90
!! Main file.

!> Main routine of IMPACT. All of the initialization is done here. Depending on the parameter task 
!! either the DNS is 
!! started, the postprocessing is called, an eigenvalue problem is solved (this feature corresponding to
!! task=3 has been migrated to a usr file by mjohn) or matrices are analyzed. The parameter is specified in
!! usr_config.f90
PROGRAM impact
  
  !USE mod_impact_core
  USE mod_vars
  USE mod_test
  USE mod_timeint
  USE usr_func

  USE HDF5
  USE MPI

  IMPLICIT NONE

  !--- Set configuration / topology --------------------------------------------------------------------------
  CALL configuration

  CALL impact_core_init

  !===========================================================================================================
  !=== Main taks =============================================================================================
  !===========================================================================================================
  !--- DNS / time integration --------------------------------------------------------------------------------
  IF (task == 1) CALL timeintegration

  !--- read and average fields -------------------------------------------------------------------------------
  IF (task == 2) CALL postprocess

  !--- solve eigenproblem (removed on 060513) ----------------------------------------------------------------
  !   IF (task == 3) CALL solve_eigenproblem(0)
  ! mjohn 060513 - see annotation above
  IF (task == 3) THEN
   if (rank .eq. 0) then
     write(*,*) 'Task 3 is no longer available. Please read information concerning update on 06 May 2013.'
   end if
  END IF

  !--- analyze matrices --------------------------------------------------------------------------------------
  IF (task == 4) CALL analyze_matrix(0)

  !--- Check Node IDs (FEM) bbecsek 2016 ---------------------------------------------------------------------
  IF (task == 5) CALL check_node_ids
  !===========================================================================================================

  CALL impact_core_finalize
  
END PROGRAM impact
