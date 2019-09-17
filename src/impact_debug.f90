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
PROGRAM impact_debug

  USE mod_impact_core
  USE mod_vars
  USE mod_test
  USE mod_timeint
  USE usr_func
  USE usr_vars

  USE HDF5
  USE MPI

  !USE mr_io_protocol, only: DistFlowMRIPadded
  USE mr_io_parallel_spacetime !, only: mr_io_read_parallel_flow_padded

  IMPLICIT NONE
  
  character(len=200) :: mri_file_path = "./bern_experimental_dataset_flow_mri.h5"
  type(DistFlowMRIPadded) :: mri_inst

  INTEGER :: gdb = 0
  do while (gdb == 0)
    call sleep(2)
  end do
  
  !--- Set configuration / topology --------------------------------------------------------------------------
  CALL configuration

  CALL impact_core_init

  ! Assign domain-padding read from config file (kalman_num_data_voxels_per_process can be used to infer the
  ! size of the data voxel grid per process)
  mri_inst%domain_padding = kalman_domain_padding
  CALL mr_io_read_parallel_flow_padded(MPI_COMM_WORLD, MPI_INFO_NULL, (/NB1,NB2,NB3/), &
       mri_file_path, mri_inst) 
  CALL h5open_f(herror) ! Required as hpc-predict-io closes HDF5 environment with h5close_f
 
  ! The indices in distributed arrays of mri_inst are aligned with the pressure grid x1p/x2p/x3p.
  ! That is for 3-dimensional MRI voxel index i = (i1,i2,i3) e.g. used to access 
  !
  !              mri_inst%mri%intensity%array(i1,i2,i3) 
  !
  ! and 3-dimensional refinement sr = (sr1,sr2,sr3) of the pressure grid relative to the MRI voxel grid
  ! all the pressure grid points with coordinates in 
  !
  !              x1p(i1*sr1+1:(i1+1)*sr1+1)
  !              x2p(i2*sr2+1:(i2+1)*sr2+1)
  !              x3p(i3*sr3+1:(i3+1)*sr3+1)
  !
  ! are within that MRI voxel (the pressure grid points with the first/last coordinate from this set always
  ! lie within the neighboring voxel as well).
 
  write(0,*) "### MRI data ('local' refers to local pressure grid) ##" 
  write(0,*) "## Number of data voxels "
  write(0,*) "   (incl. both with  "
  write(0,*) "   and without data) "
  write(0,*) "   per process:  ",kalman_num_data_voxels_per_process
  write(0,*) "## cart MPI rank:",iB(1:3,1)-(/1,1,1/)
  write(0,*) "## MRI: intensity field ##"
  write(0,*) "   local  shape: ",shape(mri_inst%mri%intensity%array)
  write(0,*) "   local lbound: ",lbound(mri_inst%mri%intensity%array)
  write(0,*) "   local ubound: ",ubound(mri_inst%mri%intensity%array)
  write(0,*) "   MRI offset:   ",mri_inst%mri%intensity%time_offset,mri_inst%mri%intensity%offset

  write(0,*) "## MRI: velocity mean field ##"
  write(0,*) "   local  shape: ",shape(mri_inst%mri%velocity_mean%array)
  write(0,*) "   local lbound: ",lbound(mri_inst%mri%velocity_mean%array)
  write(0,*) "   local ubound: ",ubound(mri_inst%mri%velocity_mean%array)
  write(0,*) "   MRI offset:   ",0,mri_inst%mri%velocity_mean%time_offset,mri_inst%mri%velocity_mean%offset

  write(0,*) "## MRI: velocity covariance field ##"
  write(0,*) "   local  shape: ",shape(mri_inst%mri%velocity_cov%array)
  write(0,*) "   local lbound: ",lbound(mri_inst%mri%velocity_cov%array)
  write(0,*) "   local ubound: ",ubound(mri_inst%mri%velocity_cov%array)
  write(0,*) "   MRI offset:   ",0,0,mri_inst%mri%velocity_cov%time_offset,mri_inst%mri%velocity_cov%offset
  
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

END PROGRAM impact_debug
