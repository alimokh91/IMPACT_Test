!       10        20        30        40        50        60        70        80        90       100       110
!        |         |         |         |         |         |         |         |         |         |         |
!2345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
!=============================================================================================================
  !===========================================================================================================
     !========================================================================================================
        !=====================================================================================================
           !==================================================================================================


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
  

  USE mod_dims
  USE mod_vars
  USE mod_setup
  USE mod_coeffs
  USE mod_lib
  USE mod_inout
  USE mod_test
  USE mod_geometry
  USE mod_timeint
  USE HDF5
  USE usr_func
  USE usr_vars
  
  USE mr_protocol, only: DistHPCPredictMRIPadded
  USE mr_io_parallel_spacetime, only: mr_io_read_parallel_hpcpredict_padded
    
  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
  character(len=200) :: mri_file_path = "./bern_experimental_dataset_hpc_predict_mri.h5"
  type(DistHPCPredictMRIPadded) :: mri_inst

  INTEGER :: gdb = 0
  do while (gdb == 0)
    call sleep(2)
  end do
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
  
  !--- Set configuration / topology --------------------------------------------------------------------------
  CALL configuration
  
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

  ! Assign domain-padding read from config file (kalman_num_data_voxels_per_process can be used to infer the
  ! size of the data voxel grid per process)
  mri_inst%domain_padding = kalman_domain_padding
  CALL mr_io_read_parallel_hpcpredict_padded(MPI_COMM_WORLD, MPI_INFO_NULL, (/NB1,NB2,NB3/), &
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
  
  
END PROGRAM impact
