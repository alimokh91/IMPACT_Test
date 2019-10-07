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
  
  type(DistFlowMRIPadded) :: mri_inst
  type(DistSpaceTimeMRI) :: mri_dest

  INTEGER :: i, ix, iy, iz, it, iv ! Helper indices for writing results

  INTEGER :: gdb = 0
  do while (gdb == 0)
    call sleep(2)
  end do
  
  !--- Set configuration / topology --------------------------------------------------------------------------
  CALL configuration

  CALL impact_core_init

  write(0,*) "Finished initialization. Reading MRI from ",kalman_mri_input_file_path," and writing to ",kalman_mri_output_file_path,"..."

  ! Assign domain-padding read from config file (kalman_num_data_voxels_per_process can be used to infer the
  ! size of the data voxel grid per process)
  mri_inst%domain_padding = kalman_domain_padding
  CALL mr_io_read_parallel_flow_padded(MPI_COMM_WORLD, MPI_INFO_NULL, (/NB1,NB2,NB3/), &
       kalman_mri_input_file_path, mri_inst)
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

  write(0,*) "## Number of refinements "
  write(0,*) "   time:         ",kalman_num_time_refinements
  write(0,*) "   space:        ",kalman_num_spatial_refinements

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

  ! From hpc-predict-io test:

  ! Allocate amount of data to be written
  ! TODO: fix time dimensions
  allocate(mri_dest%t_coordinates(size(mri_inst%mri%t_coordinates)*kalman_num_time_refinements))
  allocate(mri_dest%x_coordinates(size(mri_inst%mri%x_coordinates)*kalman_num_spatial_refinements(1)))
  allocate(mri_dest%y_coordinates(size(mri_inst%mri%y_coordinates)*kalman_num_spatial_refinements(2)))
  allocate(mri_dest%z_coordinates(size(mri_inst%mri%z_coordinates)*kalman_num_spatial_refinements(3)))
  allocate(mri_dest%voxel_feature%array(  lbound(mri_inst%mri%velocity_mean%array,1):ubound(mri_inst%mri%velocity_mean%array,1), &
                                         (lbound(mri_inst%mri%velocity_mean%array,2)-1)*kalman_num_time_refinements+1:ubound(mri_inst%mri%velocity_mean%array,2)*kalman_num_time_refinements, &
                                         (lbound(mri_inst%mri%velocity_mean%array,3)-1)*kalman_num_spatial_refinements(1)+1:ubound(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1), &
                                         (lbound(mri_inst%mri%velocity_mean%array,4)-1)*kalman_num_spatial_refinements(2)+1:ubound(mri_inst%mri%velocity_mean%array,4)*kalman_num_spatial_refinements(2), &
                                         (lbound(mri_inst%mri%velocity_mean%array,5)-1)*kalman_num_spatial_refinements(3)+1:ubound(mri_inst%mri%velocity_mean%array,5)*kalman_num_spatial_refinements(3) ))

  ! TODO: Fix all the below values
  ! time
  do i=0,size(mri_inst%mri%t_coordinates)-1  ! TODO: Fill with exact time (stepping) values
      mri_dest%t_coordinates(i*kalman_num_time_refinements+1:(i+1)*kalman_num_time_refinements) = mri_inst%mri%t_coordinates(i+1)
  end do
  mri_dest%t_dim = mri_inst%mri%t_dim*kalman_num_time_refinements

  ! geometry
  do i=0,size(mri_inst%mri%x_coordinates)-1 ! TODO: Fill with interpolation of y1p
      mri_dest%x_coordinates(i*kalman_num_spatial_refinements(1)+1:(i+1)*kalman_num_spatial_refinements(1)) = mri_inst%mri%x_coordinates(i+1)
  end do
  mri_dest%x_dim = mri_inst%mri%x_dim*kalman_num_spatial_refinements(1)

  do i=0,size(mri_inst%mri%y_coordinates)-1 ! TODO: Fill with interpolation of y2p
      mri_dest%y_coordinates(i*kalman_num_spatial_refinements(2)+1:(i+1)*kalman_num_spatial_refinements(2)) = mri_inst%mri%y_coordinates(i+1)
  end do
  mri_dest%y_dim = mri_inst%mri%y_dim*kalman_num_spatial_refinements(2)

  do i=0,size(mri_inst%mri%z_coordinates)-1 ! TODO: Fill with interpolation of y3p
      mri_dest%z_coordinates(i*kalman_num_spatial_refinements(3)+1:(i+1)*kalman_num_spatial_refinements(3)) = mri_inst%mri%z_coordinates(i+1)
  end do
  mri_dest%z_dim = mri_inst%mri%z_dim*kalman_num_spatial_refinements(3)

  ! TODO: Fill with interpolation of phi/work arrays
  ! voxel_feature
  do iz=(lbound(mri_inst%mri%velocity_mean%array,5)-1),(ubound(mri_inst%mri%velocity_mean%array,5)-1)
      do iy=(lbound(mri_inst%mri%velocity_mean%array,4)-1),(ubound(mri_inst%mri%velocity_mean%array,4)-1)
          do ix=(lbound(mri_inst%mri%velocity_mean%array,3)-1),(ubound(mri_inst%mri%velocity_mean%array,3)-1)
              do it=(lbound(mri_inst%mri%velocity_mean%array,2)-1),(ubound(mri_inst%mri%velocity_mean%array,2)-1)
                  do iv=1,3
                      mri_dest%voxel_feature%array(iv, &
                          it*kalman_num_time_refinements+1:(it+1)*kalman_num_time_refinements, &
                          ix*kalman_num_spatial_refinements(1)+1:(ix+1)*kalman_num_spatial_refinements(1), &
                          iy*kalman_num_spatial_refinements(2)+1:(iy+1)*kalman_num_spatial_refinements(2), &
                          iz*kalman_num_spatial_refinements(3)+1:(iz+1)*kalman_num_spatial_refinements(3)) &
                          = mri_inst%mri%velocity_mean%array(iv,it+1,ix+1,iy+1,iz+1)
                  end do
              end do
          end do
      end do
  end do

  ! Local hyperslab dimensions
  !TODO: Fix time dimensions, spatial dimensions can be left as is
  mri_dest%voxel_feature%time_offset = mri_inst%mri%velocity_mean%time_offset*kalman_num_time_refinements
  mri_dest%voxel_feature%time_dim = mri_inst%mri%velocity_mean%time_dim*kalman_num_time_refinements
  mri_dest%voxel_feature%offset = (/ mri_inst%mri%velocity_mean%offset(1)*kalman_num_spatial_refinements(1), &
                                     mri_inst%mri%velocity_mean%offset(2)*kalman_num_spatial_refinements(2), &
                                     mri_inst%mri%velocity_mean%offset(3)*kalman_num_spatial_refinements(3) /)
  mri_dest%voxel_feature%dims = (/ mri_inst%mri%velocity_mean%dims(1)*kalman_num_spatial_refinements(1), &
                                   mri_inst%mri%velocity_mean%dims(2)*kalman_num_spatial_refinements(2), &
                                   mri_inst%mri%velocity_mean%dims(3)*kalman_num_spatial_refinements(3) /)

  ! Write MRI-assimilated fluid simulation to output file
  call mr_io_write_parallel_spacetime(MPI_COMM_WORLD, MPI_INFO_NULL, kalman_mri_output_file_path, mri_dest)

  CALL impact_core_finalize

END PROGRAM impact_debug
