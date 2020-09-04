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

  USE mr_io_protocol!, only: DistSegmentedFlowMRIPadded
  USE mr_io_parallel_spacetime !, only: mr_io_read_parallel_flow_padded

  IMPLICIT NONE
  
  TYPE(DistSegmentedFlowMRI),   pointer :: mri
  CHARACTER(LEN=300) :: write_file
  INTEGER :: i,j,k

  INTEGER :: gdb = 1
  do while (gdb == 0)
    call sleep(2)
  end do
  
  !--- Set configuration / topology --------------------------------------------------------------------------
  CALL configuration

  CALL impact_core_init

  write(0,*) "Finished initialization. Reading MRI from ",kalman_mri_input_file_path," and writing to ",kalman_mri_output_file_path,"..."

  ! Assign domain-padding read from config file (kalman_num_data_voxels_per_process can be used to infer the
  ! size of the data voxel grid per process)
  nullify(mri_inst); allocate(mri_inst)
  mri_inst%domain_padding = kalman_domain_padding
  CALL mr_io_read_parallel_segmentedflow_padded(MPI_COMM_WORLD, MPI_INFO_NULL, (/NB1,NB2,NB3/), &
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

  write(0,*) "## cart MPI rank:",iB(1:3,1)-(/1,1,1/)
  write(0,*) "## MRI: segmentation_prob field ##"
  write(0,*) "   local  shape: ",shape(mri_inst%mri%segmentation_prob%array)
  write(0,*) "   local lbound: ",lbound(mri_inst%mri%segmentation_prob%array)
  write(0,*) "   local ubound: ",ubound(mri_inst%mri%segmentation_prob%array)
  write(0,*) "   MRI offset:   ",mri_inst%mri%intensity%time_offset,mri_inst%mri%segmentation_prob%offset


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

! Main part of IMPACT currently commented out for demonstrating MRI I/O

!  !===========================================================================================================
!  !=== Main taks =============================================================================================
!  !===========================================================================================================
!  !--- DNS / time integration --------------------------------------------------------------------------------
!  IF (task == 1) CALL timeintegration
!
!  !--- read and average fields -------------------------------------------------------------------------------
!  IF (task == 2) CALL postprocess
!
!  !--- solve eigenproblem (removed on 060513) ----------------------------------------------------------------
!  !   IF (task == 3) CALL solve_eigenproblem(0)
!  ! mjohn 060513 - see annotation above
!  IF (task == 3) THEN
!   if (rank .eq. 0) then
!     write(*,*) 'Task 3 is no longer available. Please read information concerning update on 06 May 2013.'
!   end if
!  END IF
!
!  !--- analyze matrices --------------------------------------------------------------------------------------
!  IF (task == 4) CALL analyze_matrix(0)
!
!  !--- Check Node IDs (FEM) bbecsek 2016 ---------------------------------------------------------------------
!  IF (task == 5) CALL check_node_ids
!  !===========================================================================================================


!  ! From hpc-predict-io test
  nullify(mri_flow)
  allocate(mri_flow)
  mri_flow%domain_padding = mri_inst%domain_padding
  mri_flow%domain_padding%lhs = mri_flow%domain_padding%lhs*kalman_num_spatial_refinements
  mri_flow%domain_padding%rhs = mri_flow%domain_padding%rhs*kalman_num_spatial_refinements

  !time informations
  mri_flow%mri%t_dim = mri_inst%mri%t_dim*kalman_num_time_refinements
  ALLOCATE(dtime_phases(1:intervals)); dtime_phases(1:intervals) = kalman_mri_input_attr_t_heart_cycle_period/intervals
  allocate(mri_flow%mri%t_coordinates(1:size(mri_inst%mri%t_coordinates)*kalman_num_time_refinements))
  mri_flow%mri%t_coordinates(1) = mri_inst%mri%t_coordinates(1)
  do i = 1, size(mri_inst%mri%t_coordinates)-1
     mri_flow%mri%t_coordinates(i+1) = mri_flow%mri%t_coordinates(i)+dtime_phases(i)
  end do
  mri_flow%mri%intensity%time_offset = mri_inst%mri%intensity%time_offset*kalman_num_time_refinements
  mri_flow%mri%intensity%time_dim    = mri_inst%mri%intensity%time_dim*kalman_num_time_refinements
  mri_flow%mri%segmentation_prob%time_offset = mri_inst%mri%segmentation_prob%time_offset*kalman_num_time_refinements
  mri_flow%mri%segmentation_prob%time_dim    = mri_inst%mri%segmentation_prob%time_dim*kalman_num_time_refinements
  mri_flow%mri%velocity_mean%time_offset = mri_inst%mri%velocity_mean%time_offset*kalman_num_time_refinements
  mri_flow%mri%velocity_mean%time_dim    = mri_inst%mri%velocity_mean%time_dim*kalman_num_time_refinements
  mri_flow%mri%velocity_cov%time_offset = mri_inst%mri%velocity_cov%time_offset*kalman_num_time_refinements
  mri_flow%mri%velocity_cov%time_dim    = mri_inst%mri%velocity_cov%time_dim*kalman_num_time_refinements

  !flow features
  allocate(mri_flow%mri%intensity%array( &
           (lbound(mri_inst%mri%intensity%array,1)-1)*kalman_num_time_refinements+1      :ubound(mri_inst%mri%intensity%array,1)*kalman_num_time_refinements, &
           (lbound(mri_inst%mri%intensity%array,2)-1)*kalman_num_spatial_refinements(1)+1:ubound(mri_inst%mri%intensity%array,2)*kalman_num_spatial_refinements(1), &
           (lbound(mri_inst%mri%intensity%array,3)-1)*kalman_num_spatial_refinements(2)+1:ubound(mri_inst%mri%intensity%array,3)*kalman_num_spatial_refinements(2), &
           (lbound(mri_inst%mri%intensity%array,4)-1)*kalman_num_spatial_refinements(3)+1:ubound(mri_inst%mri%intensity%array,4)*kalman_num_spatial_refinements(3) ))
  mri_flow%mri%intensity%array = 0.

  allocate(mri_flow%mri%segmentation_prob%array( &
           (lbound(mri_inst%mri%segmentation_prob%array,1)-1)*kalman_num_time_refinements+1      :ubound(mri_inst%mri%segmentation_prob%array,1)*kalman_num_time_refinements, &
           (lbound(mri_inst%mri%segmentation_prob%array,2)-1)*kalman_num_spatial_refinements(1)+1:ubound(mri_inst%mri%segmentation_prob%array,2)*kalman_num_spatial_refinements(1), &
           (lbound(mri_inst%mri%segmentation_prob%array,3)-1)*kalman_num_spatial_refinements(2)+1:ubound(mri_inst%mri%segmentation_prob%array,3)*kalman_num_spatial_refinements(2), &
           (lbound(mri_inst%mri%segmentation_prob%array,4)-1)*kalman_num_spatial_refinements(3)+1:ubound(mri_inst%mri%segmentation_prob%array,4)*kalman_num_spatial_refinements(3) ))
  mri_flow%mri%segmentation_prob%array = 0.

  allocate(mri_flow%mri%velocity_mean%array( &
            lbound(mri_inst%mri%velocity_mean%array,1)                                       :ubound(mri_inst%mri%velocity_mean%array,1), &
           (lbound(mri_inst%mri%velocity_mean%array,2)-1)*kalman_num_time_refinements+1      :ubound(mri_inst%mri%velocity_mean%array,2)*kalman_num_time_refinements, &
           (lbound(mri_inst%mri%velocity_mean%array,3)-1)*kalman_num_spatial_refinements(1)+1:ubound(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1), &
           (lbound(mri_inst%mri%velocity_mean%array,4)-1)*kalman_num_spatial_refinements(2)+1:ubound(mri_inst%mri%velocity_mean%array,4)*kalman_num_spatial_refinements(2), &
           (lbound(mri_inst%mri%velocity_mean%array,5)-1)*kalman_num_spatial_refinements(3)+1:ubound(mri_inst%mri%velocity_mean%array,5)*kalman_num_spatial_refinements(3) ))
  mri_flow%mri%velocity_mean%array = 0.

  allocate(mri_flow%mri%velocity_cov%array( &
            lbound(mri_inst%mri%velocity_cov%array,1)                                       :ubound(mri_inst%mri%velocity_cov%array,1), &
            lbound(mri_inst%mri%velocity_cov%array,2)                                       :ubound(mri_inst%mri%velocity_cov%array,2), &
           (lbound(mri_inst%mri%velocity_cov%array,3)-1)*kalman_num_time_refinements+1      :ubound(mri_inst%mri%velocity_cov%array,3)*kalman_num_time_refinements, &
           (lbound(mri_inst%mri%velocity_cov%array,4)-1)*kalman_num_spatial_refinements(1)+1:ubound(mri_inst%mri%velocity_cov%array,4)*kalman_num_spatial_refinements(1), &
           (lbound(mri_inst%mri%velocity_cov%array,5)-1)*kalman_num_spatial_refinements(2)+1:ubound(mri_inst%mri%velocity_cov%array,5)*kalman_num_spatial_refinements(2), &
           (lbound(mri_inst%mri%velocity_cov%array,6)-1)*kalman_num_spatial_refinements(3)+1:ubound(mri_inst%mri%velocity_cov%array,6)*kalman_num_spatial_refinements(3) ))
  mri_flow%mri%velocity_cov%array = 0.

  mri_flow%mri%intensity%offset = (/ mri_inst%mri%intensity%offset(1)*kalman_num_spatial_refinements(1), &
                                     mri_inst%mri%intensity%offset(2)*kalman_num_spatial_refinements(2), &
                                     mri_inst%mri%intensity%offset(3)*kalman_num_spatial_refinements(3) /)
  mri_flow%mri%intensity%dims =   (/ mri_inst%mri%intensity%dims(1)*kalman_num_spatial_refinements(1), &
                                     mri_inst%mri%intensity%dims(2)*kalman_num_spatial_refinements(2), &
                                     mri_inst%mri%intensity%dims(3)*kalman_num_spatial_refinements(3) /)

  mri_flow%mri%segmentation_prob%offset = (/ mri_inst%mri%segmentation_prob%offset(1)*kalman_num_spatial_refinements(1), &
                                             mri_inst%mri%segmentation_prob%offset(2)*kalman_num_spatial_refinements(2), &
                                             mri_inst%mri%segmentation_prob%offset(3)*kalman_num_spatial_refinements(3) /)
  mri_flow%mri%segmentation_prob%dims =   (/ mri_inst%mri%segmentation_prob%dims(1)*kalman_num_spatial_refinements(1), &
                                             mri_inst%mri%segmentation_prob%dims(2)*kalman_num_spatial_refinements(2), &
                                             mri_inst%mri%segmentation_prob%dims(3)*kalman_num_spatial_refinements(3) /)

  mri_flow%mri%velocity_mean%offset = (/ mri_inst%mri%velocity_mean%offset(1)*kalman_num_spatial_refinements(1), &
                                         mri_inst%mri%velocity_mean%offset(2)*kalman_num_spatial_refinements(2), &
                                         mri_inst%mri%velocity_mean%offset(3)*kalman_num_spatial_refinements(3) /)
  mri_flow%mri%velocity_mean%dims =   (/ mri_inst%mri%velocity_mean%dims(1)*kalman_num_spatial_refinements(1), &
                                         mri_inst%mri%velocity_mean%dims(2)*kalman_num_spatial_refinements(2), &
                                         mri_inst%mri%velocity_mean%dims(3)*kalman_num_spatial_refinements(3) /)

  mri_flow%mri%velocity_cov%offset = (/ mri_inst%mri%velocity_cov%offset(1)*kalman_num_spatial_refinements(1), &
                                        mri_inst%mri%velocity_cov%offset(2)*kalman_num_spatial_refinements(2), &
                                        mri_inst%mri%velocity_cov%offset(3)*kalman_num_spatial_refinements(3) /)
  mri_flow%mri%velocity_cov%dims =   (/ mri_inst%mri%velocity_cov%dims(1)*kalman_num_spatial_refinements(1), &
                                        mri_inst%mri%velocity_cov%dims(2)*kalman_num_spatial_refinements(2), &
                                        mri_inst%mri%velocity_cov%dims(3)*kalman_num_spatial_refinements(3) /)  

  !spatial informations
  mri_flow%mri%x_dim = mri_inst%mri%x_dim*kalman_num_spatial_refinements(1)
  mri_flow%mri%y_dim = mri_inst%mri%y_dim*kalman_num_spatial_refinements(2)
  mri_flow%mri%z_dim = mri_inst%mri%z_dim*kalman_num_spatial_refinements(3)

  allocate(mri_flow%mri%x_coordinates(1:size(mri_inst%mri%x_coordinates)*kalman_num_spatial_refinements(1)))
  allocate(mri_flow%mri%y_coordinates(1:size(mri_inst%mri%y_coordinates)*kalman_num_spatial_refinements(2)))
  allocate(mri_flow%mri%z_coordinates(1:size(mri_inst%mri%z_coordinates)*kalman_num_spatial_refinements(3)))

  mri_flow%mri%x_coordinates = 0.5*y1p( mri_flow%domain_padding%lhs(1)+1:                  &
                                       (mri_flow%domain_padding%lhs(1)+size(mri_flow%mri%x_coordinates))+0 )
  mri_flow%mri%x_coordinates = mri_flow%mri%x_coordinates + &
                               0.5*y1p( mri_flow%domain_padding%lhs(1)+2:                  &
                                       (mri_flow%domain_padding%lhs(1)+size(mri_flow%mri%x_coordinates))+1 )

  mri_flow%mri%y_coordinates = 0.5*y2p( mri_flow%domain_padding%lhs(2)+1:                  &
                                       (mri_flow%domain_padding%lhs(2)+size(mri_flow%mri%y_coordinates))+0 )
  mri_flow%mri%y_coordinates = mri_flow%mri%y_coordinates + &
                               0.5*y2p( mri_flow%domain_padding%lhs(2)+2:                  &
                                       (mri_flow%domain_padding%lhs(2)+size(mri_flow%mri%y_coordinates))+1 )

  mri_flow%mri%z_coordinates = 0.5*y3p( mri_flow%domain_padding%lhs(3)+1:                  &
                                       (mri_flow%domain_padding%lhs(3)+size(mri_flow%mri%z_coordinates))+0 )
  mri_flow%mri%z_coordinates = mri_flow%mri%z_coordinates + &
                               0.5*y3p( mri_flow%domain_padding%lhs(3)+2:                  &
                                       (mri_flow%domain_padding%lhs(3)+size(mri_flow%mri%z_coordinates))+1 )

    !=========================================================================================================
    !=== write mri_flow at the end ===========================================================================
    !=========================================================================================================
    write_file = kalman_mri_output_file_path
    CALL mr_io_write_parallel_segmentedflow(MPI_COMM_WORLD, MPI_INFO_NULL, write_file, mri_flow%mri)
    CALL h5open_f(herror)
    !=========================================================================================================

  CALL impact_core_finalize

END PROGRAM impact_debug
