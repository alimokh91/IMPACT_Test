program mr_io_test_impact_mri

  USE mod_impact_core
  USE mod_vars
  USE usr_vars
  USE mr_io_parallel_spacetime
  USE hdf5
  USE mpi

  implicit none

  type(DistFlowMRIPadded) :: mri_inst2
  type(DistSpaceTimeMRI) :: mri_dest2

  INTEGER :: i, ix, iy, iz, it, iv ! Helper indices for writing results

  integer :: lbound_v, lbound_t, lbound_x, lbound_y, lbound_z
  integer :: ubound_v, ubound_t, ubound_x, ubound_y, ubound_z

  real*8, dimension(:,:,:,:), allocatable :: pressure_voxel_velocity
!  real*8, dimension(:,:,:,:), allocatable :: pressure_voxel_mri_cell_pattern


!  ! For interactive testing
!  INTEGER :: gdb = 0
!  do while (gdb == 0)
!      call sleep(2)
!  end do

  !--- Set configuration / topology --------------------------------------------------------------------------
  CALL configuration

  CALL impact_core_init

  write(0,*) "Finished initialization. Reading MRI from ",kalman_mri_input_file_path," and writing to ",kalman_mri_output_file_path,"..."

  ! Assign domain-padding read from config file (kalman_num_data_voxels_per_process can be used to infer the
  ! size of the data voxel grid per process)
  mri_inst2%domain_padding = kalman_domain_padding
  CALL mr_io_read_parallel_flow_padded(MPI_COMM_WORLD, MPI_INFO_NULL, (/NB1,NB2,NB3/), &
       kalman_mri_input_file_path, mri_inst2)
  CALL h5open_f(herror) ! Required as hpc-predict-io closes HDF5 environment with h5close_f

  ! The indices in distributed arrays of mri_inst2 are aligned with the pressure grid x1p/x2p/x3p.
  ! That is for 3-dimensional MRI voxel index i = (i1,i2,i3) e.g. used to access
  !
  !              mri_inst2%mri%intensity%array(i1,i2,i3)
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
  write(0,*) "   local  shape: ",shape(mri_inst2%mri%intensity%array)
  write(0,*) "   local lbound: ",lbound(mri_inst2%mri%intensity%array)
  write(0,*) "   local ubound: ",ubound(mri_inst2%mri%intensity%array)
  write(0,*) "   MRI offset:   ",mri_inst2%mri%intensity%time_offset,mri_inst2%mri%intensity%offset

  write(0,*) "## MRI: velocity mean field ##"
  write(0,*) "   local  shape: ",shape(mri_inst2%mri%velocity_mean%array)
  write(0,*) "   local lbound: ",lbound(mri_inst2%mri%velocity_mean%array)
  write(0,*) "   local ubound: ",ubound(mri_inst2%mri%velocity_mean%array)
  write(0,*) "   MRI offset:   ",0,mri_inst2%mri%velocity_mean%time_offset,mri_inst2%mri%velocity_mean%offset

  write(0,*) "## MRI: velocity covariance field ##"
  write(0,*) "   local  shape: ",shape(mri_inst2%mri%velocity_cov%array)
  write(0,*) "   local lbound: ",lbound(mri_inst2%mri%velocity_cov%array)
  write(0,*) "   local ubound: ",ubound(mri_inst2%mri%velocity_cov%array)
  write(0,*) "   MRI offset:   ",0,0,mri_inst2%mri%velocity_cov%time_offset,mri_inst2%mri%velocity_cov%offset

  write(0,*) "Performing fluid simulation..."
  call sleep(1)
  write(0,*) "Done. Write results to MRI file in ",kalman_mri_output_file_path,"."
  call sleep(1)


  ! Allocate amount of data to be written - TODO: Fix time dimensions
  allocate(mri_dest2%t_coordinates(size(mri_inst2%mri%t_coordinates)*kalman_num_time_refinements))
  allocate(mri_dest2%x_coordinates(size(mri_inst2%mri%x_coordinates)*kalman_num_spatial_refinements(1)))
  allocate(mri_dest2%y_coordinates(size(mri_inst2%mri%y_coordinates)*kalman_num_spatial_refinements(2)))
  allocate(mri_dest2%z_coordinates(size(mri_inst2%mri%z_coordinates)*kalman_num_spatial_refinements(3)))
  allocate(mri_dest2%vector_feature%array(  lbound(mri_inst2%mri%velocity_mean%array,1):ubound(mri_inst2%mri%velocity_mean%array,1), &
                                         (lbound(mri_inst2%mri%velocity_mean%array,2)-1)*kalman_num_time_refinements+1:ubound(mri_inst2%mri%velocity_mean%array,2)*kalman_num_time_refinements, &
                                         (lbound(mri_inst2%mri%velocity_mean%array,3)-1)*kalman_num_spatial_refinements(1)+1:ubound(mri_inst2%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1), &
                                         (lbound(mri_inst2%mri%velocity_mean%array,4)-1)*kalman_num_spatial_refinements(2)+1:ubound(mri_inst2%mri%velocity_mean%array,4)*kalman_num_spatial_refinements(2), &
                                         (lbound(mri_inst2%mri%velocity_mean%array,5)-1)*kalman_num_spatial_refinements(3)+1:ubound(mri_inst2%mri%velocity_mean%array,5)*kalman_num_spatial_refinements(3) ))

  ! Write coordinates - TODO: Fix time dimensions values
  ! time
  do i=0,size(mri_inst2%mri%t_coordinates)-1  ! TODO: Fill with exact time (stepping) values
      mri_dest2%t_coordinates(i*kalman_num_time_refinements+1:(i+1)*kalman_num_time_refinements) = mri_inst2%mri%t_coordinates(i+1)
  end do
  mri_dest2%t_dim = mri_inst2%mri%t_dim*kalman_num_time_refinements

  ! geometry
  mri_dest2%x_coordinates = 0.5*y1p(kalman_domain_padding%lhs(1)*kalman_num_spatial_refinements(1)+1:(kalman_domain_padding%lhs(1)+size(mri_inst2%mri%x_coordinates))*kalman_num_spatial_refinements(1)+0) + &
                           0.5*y1p(kalman_domain_padding%lhs(1)*kalman_num_spatial_refinements(1)+2:(kalman_domain_padding%lhs(1)+size(mri_inst2%mri%x_coordinates))*kalman_num_spatial_refinements(1)+1)
  mri_dest2%x_dim = mri_inst2%mri%x_dim*kalman_num_spatial_refinements(1)
  mri_dest2%y_coordinates = 0.5*y2p(kalman_domain_padding%lhs(2)*kalman_num_spatial_refinements(2)+1:(kalman_domain_padding%lhs(2)+size(mri_inst2%mri%y_coordinates))*kalman_num_spatial_refinements(2)+0) + &
                           0.5*y2p(kalman_domain_padding%lhs(2)*kalman_num_spatial_refinements(2)+2:(kalman_domain_padding%lhs(2)+size(mri_inst2%mri%y_coordinates))*kalman_num_spatial_refinements(2)+1)
  mri_dest2%y_dim = mri_inst2%mri%y_dim*kalman_num_spatial_refinements(2)
  mri_dest2%z_coordinates = 0.5*y3p(kalman_domain_padding%lhs(3)*kalman_num_spatial_refinements(3)+1:(kalman_domain_padding%lhs(3)+size(mri_inst2%mri%z_coordinates))*kalman_num_spatial_refinements(3)+0) + &
                           0.5*y3p(kalman_domain_padding%lhs(3)*kalman_num_spatial_refinements(3)+2:(kalman_domain_padding%lhs(3)+size(mri_inst2%mri%z_coordinates))*kalman_num_spatial_refinements(3)+1)
  mri_dest2%z_dim = mri_inst2%mri%z_dim*kalman_num_spatial_refinements(3)

  ! Local hyperslab dimensions - TODO: Fix time dimensions
  mri_dest2%vector_feature%time_offset = mri_inst2%mri%velocity_mean%time_offset*kalman_num_time_refinements
  mri_dest2%vector_feature%time_dim = mri_inst2%mri%velocity_mean%time_dim*kalman_num_time_refinements
  mri_dest2%vector_feature%offset = (/ mri_inst2%mri%velocity_mean%offset(1)*kalman_num_spatial_refinements(1), &
                                     mri_inst2%mri%velocity_mean%offset(2)*kalman_num_spatial_refinements(2), &
                                     mri_inst2%mri%velocity_mean%offset(3)*kalman_num_spatial_refinements(3) /)
  mri_dest2%vector_feature%dims = (/ mri_inst2%mri%velocity_mean%dims(1)*kalman_num_spatial_refinements(1), &
                                   mri_inst2%mri%velocity_mean%dims(2)*kalman_num_spatial_refinements(2), &
                                   mri_inst2%mri%velocity_mean%dims(3)*kalman_num_spatial_refinements(3) /)

  ! Local hyperslab boundary indexes - TODO: Fix time dimensions
  lbound_v = lbound(mri_inst2%mri%velocity_mean%array,1)
  ubound_v = ubound(mri_inst2%mri%velocity_mean%array,1)
  lbound_t = (lbound(mri_inst2%mri%velocity_mean%array,2)-1)*kalman_num_time_refinements+1
  ubound_t = ubound(mri_inst2%mri%velocity_mean%array,2)*kalman_num_time_refinements
  lbound_x = (lbound(mri_inst2%mri%velocity_mean%array,3)-1)*kalman_num_spatial_refinements(1)+1
  ubound_x = ubound(mri_inst2%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)
  lbound_y = (lbound(mri_inst2%mri%velocity_mean%array,4)-1)*kalman_num_spatial_refinements(2)+1
  ubound_y = ubound(mri_inst2%mri%velocity_mean%array,4)*kalman_num_spatial_refinements(2)
  lbound_z = (lbound(mri_inst2%mri%velocity_mean%array,5)-1)*kalman_num_spatial_refinements(3)+1
  ubound_z = ubound(mri_inst2%mri%velocity_mean%array,5)*kalman_num_spatial_refinements(3)

  allocate(pressure_voxel_velocity(lbound_v:ubound_v, lbound_x:ubound_x, lbound_y:ubound_y, lbound_z:ubound_z))

!  allocate(pressure_voxel_mri_cell_pattern(lbound_v:ubound_v, &
!                                           kalman_num_spatial_refinements(1), &
!                                           kalman_num_spatial_refinements(2), &
!                                           kalman_num_spatial_refinements(3)))
!  ! Add row major index pattern to randomize data
!  do iz=lbound(pressure_voxel_mri_cell_pattern,3),ubound(pressure_voxel_mri_cell_pattern,3)
!      do iy=lbound(pressure_voxel_mri_cell_pattern,2),ubound(pressure_voxel_mri_cell_pattern,2)
!          do ix=lbound(pressure_voxel_mri_cell_pattern,1),ubound(pressure_voxel_mri_cell_pattern,1)
!              do iv=lbound_v,ubound_v
!                  pressure_voxel_mri_cell_pattern(iv,ix,iy,iz) = &
!                      iv + (ubound_v-lbound_v+1)*(iz + kalman_num_spatial_refinements(3)* &
!                                                 (iy + kalman_num_spatial_refinements(2)* &
!                                                  ix))
!              end do
!          end do
!      end do
!  end do

  ! Mimick time integration loop
  do it=lbound_t,ubound_t

      ! could be made time-dependent
      do iz=(lbound(mri_inst2%mri%velocity_mean%array,5)-1),(ubound(mri_inst2%mri%velocity_mean%array,5)-1)
          do iy=(lbound(mri_inst2%mri%velocity_mean%array,4)-1),(ubound(mri_inst2%mri%velocity_mean%array,4)-1)
              do ix=(lbound(mri_inst2%mri%velocity_mean%array,3)-1),(ubound(mri_inst2%mri%velocity_mean%array,3)-1)
                  do iv=1,3
                      pressure_voxel_velocity(iv, &
                          ix*kalman_num_spatial_refinements(1)+1:(ix+1)*kalman_num_spatial_refinements(1), &
                          iy*kalman_num_spatial_refinements(2)+1:(iy+1)*kalman_num_spatial_refinements(2), &
                          iz*kalman_num_spatial_refinements(3)+1:(iz+1)*kalman_num_spatial_refinements(3)) &
                          = mri_inst2%mri%velocity_mean%array(iv,(it-1)/kalman_num_time_refinements+1,ix+1,iy+1,iz+1) ! + pressure_voxel_mri_cell_pattern(iv,:,:,:)
                  end do
              end do
          end do
      end do

      mri_dest2%vector_feature%array(:, it, lbound_x:ubound_x, lbound_y:ubound_y, lbound_z:ubound_z) = pressure_voxel_velocity
  end do


  !****** To write the velocity field in IMPACT use the following commented code ******
!  ! time integration loop
!  do it=lbound_t,ubound_t
!
!      vel = ... ! compute velocity field during time integration
!      CALL interpolate_vel(.TRUE.)  ! transfer velocity to pressure grid (requires USE mod_diff): vel(:,:,:,i) --> worki(:,:,:)
!
!      pressure_voxel_velocity(1,:,:,:) &
!          = 0.125*(work1(lbound_x+1:ubound_x+1, lbound_y+1:ubound_y+1, lbound_z+1:ubound_z+1) + &
!                   work1(lbound_x-1:ubound_x+0, lbound_y+1:ubound_y+1, lbound_z+1:ubound_z+1) + &
!                   work1(lbound_x+1:ubound_x+1, lbound_y+0:ubound_y+0, lbound_z+1:ubound_z+1) + &
!                   work1(lbound_x+1:ubound_x+1, lbound_y+1:ubound_y+1, lbound_z+0:ubound_z+0) + &
!                   work1(lbound_x+1:ubound_x+1, lbound_y+0:ubound_y+0, lbound_z+0:ubound_z+0) + &
!                   work1(lbound_x+0:ubound_x+0, lbound_y+1:ubound_y+1, lbound_z+0:ubound_z+0) + &
!                   work1(lbound_x+0:ubound_x+0, lbound_y+0:ubound_y+0, lbound_z+1:ubound_z+1) + &
!                   work1(lbound_x+0:ubound_x+0, lbound_y+0:ubound_y+0, lbound_z+0:ubound_z+0) )
!
!      pressure_voxel_velocity(2,:,:,:) &
!          = 0.125*(work2(lbound_x+1:ubound_x+1, lbound_y+1:ubound_y+1, lbound_z+1:ubound_z+1) + &
!                   work2(lbound_x+0:ubound_x+0, lbound_y+1:ubound_y+1, lbound_z+1:ubound_z+1) + &
!                   work2(lbound_x+1:ubound_x+1, lbound_y+0:ubound_y+0, lbound_z+1:ubound_z+1) + &
!                   work2(lbound_x+1:ubound_x+1, lbound_y+1:ubound_y+1, lbound_z+0:ubound_z+0) + &
!                   work2(lbound_x+1:ubound_x+1, lbound_y+0:ubound_y+0, lbound_z+0:ubound_z+0) + &
!                   work2(lbound_x+0:ubound_x+0, lbound_y+1:ubound_y+1, lbound_z+0:ubound_z+0) + &
!                   work2(lbound_x+0:ubound_x+0, lbound_y+0:ubound_y+0, lbound_z+1:ubound_z+1) + &
!                   work2(lbound_x+0:ubound_x+0, lbound_y+0:ubound_y+0, lbound_z+0:ubound_z+0) )
!
!      pressure_voxel_velocity(3,:,:,:) &
!          = 0.125*(work3(lbound_x+1:ubound_x+1, lbound_y+1:ubound_y+1, lbound_z+1:ubound_z+1) + &
!                   work3(lbound_x+0:ubound_x+0, lbound_y+1:ubound_y+1, lbound_z+1:ubound_z+1) + &
!                   work3(lbound_x+1:ubound_x+1, lbound_y+0:ubound_y+0, lbound_z+1:ubound_z+1) + &
!                   work3(lbound_x+1:ubound_x+1, lbound_y+1:ubound_y+1, lbound_z+0:ubound_z+0) + &
!                   work3(lbound_x+1:ubound_x+1, lbound_y+0:ubound_y+0, lbound_z+0:ubound_z+0) + &
!                   work3(lbound_x+0:ubound_x+0, lbound_y+1:ubound_y+1, lbound_z+0:ubound_z+0) + &
!                   work3(lbound_x+0:ubound_x+0, lbound_y+0:ubound_y+0, lbound_z+1:ubound_z+1) + &
!                   work3(lbound_x+0:ubound_x+0, lbound_y+0:ubound_y+0, lbound_z+0:ubound_z+0) )
!
!      mri_dest2%vector_feature%array(:, it, lbound_x:ubound_x, lbound_y:ubound_y, lbound_z:ubound_z) = pressure_voxel_velocity
!  end do

  ! Write MRI-assimilated fluid simulation to output file
  call mr_io_write_parallel_spacetime(MPI_COMM_WORLD, MPI_INFO_NULL, kalman_mri_output_file_path, mri_dest2)

  CALL h5open_f(herror)

  CALL impact_core_finalize

end program mr_io_test_impact_mri
