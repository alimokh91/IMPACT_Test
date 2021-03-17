!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Dario De Marinis, ARTORG CVE, University of Bern (dario.demarinis@artorg.unibe.ch)            *
!* November 2017                                                                                             *
!*************************************************************************************************************

! file stats.f90
! file containing subroutines used for writing and reading statistical data to and from files.

!!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8

  SUBROUTINE open_stats
  ! (basic subroutine)

  USE mod_dims
  USE mod_vars
  USE usr_func
  USE usr_vars
  USE mod_lib !added for writing turb_statxz_ with M1 M2 M3
  USE MPI
  USE HDF5
  USE mr_io_protocol
  USE mr_io_parallel_spacetime

  IMPLICIT NONE
 
  TYPE(DistSegmentedFlowMRI), pointer :: mri
  INTEGER :: l,n,i,j,k,ii,jj,kk,id
  INTEGER :: n_data,n_procs,strd(3),pw
  INTEGER, ALLOCATABLE :: n_data_count(:)
  REAL :: x_min,x_max,y_min,y_max,z_min,z_max
  CHARACTER(LEN=300) :: write_dir
  CHARACTER(LEN=2) ::  id_char, phs
  CHARACTER(LEN=3) ::  M1_char, M2_char, M3_char
 
  if (rank.eq.0) write(0,*) "open stats... "

  intervals = kalman_num_time_refinements*intervals

  ALLOCATE(dtime_phases(1:intervals)); dtime_phases(1:intervals) = kalman_mri_input_attr_t_heart_cycle_period/intervals !dtime_out_scal 
  if (allocated(dtime_kalm_phases)) then
     do i = 1, intervals/kalman_num_time_refinements
        j = (i-1)*kalman_num_time_refinements
        dtime_phases(j+1:j+kalman_num_time_refinements) = dtime_kalm_phases(i)/kalman_num_time_refinements
     end do

  end if
  dtime_out_scal = dtime_phases(1)

  write_stats_count = 0
 
  intervals = intervals/kalman_num_time_refinements

  !===========================================================================================================
  !=== construct pressure-grid mri data-structure ============================================================
  !===========================================================================================================
  nullify(mri_flow) 

  !===========================================================================================================
  !added for writing tke_kalman with M1 M2 M3 ================================================================
  !===========================================================================================================
  IF (rank == 0 .AND. dtime_out_kalm /= 0.) THEN
     CALL num_to_string(3,restart,restart_char)
     CALL num_to_string(3,M1,M1_char)
     CALL num_to_string(3,M2,M2_char)
     CALL num_to_string(3,M3,M3_char)
     CALL system('mkdir -p kf_result')
     OPEN(63,FILE='./kf_result/tke_window_kalm_'//restart_char//'_'//M1_char//'x'//M2_char//'x'//M3_char//'.txt',STATUS='UNKNOWN')
     do n = 1,intervals
        CALL num_to_string(2,n,phs)
        CALL system('mkdir -p kf_result/phase_'//phs )
     end do
  END IF
  !===========================================================================================================

  allocate(mri_flow)
  mri_flow%domain_padding = mri_inst%domain_padding
  mri_flow%domain_padding%lhs = mri_flow%domain_padding%lhs*kalman_num_spatial_refinements
  mri_flow%domain_padding%rhs = mri_flow%domain_padding%rhs*kalman_num_spatial_refinements

  !time informations
  mri_flow%mri%t_dim = mri_inst%mri%t_dim*kalman_num_time_refinements
  allocate(mri_flow%mri%t_coordinates(1:size(mri_inst%mri%t_coordinates)*kalman_num_time_refinements))
  mri_flow%mri%t_coordinates(1) = mri_inst%mri%t_coordinates(1)
  do i = 1, size(mri_flow%mri%t_coordinates)-1
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
  !===========================================================================================================


  END SUBROUTINE open_stats
 
  SUBROUTINE close_stats
  ! (basic subroutine)
 
  USE mod_dims
  USE mod_vars
  USE usr_func
  USE usr_vars
  USE MPI
  USE HDF5
  USE mr_io_protocol
  USE mr_io_parallel_spacetime
 
  IMPLICIT NONE
 
  INTEGER :: id
 
  IF (rank == 0 .AND. dtime_out_scal /= 0.) THEN
    CLOSE(23)
    CLOSE(33)
  END IF

  IF (rank == 0 .AND. dtime_out_kalm /= 0.) THEN
     CLOSE(63)
  END IF

  DEALLOCATE(dtime_phases)

  if (associated(mri_flow)) then
     call mr_io_deallocate_dist_segmentedflow_mri_padded(mri_flow)
  end if

  END SUBROUTINE close_stats

 
  SUBROUTINE compute_stats
  ! (basic subroutine)
 
  USE mod_dims
  USE mod_vars
  USE mod_lib
  USE mod_diff
  USE mod_inout
  USE usr_func
  USE usr_vars
  USE MPI
  USE HDF5
  USE mr_io_protocol
  USE mr_io_parallel_spacetime
 
  IMPLICIT NONE
 
  TYPE(DistSegmentedFlowMRI),   pointer :: mri
  INTEGER :: i,j,k,ii,jj,kk
  INTEGER, DIMENSION(2,3) :: bounds
  REAL    :: TKE(1:2), TKE_global(1:2)
  REAL    :: vel_voxel(1:3)
  CHARACTER(LEN=8) :: count_char
  CHARACTER(LEN=2) :: id
  CHARACTER(LEN=50) :: write_file

  intervals = intervals*kalman_num_time_refinements
  phase = mod(write_stats_count,intervals) + 1
 
  IF (rank == 0) WRITE(*,'(a,i8,a,i8,a)') 'data stats repetition', write_stats_count/intervals+1, '   for phase', phase,' ...'
  !===========================================================================================================
  !=== interpolation of velocity to pressure grid points =====================================================
  !===========================================================================================================
  ! vel(:,:,:,i) --> worki(:,:,:)
  IF (task == 1) CALL interpolate_vel(.TRUE.)
  !===========================================================================================================

  !===========================================================================================================
  !=== write pressure-grid mri data-fields ===================================================================
  !===========================================================================================================
  IF (associated(mri_flow)) then

    mri => mri_flow%mri
    bounds(1,1) = lbound(mri%velocity_mean%array,3)
    bounds(2,1) = ubound(mri%velocity_mean%array,3)
    bounds(1,2) = lbound(mri%velocity_mean%array,4)
    bounds(2,2) = ubound(mri%velocity_mean%array,4)
    bounds(1,3) = lbound(mri%velocity_mean%array,5)
    bounds(2,3) = ubound(mri%velocity_mean%array,5)

    !--------------------------------------------------------------------------------------------------------
    !--- <u_i u_j>(t_k) = <u'_i u'_j>(t_k) + <u_i>(t_k) <u_j>(t_k) ------------------------------------------
    !--------------------------------------------------------------------------------------------------------
    mri%velocity_cov%array(1,1,phase,:,:,:) = mri%velocity_cov%array(1,1,phase,:,:,:) + &
                                              mri%velocity_mean%array(1,phase,:,:,:)*mri%velocity_mean%array(1,phase,:,:,:)
    mri%velocity_cov%array(2,2,phase,:,:,:) = mri%velocity_cov%array(2,2,phase,:,:,:) + &
                                              mri%velocity_mean%array(2,phase,:,:,:)*mri%velocity_mean%array(2,phase,:,:,:)
    mri%velocity_cov%array(3,3,phase,:,:,:) = mri%velocity_cov%array(3,3,phase,:,:,:) + &
                                              mri%velocity_mean%array(3,phase,:,:,:)*mri%velocity_mean%array(3,phase,:,:,:)
    mri%velocity_cov%array(1,2,phase,:,:,:) = mri%velocity_cov%array(1,2,phase,:,:,:) + &
                                              mri%velocity_mean%array(1,phase,:,:,:)*mri%velocity_mean%array(2,phase,:,:,:)
    mri%velocity_cov%array(1,3,phase,:,:,:) = mri%velocity_cov%array(1,3,phase,:,:,:) + &
                                              mri%velocity_mean%array(1,phase,:,:,:)*mri%velocity_mean%array(3,phase,:,:,:)
    mri%velocity_cov%array(2,3,phase,:,:,:) = mri%velocity_cov%array(2,3,phase,:,:,:) + &
                                              mri%velocity_mean%array(2,phase,:,:,:)*mri%velocity_mean%array(3,phase,:,:,:)
    !--------------------------------------------------------------------------------------------------------
    !--- <u_i>(t_k) -----------------------------------------------------------------------------------------
    !--- <u_i u_j>(t_k) -------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------
    mri%velocity_mean%array (:,phase,:,:,:) = mri%velocity_mean%array (:,phase,:,:,:)*(write_stats_count/intervals)
    mri%velocity_cov%array(:,:,phase,:,:,:) = mri%velocity_cov%array(:,:,phase,:,:,:)*(write_stats_count/intervals)
    DO k = bounds(1,3), bounds(2,3)
      DO j = bounds(1,2), bounds(2,2)
        DO i = bounds(1,1), bounds(2,1)
          vel_voxel(1:3) = 0.0 
          DO kk = k,k+1
            DO jj = j,j+1
              DO ii = i,i+1
                vel_voxel(1) = vel_voxel(1) + 0.125*work1(ii,jj,kk) 
                vel_voxel(2) = vel_voxel(2) + 0.125*work2(ii,jj,kk) 
                vel_voxel(3) = vel_voxel(3) + 0.125*work3(ii,jj,kk) 
              END DO
            END DO
          END DO
          mri%velocity_mean%array(1,phase,i,j,k) = mri%velocity_mean%array(1,phase,i,j,k) + vel_voxel(1)
          mri%velocity_mean%array(2,phase,i,j,k) = mri%velocity_mean%array(2,phase,i,j,k) + vel_voxel(2)
          mri%velocity_mean%array(3,phase,i,j,k) = mri%velocity_mean%array(3,phase,i,j,k) + vel_voxel(3)

          mri%velocity_cov%array(1,1,phase,i,j,k) = mri%velocity_cov%array(1,1,phase,i,j,k) + vel_voxel(1)*vel_voxel(1)
          mri%velocity_cov%array(2,2,phase,i,j,k) = mri%velocity_cov%array(2,2,phase,i,j,k) + vel_voxel(2)*vel_voxel(2)
          mri%velocity_cov%array(3,3,phase,i,j,k) = mri%velocity_cov%array(3,3,phase,i,j,k) + vel_voxel(3)*vel_voxel(3)
          mri%velocity_cov%array(1,2,phase,i,j,k) = mri%velocity_cov%array(1,2,phase,i,j,k) + vel_voxel(1)*vel_voxel(2)
          mri%velocity_cov%array(1,3,phase,i,j,k) = mri%velocity_cov%array(1,3,phase,i,j,k) + vel_voxel(1)*vel_voxel(3)
          mri%velocity_cov%array(2,3,phase,i,j,k) = mri%velocity_cov%array(2,3,phase,i,j,k) + vel_voxel(2)*vel_voxel(3)
        END DO
      END DO
    END DO
    mri%velocity_mean%array (:,phase,:,:,:) = mri%velocity_mean%array (:,phase,:,:,:)/(write_stats_count/intervals+1)
    mri%velocity_cov%array(:,:,phase,:,:,:) = mri%velocity_cov%array(:,:,phase,:,:,:)/(write_stats_count/intervals+1)
    !-----------------------------------------------------------------------------------------------------------
    !--- <u'_i u'_j>(t_k) --------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------------------
    mri%velocity_cov%array(1,1,phase,:,:,:) = mri%velocity_cov%array(1,1,phase,:,:,:) - &
                                              mri%velocity_mean%array(1,phase,:,:,:)*mri%velocity_mean%array(1,phase,:,:,:)
    mri%velocity_cov%array(2,2,phase,:,:,:) = mri%velocity_cov%array(2,2,phase,:,:,:) - &
                                              mri%velocity_mean%array(2,phase,:,:,:)*mri%velocity_mean%array(2,phase,:,:,:)
    mri%velocity_cov%array(3,3,phase,:,:,:) = mri%velocity_cov%array(3,3,phase,:,:,:) - &
                                              mri%velocity_mean%array(3,phase,:,:,:)*mri%velocity_mean%array(3,phase,:,:,:)
    mri%velocity_cov%array(1,2,phase,:,:,:) = mri%velocity_cov%array(1,2,phase,:,:,:) - &
                                              mri%velocity_mean%array(1,phase,:,:,:)*mri%velocity_mean%array(2,phase,:,:,:)
    mri%velocity_cov%array(1,3,phase,:,:,:) = mri%velocity_cov%array(1,3,phase,:,:,:) - &
                                              mri%velocity_mean%array(1,phase,:,:,:)*mri%velocity_mean%array(3,phase,:,:,:)
    mri%velocity_cov%array(2,3,phase,:,:,:) = mri%velocity_cov%array(2,3,phase,:,:,:) - &
                                              mri%velocity_mean%array(2,phase,:,:,:)*mri%velocity_mean%array(3,phase,:,:,:)
    mri%velocity_cov%array(2,1,phase,:,:,:) = mri%velocity_cov%array(1,2,phase,:,:,:)
    mri%velocity_cov%array(3,1,phase,:,:,:) = mri%velocity_cov%array(1,3,phase,:,:,:)
    mri%velocity_cov%array(3,2,phase,:,:,:) = mri%velocity_cov%array(2,3,phase,:,:,:)
    !===========================================================================================================
    
    !========================================================================================================
    !=== integral energy for the Kalman window ==============================================================
    !========================================================================================================
    TKE = 0.
    DO k = bounds(1,3), bounds(2,3)
      DO j = bounds(1,2), bounds(2,2)
        DO i = bounds(1,1), bounds(2,1)
          !--------------------------------------------------------------------------------------------------------
          !--- <u_i> <u_i>(t_k) -----------------------------------------------------------------------------------
          !--------------------------------------------------------------------------------------------------------
          TKE(1) = TKE(1) + (mri%velocity_mean%array(1,phase,i,j,k)**2 + mri%velocity_mean%array(2,phase,i,j,k)**2 + &
                             mri%velocity_mean%array(3,phase,i,j,k)**2 )
          !--------------------------------------------------------------------------------------------------------
          !--- <u'_i u'_i>(t_k) -----------------------------------------------------------------------------------
          !--------------------------------------------------------------------------------------------------------
          TKE(2) = TKE(2) + (mri%velocity_cov%array(1,1,phase,i,j,k) + mri%velocity_cov%array(2,2,phase,i,j,k) + &
                             mri%velocity_cov%array(3,3,phase,i,j,k) )
        END DO
      END DO
    END DO
    CALL MPI_ALLREDUCE(TKE,TKE_global,2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
    TKE_global = (mri%x_coordinates(2)-mri%x_coordinates(1)) * (mri%y_coordinates(2)-mri%y_coordinates(1)) * &
                 (mri%z_coordinates(2)-mri%z_coordinates(1)) * TKE_global/2.

    IF (rank == 0) THEN
      WRITE(63,'(3E25.17)') time, TKE_global(1:2)
      CALL flush(63)
    END IF

    !=========================================================================================================
    !=== write mri_flow at the end of the current cycle ======================================================
    !=========================================================================================================
    !if (phase.eq.intervals) then
    !   mri%velocity_mean%array = mri%velocity_mean%array*U_ref
    !   mri%velocity_cov%array  = mri%velocity_cov%array*U_ref*U_ref
    !   CALL num_to_string(8,write_stats_count/intervals+1,count_char)
    !   write_file = kalman_mri_output_file_path
    !   !write_file = trim(kalman_mri_output_file_path)//'_cycle.'//count_char//'.h5' 
    !   CALL mr_io_write_parallel_segmentedflow(MPI_COMM_WORLD, MPI_INFO_NULL, write_file, mri)
    !   CALL h5open_f(herror)
    !   mri%velocity_mean%array = mri%velocity_mean%array/U_ref
    !   mri%velocity_cov%array  = mri%velocity_cov%array/U_ref/U_ref
    !end if 
    !=========================================================================================================
 
  end if

  intervals = intervals/kalman_num_time_refinements

  write_stats_count = write_stats_count + 1
  dtime_out_scal = dtime_phases(phase)
  dtime_out_vect = dtime_out_scal
  time_out_scal = time_out_scal + dtime_out_scal
  write_out_scal = .FALSE.
 
  RETURN
 
  END SUBROUTINE compute_stats

 
  SUBROUTINE write_restart_stats
 
  USE mod_dims
  USE mod_vars
  USE mod_lib
  USE mod_diff
  USE mod_inout
  USE usr_func
  USE usr_vars
  USE MPI
  USE HDF5
  USE mr_io_protocol
  USE mr_io_parallel_spacetime
 
  IMPLICIT NONE
 
  !--- write time-dependent data after time integration for next restart ---
  ! this is just an example/template

  CHARACTER(LEN=8) :: count_char
  CHARACTER(LEN=3) :: next_restart_char
  CHARACTER(LEN=2) :: id,phase_char
  INTEGER :: i
  CHARACTER(LEN=300) :: write_file
 
  IF (rank == 0) WRITE(*,'(a,i3,a)') 'writing stats data for restart',restart,' ...'
 
  !========================================================================================================
  !=== neue Restart-Nr. als String fuer File-Namen ========================================================
  !========================================================================================================
  CALL num_to_string(3,restart,next_restart_char)
  !========================================================================================================

  IF (write_restart_yes.EQV..FALSE.) RETURN

  IF (associated(mri_flow)) then
    mri_flow%mri%velocity_mean%array = mri_flow%mri%velocity_mean%array*U_ref
    mri_flow%mri%velocity_cov%array  = mri_flow%mri%velocity_cov%array*U_ref*U_ref
    write_file = kalman_mri_output_file_path
    CALL mr_io_write_parallel_segmentedflow(MPI_COMM_WORLD, MPI_INFO_NULL, write_file, mri_flow%mri)
    CALL h5open_f(herror)
  END IF

  RETURN
 
  END SUBROUTINE write_restart_stats


  ! subroutine that reads data from previous restart before starting time integration 
  SUBROUTINE read_restart_stats
 
  USE mod_dims
  USE mod_vars
  USE mod_lib
  USE mod_diff
  USE mod_inout
  USE usr_func
  USE usr_vars
  USE MPI
  USE HDF5
  USE mr_io_protocol
  USE mr_io_parallel_spacetime
 
  !--- read time-dependent data from previous restart before time integration starts ---

  IMPLICIT NONE
 
  CHARACTER(LEN=2) :: id,phase_char
  INTEGER :: i
  CHARACTER(LEN=300) :: read_file
 
  IF (rank == 0) WRITE(*,'(a,i3,a)') 'reading stats data for restart',restart,' ...'
 
  !========================================================================================================
  !=== neue Restart-Nr. als String fuer File-Namen ========================================================
  !========================================================================================================
  CALL num_to_string(3,restart,restart_char)
  !========================================================================================================

  IF (associated(mri_flow)) then

    read_file = kalman_mri_output_file_path

    call mr_io_deallocate_dist_segmentedflow_mri_padded(mri_flow)
    nullify(mri_flow); allocate(mri_flow)
    mri_flow%domain_padding%lhs = mri_inst%domain_padding%lhs*kalman_num_spatial_refinements
    mri_flow%domain_padding%rhs = mri_inst%domain_padding%rhs*kalman_num_spatial_refinements
    CALL mr_io_read_parallel_segmentedflow_padded(MPI_COMM_WORLD, MPI_INFO_NULL, (/NB1,NB2,NB3/), read_file, mri_flow)
    CALL h5open_f(herror)

    mri_flow%mri%velocity_mean%array = mri_flow%mri%velocity_mean%array/U_ref
    mri_flow%mri%velocity_cov%array  = mri_flow%mri%velocity_cov%array/U_ref/U_ref

  END IF

  RETURN
 
  END SUBROUTINE read_restart_stats
