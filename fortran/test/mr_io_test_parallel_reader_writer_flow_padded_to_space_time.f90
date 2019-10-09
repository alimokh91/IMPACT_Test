program mr_io_test_parallel_reader_writer_flow_padded
  
    use mr_io_parallel_spacetime
    use mr_io_test_arg_parser
    use mpi

    implicit none

    INTEGER err
!    character(len=100) :: path = "mr_io_test_parallel_flow_with_padding.h5"

    type(DistFlowMRIPadded) :: mri_dest_padded
    type(DistSpacetimeMRI) :: mri_dest


    integer, dimension(5) :: voxel_feature_shape
    integer, dimension(5) :: voxel_feature_shape_for_loop

    integer :: i, ix, iy, iz, it, iv

!    integer :: gdb = 0
!
!    do while (gdb == 0)
!      call sleep(1)
!    end do

    call MPI_Init(err)

    call mr_io_test_parse_args_parallel_reader_writer_padded_to_st()
    
!    print *, domain_padding_lhs
!    print *, domain_padding_rhs

    mri_dest_padded%domain_padding%lhs = domain_padding_lhs
    mri_dest_padded%domain_padding%rhs = domain_padding_rhs

!    print *, "Padding... (lhs/rhs)"
!    print *, mri_dest_padded%domain_padding%lhs
!    print *, mri_dest_padded%domain_padding%rhs
!    call flush()

    call mr_io_read_parallel_flow_padded(MPI_COMM_WORLD, MPI_INFO_NULL, mr_io_test_mpi_cart_dims, path, mri_dest_padded)

    allocate(mri_dest%t_coordinates(size(mri_dest_padded%mri%t_coordinates)*simulation_time_refinement))
    allocate(mri_dest%x_coordinates(size(mri_dest_padded%mri%x_coordinates)*simulation_spatial_refinement(1)))
    allocate(mri_dest%y_coordinates(size(mri_dest_padded%mri%y_coordinates)*simulation_spatial_refinement(2)))
    allocate(mri_dest%z_coordinates(size(mri_dest_padded%mri%z_coordinates)*simulation_spatial_refinement(3)))
    allocate(mri_dest%voxel_feature%array(  lbound(mri_dest_padded%mri%velocity_mean%array,1):ubound(mri_dest_padded%mri%velocity_mean%array,1), &
                                           (lbound(mri_dest_padded%mri%velocity_mean%array,2)-1)*simulation_time_refinement+1:ubound(mri_dest_padded%mri%velocity_mean%array,2)*simulation_time_refinement, &
                                           (lbound(mri_dest_padded%mri%velocity_mean%array,3)-1)*simulation_spatial_refinement(1)+1:ubound(mri_dest_padded%mri%velocity_mean%array,3)*simulation_spatial_refinement(1), &
                                           (lbound(mri_dest_padded%mri%velocity_mean%array,4)-1)*simulation_spatial_refinement(2)+1:ubound(mri_dest_padded%mri%velocity_mean%array,4)*simulation_spatial_refinement(2), &
                                           (lbound(mri_dest_padded%mri%velocity_mean%array,5)-1)*simulation_spatial_refinement(3)+1:ubound(mri_dest_padded%mri%velocity_mean%array,5)*simulation_spatial_refinement(3) ))

    ! time
    do i=0,size(mri_dest_padded%mri%t_coordinates)-1
        mri_dest%t_coordinates(i*simulation_time_refinement+1:(i+1)*simulation_time_refinement) = mri_dest_padded%mri%t_coordinates(i+1)
    end do
    mri_dest%t_dim = mri_dest_padded%mri%t_dim*simulation_time_refinement

    ! geometry
    do i=0,size(mri_dest_padded%mri%x_coordinates)-1
        mri_dest%x_coordinates(i*simulation_spatial_refinement(1)+1:(i+1)*simulation_spatial_refinement(1)) = mri_dest_padded%mri%x_coordinates(i+1)
    end do
    mri_dest%x_dim = mri_dest_padded%mri%x_dim*simulation_spatial_refinement(1)

    do i=0,size(mri_dest_padded%mri%y_coordinates)-1
        mri_dest%y_coordinates(i*simulation_spatial_refinement(2)+1:(i+1)*simulation_spatial_refinement(2)) = mri_dest_padded%mri%y_coordinates(i+1)
    end do
    mri_dest%y_dim = mri_dest_padded%mri%y_dim*simulation_spatial_refinement(2)

    do i=0,size(mri_dest_padded%mri%z_coordinates)-1
        mri_dest%z_coordinates(i*simulation_spatial_refinement(3)+1:(i+1)*simulation_spatial_refinement(3)) = mri_dest_padded%mri%z_coordinates(i+1)
    end do
    mri_dest%z_dim = mri_dest_padded%mri%z_dim*simulation_spatial_refinement(3)

    ! voxel_feature
    do iz=(lbound(mri_dest_padded%mri%velocity_mean%array,5)-1),(ubound(mri_dest_padded%mri%velocity_mean%array,5)-1)
        do iy=(lbound(mri_dest_padded%mri%velocity_mean%array,4)-1),(ubound(mri_dest_padded%mri%velocity_mean%array,4)-1)
            do ix=(lbound(mri_dest_padded%mri%velocity_mean%array,3)-1),(ubound(mri_dest_padded%mri%velocity_mean%array,3)-1)
                do it=(lbound(mri_dest_padded%mri%velocity_mean%array,2)-1),(ubound(mri_dest_padded%mri%velocity_mean%array,2)-1)
                    do iv=1,3
                        mri_dest%voxel_feature%array(iv, &
                                             it*simulation_time_refinement+1:(it+1)*simulation_time_refinement, &
                                             ix*simulation_spatial_refinement(1)+1:(ix+1)*simulation_spatial_refinement(1), &
                                             iy*simulation_spatial_refinement(2)+1:(iy+1)*simulation_spatial_refinement(2), &
                                             iz*simulation_spatial_refinement(3)+1:(iz+1)*simulation_spatial_refinement(3)) &
                                             = mri_dest_padded%mri%velocity_mean%array(iv,it+1,ix+1,iy+1,iz+1)
                    end do
                end do
            end do
        end do
    end do

    mri_dest%voxel_feature%time_offset = mri_dest_padded%mri%velocity_mean%time_offset*simulation_time_refinement
    mri_dest%voxel_feature%time_dim = mri_dest_padded%mri%velocity_mean%time_dim*simulation_time_refinement
    mri_dest%voxel_feature%offset = (/ mri_dest_padded%mri%velocity_mean%offset(1)*simulation_spatial_refinement(1), &
                                       mri_dest_padded%mri%velocity_mean%offset(2)*simulation_spatial_refinement(2), &
                                       mri_dest_padded%mri%velocity_mean%offset(3)*simulation_spatial_refinement(3) /)
    mri_dest%voxel_feature%dims = (/ mri_dest_padded%mri%velocity_mean%dims(1)*simulation_spatial_refinement(1), &
                                     mri_dest_padded%mri%velocity_mean%dims(2)*simulation_spatial_refinement(2), &
                                     mri_dest_padded%mri%velocity_mean%dims(3)*simulation_spatial_refinement(3) /)

    print *, SpaceTimeMRI_group_name

    print *, shape(mri_dest%t_coordinates)
    print *, mri_dest%t_coordinates

    print *, shape(mri_dest%x_coordinates)
    print *, mri_dest%x_coordinates

    print *, shape(mri_dest%y_coordinates)
    print *, mri_dest%y_coordinates

    print *, shape(mri_dest%z_coordinates)
    print *, mri_dest%z_coordinates

    print *, mri_dest%voxel_feature%dims
    print *, mri_dest%voxel_feature%offset
    print *, voxel_feature_shape(3:5)

    print *, mri_dest%voxel_feature%time_dim
    print *, mri_dest%voxel_feature%time_offset
    print *, voxel_feature_shape(2)

    print *, voxel_feature_shape(1)
    print *, mri_dest%voxel_feature%array

    call mr_io_write_parallel_spacetime(MPI_COMM_WORLD, MPI_INFO_NULL, out_path, mri_dest)

    call MPI_Finalize(err)
        
end program mr_io_test_parallel_reader_writer_flow_padded
