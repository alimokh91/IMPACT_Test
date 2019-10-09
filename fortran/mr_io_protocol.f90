module mr_io_protocol

use hdf5

implicit none
public

type DomainPadding
     integer, dimension(3) :: lhs
     integer, dimension(3) :: rhs
end type

! ************************ SpatialMRI ************************

character(len=100) :: SpatialMRI_group_name = "spatial-mri"

type SpatialMRI
     ! voxel_feature
     real*8, dimension(:,:,:), allocatable :: voxel_feature
     integer, dimension(3) :: voxel_feature_dims
end type

type DistSpatialFeature
     real*8, dimension(:,:,:), allocatable :: array 
     integer, dimension(3) :: offset = (/ -1, -1, -1 /) ! file_offset
     integer, dimension(3) :: dims = (/ -1, -1, -1 /) ! file_dim     
end type

type DistSpatialMRI
     ! voxel_feature
     type(DistSpatialFeature) :: voxel_feature
end type


! ************************ SpaceTimeMRI ************************

character(len=100) :: SpaceTimeMRI_group_name = "space-time-mri"

type SpaceTimeMRI
     ! time
     real*8, dimension(:), allocatable :: t_coordinates
     integer :: t_dim

     ! geometry
     real*8, dimension(:), allocatable :: x_coordinates
     integer :: x_dim
     real*8, dimension(:), allocatable :: y_coordinates
     integer :: y_dim
     real*8, dimension(:), allocatable :: z_coordinates
     integer :: z_dim

     ! voxel_feature
     real*8, dimension(:,:,:,:,:), allocatable :: voxel_feature
     integer, dimension(5) :: voxel_feature_dims
end type

type DistSpacetimeScalarFeature ! (t, x, y, z)
     real*8, dimension(:,:,:,:), allocatable :: array
     integer :: time_offset = -1
     integer :: time_dim = -1
     integer, dimension(3) :: offset = (/ -1, -1, -1 /) ! file_offset
     integer, dimension(3) :: dims = (/ -1, -1, -1 /) ! file_dim
end type

type DistSpacetimeFeature ! (vect-comp, t, x, y, z)
     real*8, dimension(:,:,:,:,:), allocatable :: array
     integer :: time_offset = -1
     integer :: time_dim = -1
     integer, dimension(3) :: offset = (/ -1, -1, -1 /) ! file_offset
     integer, dimension(3) :: dims = (/ -1, -1, -1 /) ! file_dim
end type

type DistSpacetimeMatrixFeature ! (vect-comp, t, x, y, z)
     real*8, dimension(:,:,:,:,:,:), allocatable :: array
     integer :: time_offset = -1
     integer :: time_dim = -1
     integer, dimension(3) :: offset = (/ -1, -1, -1 /) ! file_offset
     integer, dimension(3) :: dims = (/ -1, -1, -1 /) ! file_dim
end type

type DistSpacetimeMRI
     ! time
     real*8, dimension(:), allocatable :: t_coordinates
     integer :: t_dim

     ! geometry
     real*8, dimension(:), allocatable :: x_coordinates
     integer :: x_dim
     real*8, dimension(:), allocatable :: y_coordinates
     integer :: y_dim
     real*8, dimension(:), allocatable :: z_coordinates
     integer :: z_dim

     ! voxel_feature
     type(DistSpacetimeFeature) :: voxel_feature
end type


type DistSpacetimeMRIPadded
     type(DistSpacetimeMRI) :: mri
     type(DomainPadding) :: domain_padding
end type


! ************************ FlowMRI ************************

character(len=100) :: FlowMRI_group_name = "flow-mri"

type FlowMRI
     ! time
     real*8, dimension(:), allocatable :: t_coordinates
     integer :: t_dim

     ! geometry
     real*8, dimension(:), allocatable :: x_coordinates
     integer :: x_dim
     real*8, dimension(:), allocatable :: y_coordinates
     integer :: y_dim
     real*8, dimension(:), allocatable :: z_coordinates
     integer :: z_dim

     ! velocity mean and covariance
     real*8, dimension(:,:,:,:), allocatable :: intensity
     integer, dimension(4) :: intensity_dims
     real*8, dimension(:,:,:,:,:), allocatable :: velocity_mean
     integer, dimension(5) :: velocity_mean_dims
     real*8, dimension(:,:,:,:,:,:), allocatable :: velocity_cov
     integer, dimension(6) :: velocity_cov_dims
end type

type DistFlowMRI
     ! time
     real*8, dimension(:), allocatable :: t_coordinates
     integer :: t_dim

     ! geometry
     real*8, dimension(:), allocatable :: x_coordinates
     integer :: x_dim
     real*8, dimension(:), allocatable :: y_coordinates
     integer :: y_dim
     real*8, dimension(:), allocatable :: z_coordinates
     integer :: z_dim

     ! velocity mean and covariance
     type(DistSpacetimeScalarFeature) :: intensity
     type(DistSpacetimeFeature) :: velocity_mean
     type(DistSpacetimeMatrixFeature) :: velocity_cov
end type

type DistFlowMRIPadded
     type(DistFlowMRI) :: mri
     type(DomainPadding) :: domain_padding
end type


! ************************ SegmentedFlowMRI ************************

character(len=100) :: SegmentedFlowMRI_group_name = "segmented-flow-mri"

type SegmentedFlowMRI
     ! time
     real*8, dimension(:), allocatable :: t_coordinates
     integer :: t_dim

     ! geometry
     real*8, dimension(:), allocatable :: x_coordinates
     integer :: x_dim
     real*8, dimension(:), allocatable :: y_coordinates
     integer :: y_dim
     real*8, dimension(:), allocatable :: z_coordinates
     integer :: z_dim

     ! velocity mean and covariance
     real*8, dimension(:,:,:,:), allocatable :: intensity
     integer, dimension(4) :: intensity_dims
     real*8, dimension(:,:,:,:,:), allocatable :: velocity_mean
     integer, dimension(5) :: velocity_mean_dims
     real*8, dimension(:,:,:,:,:,:), allocatable :: velocity_cov
     integer, dimension(6) :: velocity_cov_dims
     real*8, dimension(:,:,:,:), allocatable :: segmentation_prob
     integer, dimension(4) :: segmentation_prob_dims
end type

type DistSegmentedFlowMRI
     ! time
     real*8, dimension(:), allocatable :: t_coordinates
     integer :: t_dim

     ! geometry
     real*8, dimension(:), allocatable :: x_coordinates
     integer :: x_dim
     real*8, dimension(:), allocatable :: y_coordinates
     integer :: y_dim
     real*8, dimension(:), allocatable :: z_coordinates
     integer :: z_dim

     ! velocity mean and covariance
     type(DistSpacetimeScalarFeature) :: intensity
     type(DistSpacetimeFeature) :: velocity_mean
     type(DistSpacetimeMatrixFeature) :: velocity_cov
     type(DistSpacetimeScalarFeature) :: segmentation_prob
end type



contains

! Deallocation subroutines

subroutine mr_io_deallocate_spatial_mri(mri)

    implicit none
    type(SpatialMRI), intent(inout) :: mri
    deallocate(mri%voxel_feature)

end subroutine mr_io_deallocate_spatial_mri

subroutine mr_io_deallocate_dist_spatial_mri(mri)

    implicit none
    type(DistSpatialMRI), intent(inout) :: mri
    deallocate(mri%voxel_feature%array)

end subroutine mr_io_deallocate_dist_spatial_mri


subroutine mr_io_deallocate_spacetime_mri(mri)

    implicit none
    type(SpaceTimeMRI), intent(inout) :: mri
    deallocate(mri%x_coordinates)
    deallocate(mri%y_coordinates)
    deallocate(mri%z_coordinates)
    deallocate(mri%t_coordinates)
    deallocate(mri%voxel_feature)

end subroutine mr_io_deallocate_spacetime_mri

subroutine mr_io_deallocate_dist_spacetime_mri(mri)

    implicit none
    type(DistSpacetimeMRI), intent(inout) :: mri
    deallocate(mri%x_coordinates)
    deallocate(mri%y_coordinates)
    deallocate(mri%z_coordinates)
    deallocate(mri%t_coordinates)
    deallocate(mri%voxel_feature%array)

end subroutine mr_io_deallocate_dist_spacetime_mri

subroutine mr_io_deallocate_flow_mri(mri)

    implicit none
    type(FlowMRI), intent(inout) :: mri
    deallocate(mri%x_coordinates)
    deallocate(mri%y_coordinates)
    deallocate(mri%z_coordinates)
    deallocate(mri%t_coordinates)
    deallocate(mri%intensity)
    deallocate(mri%velocity_mean)
    deallocate(mri%velocity_cov)

end subroutine mr_io_deallocate_flow_mri

subroutine mr_io_deallocate_dist_flow_mri(mri)

    implicit none
    type(DistFlowMRI), intent(inout) :: mri
    deallocate(mri%x_coordinates)
    deallocate(mri%y_coordinates)
    deallocate(mri%z_coordinates)
    deallocate(mri%t_coordinates)
    deallocate(mri%intensity%array)
    deallocate(mri%velocity_mean%array)
    deallocate(mri%velocity_cov%array)

end subroutine mr_io_deallocate_dist_flow_mri

subroutine mr_io_deallocate_dist_flow_mri_padded(mri)

    implicit none
    type(DistFlowMRIPadded), intent(inout) :: mri
    call mr_io_deallocate_dist_flow_mri(mri%mri)

end subroutine mr_io_deallocate_dist_flow_mri_padded

subroutine mr_io_deallocate_segmentedflow_mri(mri)

    implicit none
    type(SegmentedFlowMRI), intent(inout) :: mri
    deallocate(mri%x_coordinates)
    deallocate(mri%y_coordinates)
    deallocate(mri%z_coordinates)
    deallocate(mri%t_coordinates)
    deallocate(mri%intensity)
    deallocate(mri%velocity_mean)
    deallocate(mri%velocity_cov)
    deallocate(mri%segmentation_prob)

end subroutine mr_io_deallocate_segmentedflow_mri

subroutine mr_io_deallocate_dist_segmentedflow_mri(mri)

    implicit none
    type(DistSegmentedFlowMRI), intent(inout) :: mri
    deallocate(mri%x_coordinates)
    deallocate(mri%y_coordinates)
    deallocate(mri%z_coordinates)
    deallocate(mri%t_coordinates)
    deallocate(mri%intensity%array)
    deallocate(mri%velocity_mean%array)
    deallocate(mri%velocity_cov%array)
    deallocate(mri%segmentation_prob%array)

end subroutine mr_io_deallocate_dist_segmentedflow_mri


end module mr_io_protocol
