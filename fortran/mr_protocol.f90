module mr_protocol

use hdf5

implicit none


! ************************ SpatialMRI ************************

character(len=11) :: SpatialMRI_group_name = "spatial-mri"

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

character(len=14) :: SpaceTimeMRI_group_name = "space-time-mri"

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


! ************************ HPCPredictMRI ************************

character(len=16) :: HPCPredictMRI_group_name = "hpc-predict-mri"

type HPCPredictMRI
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

type DistHPCPredictMRI
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

! ************************ SegmentedHPCPredictMRI ************************

character(len=26) :: SegmentedHPCPredictMRI_group_name = "segmented-hpc-predict-mri"

type SegmentedHPCPredictMRI
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

type DistSegmentedHPCPredictMRI
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

subroutine mr_io_deallocate_hpcpredict_mri(mri)

    implicit none
    type(HPCPredictMRI), intent(inout) :: mri
    deallocate(mri%x_coordinates)
    deallocate(mri%y_coordinates)
    deallocate(mri%z_coordinates)
    deallocate(mri%t_coordinates)
    deallocate(mri%intensity)
    deallocate(mri%velocity_mean)
    deallocate(mri%velocity_cov)

end subroutine mr_io_deallocate_hpcpredict_mri

subroutine mr_io_deallocate_dist_hpcpredict_mri(mri)

    implicit none
    type(DistHPCPredictMRI), intent(inout) :: mri
    deallocate(mri%x_coordinates)
    deallocate(mri%y_coordinates)
    deallocate(mri%z_coordinates)
    deallocate(mri%t_coordinates)
    deallocate(mri%intensity%array)
    deallocate(mri%velocity_mean%array)
    deallocate(mri%velocity_cov%array)

end subroutine mr_io_deallocate_dist_hpcpredict_mri

subroutine mr_io_deallocate_segmentedhpcpredict_mri(mri)

    implicit none
    type(SegmentedHPCPredictMRI), intent(inout) :: mri
    deallocate(mri%x_coordinates)
    deallocate(mri%y_coordinates)
    deallocate(mri%z_coordinates)
    deallocate(mri%t_coordinates)
    deallocate(mri%intensity)
    deallocate(mri%velocity_mean)
    deallocate(mri%velocity_cov)
    deallocate(mri%segmentation_prob)

end subroutine mr_io_deallocate_segmentedhpcpredict_mri

subroutine mr_io_deallocate_dist_segmentedhpcpredict_mri(mri)

    implicit none
    type(DistSegmentedHPCPredictMRI), intent(inout) :: mri
    deallocate(mri%x_coordinates)
    deallocate(mri%y_coordinates)
    deallocate(mri%z_coordinates)
    deallocate(mri%t_coordinates)
    deallocate(mri%intensity%array)
    deallocate(mri%velocity_mean%array)
    deallocate(mri%velocity_cov%array)
    deallocate(mri%segmentation_prob%array)

end subroutine mr_io_deallocate_dist_segmentedhpcpredict_mri


end module mr_protocol
