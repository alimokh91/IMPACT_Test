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

type DistSpacetimeFeature ! (vect-comp, t, x, y, z)
     real*8, dimension(:,:,:,:,:), allocatable :: array
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
     real*8, dimension(:,:,:,:,:), allocatable :: velocity_mean
     integer, dimension(5) :: velocity_mean_dims
     real*8, dimension(:,:,:,:,:), allocatable :: velocity_cov
     integer, dimension(5) :: velocity_cov_dims
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
     type(DistSpacetimeFeature) :: velocity_mean
     type(DistSpacetimeFeature) :: velocity_cov
end type


end module mr_protocol
