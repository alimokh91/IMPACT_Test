module mr_protocol

use hdf5

implicit none

! 
character(len=11) :: SpatialMRI_group_name = "spatial-mri"

type SpatialMRI
     ! voxel_feature
     real*8, dimension(:,:,:), allocatable :: voxel_feature
     integer, dimension(3) :: voxel_feature_dims
end type

character(len=14) :: SpaceTimeMRI_group_name = "space-time-mri"

type SpaceTimeMRI
     ! voxel_feature
     real*8, dimension(:,:,:,:,:), allocatable :: voxel_feature
     integer, dimension(5) :: voxel_feature_dims
end type

end module mr_protocol
