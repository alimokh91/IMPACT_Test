module mr_protocol

use hdf5

implicit none

! 
character(len=6) :: MRI_group_name = "3d-mri"

type MRI
     ! voxel_feature
     real*8, dimension(:,:,:), allocatable :: voxel_feature
     integer, dimension(3) :: voxel_feature_dims
end type

end module mr_protocol
