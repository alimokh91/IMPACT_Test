module mr_protocol

use hdf5

implicit none

! 
character(len=3) :: MRI_group_name = "mri"

type MRI
     ! voxel_feature
     real*8, dimension(:,:), allocatable :: voxel_feature
     integer, dimension(2) :: voxel_feature_dims
end type

end module mr_protocol
