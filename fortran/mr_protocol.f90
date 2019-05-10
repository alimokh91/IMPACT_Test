module mr_protocol

use hdf5

implicit none

! 
type MRI
     ! voxel_feature
     real, dimension(:,:), allocatable :: voxel_feature
     integer, dimension(2) :: voxel_feature_dims
end type

end module mr_protocol
