module mr

use hdf5

implicit none

type MRI
     integer, dimension(2,2) :: voxel_feature
end type

integer(HSIZE_T), dimension(2) :: voxel_feature_dims = (/2,2/)

end module mr
