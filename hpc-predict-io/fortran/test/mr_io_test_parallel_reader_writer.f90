program mr_io_test_parallel_reader_writer
  
    use mr_io_parallel
    use mr_io_test_arg_parser
    use mpi

    implicit none

    INTEGER err
!    character(len=30) :: in_path = "mr_io_test_parallel_in.h5"
!    character(len=30) :: out_path = "mr_io_test_parallel_out.h5"
    type(DistSpatialMRI) :: mri_dest    

!    INTEGER :: gdb = 0
!    integer rank
!
!    call MPI_Init(err)
!
!    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, err)
!    if (rank == 0) then
!    print *, "Rank 0 is running on process ID is ", getpid()
!    do while (gdb == 0)
!        call sleep(2)
!    end do
!    end if
!
!    ! End of debugging

    call MPI_Init(err)
    
    call mr_io_test_parse_args_parallel_reader_writer()

    call mr_io_read_parallel_spatial(MPI_COMM_WORLD, MPI_INFO_NULL, mr_io_test_mpi_cart_dims, in_path, mri_dest)
    print *, "Array properties:"
    print *, mri_dest%scalar_feature%dims
    print *, mri_dest%scalar_feature%offset
    print *, shape(mri_dest%scalar_feature%array)
    call mr_io_write_parallel_spatial(MPI_COMM_WORLD, MPI_INFO_NULL, out_path, mri_dest)

    call MPI_Finalize(err)
        
end program mr_io_test_parallel_reader_writer
