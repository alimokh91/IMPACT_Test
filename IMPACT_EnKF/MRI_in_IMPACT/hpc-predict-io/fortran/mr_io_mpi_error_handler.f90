!module mr_io_mpi_error_handler

#ifdef __GFORTRAN__

#define mr_io_handle_error(err) if(err /= 0) CALL mr_io_handle_hdf5_error(err);
#define mr_io_handle_arg_error(err) if(err /= 0) CALL mr_io_handle_argument_error(err);
#define mr_io_sget_simple_extent_dims_handle_error(err) if(err == -1) CALL mr_io_handle_hdf5_error(err);

#define mr_io_handle_mpi_error(err) if(err /= 0) CALL mr_io_handle_mpi_error_(err);

subroutine mr_io_handle_mpi_error_(error)

  implicit none

  INTEGER     ::   error          ! Error flag

  if(error /= 0) then
    write(*,*) "MPI error ",error,"- printing backtrace and aborting..."
    call backtrace
    !call abort
  end if

end subroutine mr_io_handle_mpi_error_

#endif

!end module mr_io_mpi_error_handler
