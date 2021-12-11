module mr_io_locking_utils

use hdf5

!#ifdef f2003
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit
!#else
!#define stdin  5
!#define stdout 6
!#define stderr 0
!#endif

!#define stdin  5
!#define stdout 6
!#define stderr 0

implicit none
private
public :: mr_io_h5_reader_open_f
public :: mr_io_h5_reader_close_f
public :: mr_io_h5_writer_open_f
public :: mr_io_h5_writer_close_f

public :: mr_io_h5_parallel_reader_open_f
public :: mr_io_h5_parallel_reader_close_f
public :: mr_io_h5_parallel_writer_open_f
public :: mr_io_h5_parallel_writer_close_f

public :: stdin
public :: stdout
public :: stderr
public :: mr_io_lock_stdout_stderr
public :: mr_io_unlock_stdout_stderr

! Declare interface for POSIX fsync
interface
function fsync (fd) bind(c,name="fsync")
use iso_c_binding, only: c_int
  implicit none
  integer(c_int), value :: fd
  integer(c_int) :: fsync
end function fsync

! POSIX-system-level I/O
function acquire_shared_lock_blocking (fd) bind(c,name="acquire_shared_lock_blocking")
use iso_c_binding, only: c_int
  implicit none
  integer(c_int), value :: fd
  integer(c_int) :: acquire_shared_lock_blocking
end function acquire_shared_lock_blocking

function acquire_exclusive_lock_blocking (fd) bind(c,name="acquire_exclusive_lock_blocking")
use iso_c_binding, only: c_int
  implicit none
  integer(c_int), value :: fd
  integer(c_int) :: acquire_exclusive_lock_blocking
end function acquire_exclusive_lock_blocking

function release_lock (fd) bind(c,name="release_lock")
use iso_c_binding, only: c_int
  implicit none
  integer(c_int), value :: fd
  integer(c_int) :: release_lock
end function release_lock

! Sequential HDF5-I/O
function mr_io_h5_reader_open(name) bind(c,name="mr_io_h5_reader_open")
use iso_c_binding, only: c_char, c_long
  character(kind=c_char), dimension(*) :: name
  integer(c_long) :: mr_io_h5_reader_open
end function mr_io_h5_reader_open

function mr_io_h5_reader_close(h5_file) bind(c,name="mr_io_h5_reader_close")
use iso_c_binding, only: c_long, c_int
  integer(c_long), value :: h5_file
  integer(c_int) :: mr_io_h5_reader_close
end function mr_io_h5_reader_close

function mr_io_h5_writer_open(name) bind(c,name="mr_io_h5_writer_open")
use iso_c_binding, only: c_char, c_long
  character(kind=c_char), dimension(*) :: name
  integer(c_long) :: mr_io_h5_writer_open
end function mr_io_h5_writer_open

function mr_io_h5_writer_close(h5_file) bind(c,name="mr_io_h5_writer_close")
use iso_c_binding, only: c_long, c_int
  integer(c_long), value :: h5_file
  integer(c_int) :: mr_io_h5_writer_close
end function mr_io_h5_writer_close

! Parallel HDF5-I/O
function mr_io_h5_parallel_reader_open(mr_io_mpi_comm, mr_io_mpi_info, name) bind(c,name="mr_io_h5_parallel_reader_open")
use iso_c_binding, only: c_char, c_int, c_long
  character(kind=c_char), dimension(*) :: name
  integer(c_int), value :: mr_io_mpi_comm
  integer(c_int), value :: mr_io_mpi_info
  integer(c_long) :: mr_io_h5_parallel_reader_open
end function mr_io_h5_parallel_reader_open

function mr_io_h5_parallel_reader_close(h5_file) bind(c,name="mr_io_h5_parallel_reader_close")
use iso_c_binding, only: c_long, c_int
  integer(c_long), value :: h5_file
  integer(c_int) :: mr_io_h5_parallel_reader_close
end function mr_io_h5_parallel_reader_close

function mr_io_h5_parallel_writer_open(mr_io_mpi_comm, mr_io_mpi_info, name) bind(c,name="mr_io_h5_parallel_writer_open")
use iso_c_binding, only: c_char, c_int, c_long
  character(kind=c_char), dimension(*) :: name
  integer(c_int), value :: mr_io_mpi_comm
  integer(c_int), value :: mr_io_mpi_info
  integer(c_long) :: mr_io_h5_parallel_writer_open
end function mr_io_h5_parallel_writer_open

function mr_io_h5_parallel_writer_close(h5_file) bind(c,name="mr_io_h5_parallel_writer_close")
use iso_c_binding, only: c_long, c_int
  integer(c_long), value :: h5_file
  integer(c_int) :: mr_io_h5_parallel_writer_close
end function mr_io_h5_parallel_writer_close


end interface

contains

! TODO: Write stdout/stderr locking functions for tests implementation

function mr_io_lock_stdout_stderr()
  integer :: st, mr_io_lock_stdout_stderr

  st = acquire_exclusive_lock_blocking(Fnum(stdout))
  st = acquire_exclusive_lock_blocking(Fnum(stderr))

  mr_io_lock_stdout_stderr = st
end function mr_io_lock_stdout_stderr

function mr_io_unlock_stdout_stderr()
  integer :: st, mr_io_unlock_stdout_stderr

  flush(stdout)
  st = fsync(Fnum(stdout))
  st = release_lock(Fnum(stdout))

  flush(stderr)
  st = fsync(Fnum(stderr))
  st = release_lock(Fnum(stderr))

  mr_io_unlock_stdout_stderr = st
end function mr_io_unlock_stdout_stderr


! ***** Sequential reading/writing of HDF5 *****

function mr_io_h5_reader_open_f(path)

  use iso_c_binding, only: c_null_char
  implicit none

  character(len=*), intent(in) :: path
  character(len=len(path)+1) :: c_filename
  integer(hid_t) :: mr_io_h5_reader_open_f    ! File identifier

  c_filename = trim(path)//c_null_char
  mr_io_h5_reader_open_f = mr_io_h5_reader_open(c_filename);

end function mr_io_h5_reader_open_f


function mr_io_h5_reader_close_f(h5_file)
  implicit none

  integer(hid_t), intent(in) :: h5_file
  integer :: mr_io_h5_reader_close_f

  mr_io_h5_reader_close_f = mr_io_h5_reader_close(h5_file)

end function mr_io_h5_reader_close_f


function mr_io_h5_writer_open_f(path)

  use iso_c_binding, only: c_null_char
  implicit none

  character(len=*), intent(in) :: path
  character(len=len(path)+1) :: c_filename
  integer(hid_t) :: mr_io_h5_writer_open_f    ! File identifier

  c_filename = trim(path)//c_null_char
  mr_io_h5_writer_open_f = mr_io_h5_writer_open(c_filename);

end function mr_io_h5_writer_open_f


function mr_io_h5_writer_close_f(h5_file)
  implicit none

  integer(hid_t), intent(in) :: h5_file
  integer :: mr_io_h5_writer_close_f

  mr_io_h5_writer_close_f = mr_io_h5_writer_close(h5_file)

end function mr_io_h5_writer_close_f


! ***** Parallel reading/writing of HDF5 *****

function mr_io_h5_parallel_reader_open_f(mr_io_mpi_comm, mr_io_mpi_info, path)
  
  use iso_c_binding, only: c_null_char
  implicit none

  integer :: mr_io_mpi_comm
  integer :: mr_io_mpi_info
  character(len=*), intent(in) :: path
  character(len=len(path)+1) :: c_filename
  integer(hid_t) :: mr_io_h5_parallel_reader_open_f    ! File identifier
  !write(0,*) "Entering mr_io_h5_parallel_reader_open_f - c-f90 binding"
  c_filename = trim(path)//c_null_char
  mr_io_h5_parallel_reader_open_f = mr_io_h5_parallel_reader_open(mr_io_mpi_comm, mr_io_mpi_info, c_filename);

end function mr_io_h5_parallel_reader_open_f


function mr_io_h5_parallel_reader_close_f(h5_file)
  implicit none

  integer(hid_t), intent(in) :: h5_file
  integer :: mr_io_h5_parallel_reader_close_f

  mr_io_h5_parallel_reader_close_f = mr_io_h5_parallel_reader_close(h5_file)

end function mr_io_h5_parallel_reader_close_f


function mr_io_h5_parallel_writer_open_f(mr_io_mpi_comm, mr_io_mpi_info, path)

  use iso_c_binding, only: c_null_char
  implicit none

  integer :: mr_io_mpi_comm
  integer :: mr_io_mpi_info
  character(len=*), intent(in) :: path
  character(len=len(path)+1) :: c_filename
  integer(hid_t) :: mr_io_h5_parallel_writer_open_f    ! File identifier

  c_filename = trim(path)//c_null_char
  mr_io_h5_parallel_writer_open_f = mr_io_h5_parallel_writer_open(mr_io_mpi_comm, mr_io_mpi_info, c_filename);

end function mr_io_h5_parallel_writer_open_f


function mr_io_h5_parallel_writer_close_f(h5_file)
  implicit none

  integer(hid_t), intent(in) :: h5_file
  integer :: mr_io_h5_parallel_writer_close_f

  mr_io_h5_parallel_writer_close_f = mr_io_h5_parallel_writer_close(h5_file)

end function mr_io_h5_parallel_writer_close_f



end module mr_io_locking_utils
