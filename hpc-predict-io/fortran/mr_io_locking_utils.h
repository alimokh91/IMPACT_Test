#include "hdf5.h"

// POSIX-level file locking functions
int acquire_shared_lock_blocking(int fd);
int release_lock(int fd);
int acquire_exclusive_lock_blocking(int fd);

// C standard stream I/O
FILE * mr_io_stream_reader_open(char *name);
void mr_io_stream_reader_close(FILE *reader_stream);

FILE * mr_io_stream_writer_open(char *name);
void mr_io_stream_writer_close(FILE *writer_stream);

// HDF5 I/O
hid_t mr_io_h5_reader_open(char *name);
herr_t mr_io_h5_reader_close(hid_t reader_file);

hid_t mr_io_h5_writer_open(char *name);
herr_t mr_io_h5_writer_close(hid_t writer_file);

// Parallel HDF5 I/O
hid_t mr_io_h5_parallel_reader_open(MPI_Comm mr_io_mpi_comm, MPI_Info mr_io_mpi_info, char *name);
herr_t mr_io_h5_parallel_reader_close(hid_t h5_file);
hid_t mr_io_h5_parallel_writer_open(MPI_Comm mr_io_mpi_comm, MPI_Info mr_io_mpi_info, char *name);
herr_t mr_io_h5_parallel_writer_close(hid_t h5_file);
