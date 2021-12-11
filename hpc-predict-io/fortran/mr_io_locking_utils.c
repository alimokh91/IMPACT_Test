#include <stdio.h>
#include <sys/file.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "mr_io_locking_utils.h"

// POSIX-level file locking functions
int acquire_shared_lock_blocking(int fd);
int release_lock(int fd);
int acquire_exclusive_lock_blocking(int fd);

//int acquire_shared_lock_blocking_(int *fd) { return acquire_shared_lock_blocking(*fd); }
//int release_lock_(int *fd) { return release_lock(*fd); }
//int acquire_exclusive_lock_blocking_(int *fd) { return acquire_exclusive_lock_blocking(*fd); }

int modify_lock(int fd, int op, char *error_msg) {
    int ret = flock(fd, op);
    if (ret == -1) {
        printf("%s (errno=%d, %s)", error_msg, errno, strerror(errno)); fflush(stdout);
        exit(-1);
    }
    return ret;
}

int acquire_shared_lock_blocking(int fd) {
    // Ensure writer before reader access
    off_t sz = 0;
    while ( (sz = lseek(fd, 0L, SEEK_END)) == 0) {}
    if (sz == (off_t) -1) {
        printf("Failed to lseek the end of file with descriptor %d (errno=%d, %s)", fd, errno, strerror(errno)); fflush(stdout);
        exit(-1);
    }
    lseek(fd, 0L, SEEK_SET);

    int ret = modify_lock(fd, LOCK_SH /*| LOCK_NB*/, "Error acquiring shared lock");
    //printf("PID %ld: requesting %s\n", (long) getpid(), "LOCK_SH"); fflush(stdout);
    return ret;
}

int release_lock(int fd) {
    //printf("PID %ld: releasing LOCK\n", (long) getpid()); fflush(stdout);
    return modify_lock(fd, LOCK_UN, "Error releasing lock");
}

int acquire_exclusive_lock_blocking(int fd) {
    int ret = modify_lock(fd, LOCK_EX /*| LOCK_NB*/, "Error acquiring shared lock");
    //printf("PID %ld: requesting %s\n", (long) getpid(), "LOCK_EX"); fflush(stdout);
    return ret;
}

// POSIX-level I/O (with file descriptors)

int mr_io_file_reader_open(char *fname);
void mr_io_file_reader_close(int fd);

int mr_io_file_writer_open(char *fname);
void mr_io_file_writer_close(int fd);

int mr_io_file_reader_open(char *fname){
    // Make sure input file exists before accessing it
    int access_ret = -1;
    while( (access_ret = access(fname, F_OK)) != 0) {
        if (errno != ENOENT) {
            printf("Failed to wait for file %s to be created (errno=%d, %s).\n", fname, errno, strerror(errno));
            exit(-1);
        }
    };

    int fd_in = open(fname, O_RDONLY);               /* Open file to be locked */
    if (fd_in == -1) {
        printf("Failed to open file %s.\n", fname);
        exit(-1);
    }
    acquire_shared_lock_blocking(fd_in);
    return fd_in;
}

void mr_io_file_reader_close(int reader_fd){
    release_lock(reader_fd);
    close(reader_fd);
}

int mr_io_file_writer_open(char *fname) {
    int fd_out = open(fname, O_WRONLY | O_CREAT | O_EXCL, 0600);               /* Open file to be locked */
    if (fd_out == -1) {
        printf("Failed to open file %s.\n", fname);
        exit(-1);
    }
    acquire_exclusive_lock_blocking(fd_out);
    return fd_out;
}

void mr_io_file_writer_close(int writer_fd) {
    fsync(writer_fd);
    release_lock(writer_fd);
    close(writer_fd);
}

// C standard stream I/O
//typedef struct {
//    char *name;
//    FILE *stream;
//} mr_io_stream;
//
//void mr_io_stream_reader_open(mr_io_stream *reader_stream) {
//    int reader_fd = mr_io_file_reader_open(reader_stream->name);
//    reader_stream->stream = fdopen(reader_fd, "r");
//}
//
//void mr_io_stream_reader_close(mr_io_stream *reader_stream) {
//    mr_io_file_reader_close(fileno(reader_stream->stream));
//}
//
//void mr_io_stream_writer_open(mr_io_stream *writer_stream) {
//    int writer_fd = mr_io_file_writer_open(writer_stream->name);
//    writer_stream->stream = fdopen(writer_fd, "w");
//}
//
//void mr_io_stream_writer_close(mr_io_stream *writer_stream) {
//    fflush(writer_stream->stream);
//    mr_io_file_writer_close(fileno(writer_stream->stream));
//}

FILE * mr_io_stream_reader_open(char *name) {
    return fdopen(mr_io_file_reader_open(name), "r");
}

void mr_io_stream_reader_close(FILE *stream) {
    mr_io_file_reader_close(fileno(stream));
}

FILE * mr_io_stream_writer_open(char *name) {
    return fdopen(mr_io_file_writer_open(name), "w");
}

void mr_io_stream_writer_close(FILE *stream) {
    fflush(stream);
    mr_io_file_writer_close(fileno(stream));
}


// HDF5 I/O

void h5_print_file_driver(hid_t plist_id) {
    // debugging
    printf("H5FD (file driver): %d",H5Pget_driver( plist_id ));
#ifdef H5FD_SEC2
    printf("POSIX: %d", 	H5FD_SEC2);
#endif
#ifdef H5FD_DIRECT
    printf("Direct: %d",	 	H5FD_DIRECT);
#endif
#ifdef H5FD_LOG
    printf("Log: %d",	 	H5FD_LOG);
#endif
#ifdef H5FD_WINDOWS
    printf("Windows: %d",	 	H5FD_WINDOWS);
#endif
#ifdef H5FD_STDIO
    printf("STDIO: %d",	 	H5FD_STDIO);
#endif
#ifdef H5FD_CORE
    printf("Memory: %d",	 	H5FD_CORE);
#endif
#ifdef H5FD_FAMILY
    printf("Family: %d",	 	H5FD_FAMILY);
#endif
#ifdef H5FD_MULTI
    printf("Multi: %d",	 	H5FD_MULTI);
#endif
#ifdef H5FD_SPLIT
    printf("Split: %d",	 	H5FD_SPLIT);
#endif
#ifdef H5FD_MPIO
    printf("Parallel: %d",	 	H5FD_MPIO);
#endif
    fflush(stdout);
    fsync(fileno(stdout));
    // debugging finished
}

//typedef struct {
//    char *name;
//    hid_t h5_file;
//} mr_io_h5;
//
//void mr_io_h5_reader_open(mr_io_h5 *reader_file){
//    // Need to block on file until writing is finished
//    int reader_fd = mr_io_file_reader_open(reader_file->name);
//
//    // Open file (acquires shared lock in non-blocking manner).
//    reader_file->h5_file = H5Fopen(reader_file->name, H5F_ACC_RDONLY, H5P_DEFAULT);
//    mr_io_file_reader_close(reader_fd);
//}
//
//herr_t mr_io_h5_reader_close(mr_io_h5 *reader_file){
//    return H5Fclose(reader_file->h5_file);
//}
//
//void mr_io_h5_writer_open(mr_io_h5 *writer_file){
//    writer_file->h5_file = H5Fcreate(writer_file->name, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
//}
//
//herr_t mr_io_h5_writer_close(mr_io_h5 *writer_file){
////    int *file_handle = NULL;
////    H5Fget_vfd_handle(mr_io_file_out.file, fapl??, &file_handle?? );
////    int fd_out = ???;
//    herr_t status = H5Fflush(writer_file->h5_file, H5F_SCOPE_GLOBAL);  // no fsync?
//    return H5Fclose(writer_file->h5_file);
//}

hid_t mr_io_h5_reader_open(char *name){
    // Need to block on file until writing is finished
    int reader_fd = mr_io_file_reader_open(name);

    // Open file (acquires shared lock in non-blocking manner).
    hid_t h5_file = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);
    mr_io_file_reader_close(reader_fd);
    return h5_file;
}

herr_t mr_io_h5_reader_close(hid_t h5_file){
    return H5Fclose(h5_file);
}

hid_t mr_io_h5_writer_open(char *name){
    hid_t h5_file = H5Fcreate(name, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
//    h5_print_file_driver(H5P_DEFAULT);
    return h5_file;
}

herr_t mr_io_h5_writer_close(hid_t h5_file){
//    int *file_handle = NULL;
//    H5Fget_vfd_handle(h5_file, fapl??, &file_handle?? );
//    int fd_out = ???;
    herr_t status = H5Fflush(h5_file, H5F_SCOPE_GLOBAL);  // no fsync?
    return H5Fclose(h5_file);
}


void mr_io_h5_register_reader( hid_t h5_file,int reader_fd);
int mr_io_h5_deregister_reader(hid_t h5_file);
void mr_io_h5_register_writer( hid_t h5_file, MPI_Comm mr_io_mpi_comm, int writer_fd);
int mr_io_h5_deregister_writer(hid_t h5_file);
MPI_Comm mr_io_h5_get_comm( hid_t h5_file);

hid_t mr_io_h5_parallel_reader_open(MPI_Comm mr_io_mpi_comm, MPI_Info mr_io_mpi_info, char *name){

    // Need to block on file until writing is finished
    int reader_fd = mr_io_file_reader_open(name);

    // Open file (acquires shared lock in non-blocking manner).
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, mr_io_mpi_comm, mr_io_mpi_info);
    hid_t h5_file = H5Fopen(name, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);

    mr_io_h5_register_reader(h5_file, reader_fd);

    return h5_file;
}

herr_t mr_io_h5_parallel_reader_close(hid_t h5_file){
    herr_t status = H5Fclose(h5_file);

    int reader_fd = mr_io_h5_deregister_reader(h5_file);
    mr_io_file_reader_close(reader_fd);

    return status;
}

hid_t mr_io_h5_parallel_writer_open(MPI_Comm mr_io_mpi_comm, MPI_Info mr_io_mpi_info, char *name){

    // Need to create file on rank 0
    int rank = -1;
    // MPI_Comm_rank(mr_io_mpi_comm, &rank);
    if(MPI_SUCCESS != MPI_Comm_rank(mr_io_mpi_comm, &rank)) {
        printf("Failed to get rank on MPI communicator %d (errno=%d, %s)", mr_io_mpi_comm); fflush(stdout);
        exit(-1);
    }

    int writer_fd;
    if (rank == 0) {
        writer_fd = mr_io_file_writer_open(name);
    }

    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, mr_io_mpi_comm, mr_io_mpi_info);
    hid_t h5_file = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

//    h5_print_file_driver(plist_id);
    H5Pclose(plist_id);

    mr_io_h5_register_writer(h5_file, mr_io_mpi_comm, writer_fd);

    return h5_file;
}

herr_t mr_io_h5_parallel_writer_close(hid_t h5_file){
//    int *file_handle = NULL;
//    H5Fget_vfd_handle(h5_file, fapl??, &file_handle?? );
//    int fd_out = ???;
    herr_t status = H5Fflush(h5_file, H5F_SCOPE_GLOBAL);  // no fsync?
    status = H5Fclose(h5_file);

    int rank = -1;
    MPI_Comm mr_io_mpi_comm = mr_io_h5_get_comm(h5_file);
    if(MPI_SUCCESS != MPI_Comm_rank(mr_io_mpi_comm, &rank)) {
        printf("Failed to get rank on MPI communicator %d (errno=%d, %s)", mr_io_mpi_comm); fflush(stdout);
        exit(-1);
    }

    if (rank == 0) {
        int writer_fd = mr_io_h5_deregister_writer(h5_file);
        mr_io_file_writer_close(writer_fd);
    }

    return status;
}

