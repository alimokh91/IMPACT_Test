#include <set>
#include <map>
#include "hdf5.h"

extern "C" {
    void mr_io_h5_register_reader( hid_t h5_file,int reader_fd);
    int mr_io_h5_deregister_reader(hid_t h5_file);
    void mr_io_h5_register_writer( hid_t h5_file, MPI_Comm mr_io_mpi_comm, int writer_fd);
    int mr_io_h5_deregister_writer(hid_t h5_file);
    MPI_Comm mr_io_h5_get_comm( hid_t h5_file);
}


struct writer_fd_t {
    MPI_Comm comm;
    int fd;
};

static std::map<hid_t, writer_fd_t> h5_to_fd_writer;
static std::map<hid_t, int> h5_to_fd_reader;

void mr_io_h5_register_writer(hid_t h5_file, MPI_Comm mr_io_mpi_comm, int writer_fd) {
    struct writer_fd_t writer_fd_val = { .comm = mr_io_mpi_comm, .fd = writer_fd };
    h5_to_fd_writer[h5_file] = writer_fd_val;
}

MPI_Comm mr_io_h5_get_comm(hid_t h5_file) {
    return h5_to_fd_writer.at(h5_file).comm;
}

int mr_io_h5_deregister_writer(hid_t h5_file){
    int writer_fd = h5_to_fd_writer.at(h5_file).fd;
    h5_to_fd_writer.erase(h5_file);
    return writer_fd;
}

void mr_io_h5_register_reader(hid_t h5_file, int reader_fd) {
    h5_to_fd_reader[h5_file] = reader_fd;
}

int mr_io_h5_deregister_reader(hid_t h5_file){
    int reader_fd = h5_to_fd_reader.at(h5_file);
    h5_to_fd_reader.erase(h5_file);
    return reader_fd;
}
