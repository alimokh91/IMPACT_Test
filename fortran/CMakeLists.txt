set(HPCPREDICT_IO_SOURCES
  mr_io.f90
  mr_io_parallel.f90
  mr_io_parallel_spacetime.f90
  mr_protocol.f90
)

find_package(HDF5 COMPONENTS Fortran REQUIRED)

add_library(hpc-predict-io ${HPCPREDICT_IO_SOURCES})
target_compile_options(hpc-predict-io PUBLIC -cpp)
target_include_directories(hpc-predict-io PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(hpc-predict-io PUBLIC ${HDF5_LIBRARIES})


foreach(file mr_io_example_space_time mr_io_example_spatial
    mr_io_test_parallel_reader 
    mr_io_test_parallel_reader_space_time
    mr_io_test_parallel_reader_hpc_predict
    mr_io_test_parallel_reader_writer
    mr_io_test_parallel_reader_writer_space_time 
    mr_io_test_parallel_reader_writer_hpc_predict
    mr_io_test_reader
    mr_io_test_reader_space_time 
    mr_io_test_reader_hpc_predict
    mr_io_test_reader_writer
    mr_io_test_reader_writer_space_time
    mr_io_test_reader_writer_hpc_predict
)
  add_executable(${file} ${file}.f90)
  target_link_libraries(${file} hpc-predict-io)
endforeach()
