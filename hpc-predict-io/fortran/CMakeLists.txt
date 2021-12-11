set(HPCPREDICT_IO_SOURCES
  mr_io_protocol.f90
  mr_io_hyperslab_enum.f90
  mr_io.f90
  mr_io_mpi_error_handler.f90
  mr_io_parallel.f90
  mr_io_parallel_spacetime.f90
)

find_package(HDF5 COMPONENTS Fortran REQUIRED)

add_library(hpc-predict-io ${HPCPREDICT_IO_SOURCES})
target_compile_options(hpc-predict-io PUBLIC -cpp)
target_include_directories(hpc-predict-io PUBLIC ${HDF5_INCLUDE_DIRS}
    INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/fmod)
target_link_libraries(hpc-predict-io PUBLIC ${HDF5_LIBRARIES})
set_target_properties(hpc-predict-io PROPERTIES Fortran_MODULE_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}/fmod)
install(TARGETS hpc-predict-io
        ARCHIVE DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fmod/ DESTINATION include)

add_subdirectory(test)
        
