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
target_compile_options(hpc-predict-io PUBLIC -cpp -g -fbacktrace)
target_include_directories(hpc-predict-io PUBLIC ${HDF5_INCLUDE_DIRS}
    INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/fmod)
target_link_libraries(hpc-predict-io PUBLIC ${HDF5_LIBRARIES})
set_target_properties(hpc-predict-io PROPERTIES Fortran_MODULE_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}/fmod)
install(TARGETS hpc-predict-io
        ARCHIVE DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fmod/ DESTINATION include)


foreach(file 
    mr_io_example_space_time 
    mr_io_example_spatial
)
  add_executable(${file} ${file}.f90)
  target_compile_options(${file} PUBLIC -cpp -g -fbacktrace)
  target_link_libraries(${file} hpc-predict-io)
  install(TARGETS ${file}
          RUNTIME DESTINATION bin)
endforeach()


foreach(file 
    mr_io_test_reader
    mr_io_test_reader_space_time
    mr_io_test_reader_hpc_predict
    mr_io_test_reader_segmented_hpc_predict
    mr_io_test_reader_writer
    mr_io_test_reader_writer_space_time
    mr_io_test_reader_writer_hpc_predict
    mr_io_test_reader_writer_segmented_hpc_predict
    mr_io_test_parallel_reader
    mr_io_test_parallel_reader_space_time
    mr_io_test_parallel_reader_hpc_predict
    mr_io_test_parallel_reader_hpc_predict_padded
    mr_io_test_parallel_reader_segmented_hpc_predict
    mr_io_test_parallel_reader_writer
    mr_io_test_parallel_reader_writer_space_time
    mr_io_test_parallel_reader_writer_hpc_predict
    mr_io_test_parallel_reader_writer_hpc_predict_padded
    mr_io_test_parallel_reader_writer_segmented_hpc_predict
)
  add_executable(${file} ${file}.f90 mr_io_test_arg_parser.f90)
  target_compile_options(${file} PUBLIC -cpp -g -fbacktrace)
  target_link_libraries(${file} hpc-predict-io)
  install(TARGETS ${file}
          RUNTIME DESTINATION bin)
endforeach()


# FIXME: Do some refactoring on these tests
set( IMPACT_LIBRARIES libimpact.a )
find_file(IMPACT_LIB ${IMPACT_LIBRARIES} PATHS ENV IMPACT_DIR NO_DEFAULT_PATH)

if (EXISTS $ENV{IMPACT_DIR})
    message("Found IMPACT... building IMPACT startup tests.")
	
    set( IMPACT_DEPS libm.so libdl.so libz.so libblas.so liblapack.so )
    
    add_executable(mr_io_test_impact_input mr_io_test_impact_input.f90)
	target_link_libraries(mr_io_test_impact_input PRIVATE ${IMPACT_LIBRARIES} ${IMPACT_DEPS} ${HDF5_LIBRARIES} hpc-predict-io)
	target_link_libraries(mr_io_test_impact_input PRIVATE -L$ENV{IMPACT_DIR} )
	target_include_directories(mr_io_test_impact_input PRIVATE $ENV{IMPACT_DIR} )
	target_link_libraries(mr_io_test_impact_input PRIVATE -L/home/lukasd/src/build_whale/opt_scientific/libs/petsc/lib )
	target_include_directories(mr_io_test_impact_input PRIVATE /home/lukasd/src/build_whale/opt_scientific/libs/petsc/include )
	
	install(TARGETS mr_io_test_impact_input
	        RUNTIME DESTINATION bin)
else ()
    message("Could not find IMPACT... skipping IMPACT startup tests.")
endif ()
        
