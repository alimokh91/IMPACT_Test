find_package(HDF5 COMPONENTS Fortran REQUIRED)


foreach(file 
    mr_io_example_space_time 
    mr_io_example_spatial
)
  add_executable(${file} ${file}.f90)
  target_compile_options(${file} PUBLIC -cpp -g -fbacktrace)
  target_link_libraries(${file} hpc-predict-io)
  set_target_properties(${file} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
  install(TARGETS ${file}
          RUNTIME DESTINATION bin/examples)
endforeach()


foreach(file 
    mr_io_test_reader
    mr_io_test_reader_space_time
    mr_io_test_reader_flow
    mr_io_test_reader_segmented_flow
    mr_io_test_reader_writer
    mr_io_test_reader_writer_space_time
    mr_io_test_reader_writer_flow
    mr_io_test_reader_writer_segmented_flow
    mr_io_test_parallel_reader
    mr_io_test_parallel_reader_space_time
    mr_io_test_parallel_reader_flow
    mr_io_test_parallel_reader_flow_padded
    mr_io_test_parallel_reader_segmented_flow
    mr_io_test_parallel_reader_writer
    mr_io_test_parallel_reader_writer_space_time
    mr_io_test_parallel_reader_writer_flow
    mr_io_test_parallel_reader_writer_flow_padded
    mr_io_test_parallel_reader_writer_flow_padded_to_space_time
    mr_io_test_parallel_reader_writer_segmented_flow
)
  add_executable(${file} ${file}.f90 mr_io_test_arg_parser.f90)
  target_compile_options(${file} PUBLIC -cpp -g -fbacktrace -ffree-line-length-0)
  target_link_libraries(${file} hpc-predict-io)
  set_target_properties(${file} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
  install(TARGETS ${file}
          RUNTIME DESTINATION bin/test)
endforeach()


# FIXME: Do some refactoring on these tests
set( IMPACT_LIBRARIES libimpact.a )
find_file(IMPACT_LIB ${IMPACT_LIBRARIES} PATHS ENV IMPACT_DIR NO_DEFAULT_PATH)

if (IMPACT_LIB)
    message("Found IMPACT... building IMPACT startup tests.")
	
    set( IMPACT_DEPS libm.so libdl.so libz.so libblas.so liblapack.so )
    
    foreach(file 
        mr_io_test_impact_input
        mr_io_test_impact_mri
        )
	    add_executable(${file} ${file}.f90)
        target_compile_options(${file} PUBLIC -cpp -g -fbacktrace -ffree-line-length-0)
	    target_link_libraries(${file} PRIVATE ${IMPACT_LIBRARIES} ${IMPACT_DEPS} ${HDF5_LIBRARIES} hpc-predict-io)
	    target_link_libraries(${file} PRIVATE -L$ENV{IMPACT_DIR} )
	    target_include_directories(${file} PRIVATE $ENV{IMPACT_DIR} ${HDF5_INCLUDES} )
	    set_target_properties(${file} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)	
	    install(TARGETS ${file}
		        RUNTIME DESTINATION bin/test)
	endforeach()
else ()
    message("Could not find IMPACT... skipping IMPACT startup tests.")
endif ()
        
