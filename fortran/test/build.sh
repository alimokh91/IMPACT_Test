#!/bin/bash

set -euxo pipefail

export IMPACT_DIR=/home/lukasd/src/IMPACT/prog
export HDF5_DIR=/home/lukasd/src/build_whale/opt_scientific/libs/petsc
export MPI_DIR=/usr/include/mpich

mpifort -ffree-line-length-0 -g -c -fbacktrace mr_io_test_arg_parser.f90 -I.. -I${HDF5_DIR}/include
mpifort -ffree-line-length-0 -g -c -cpp -fbacktrace mr_io_test_impact_input.f90 -I.. -I${IMPACT_DIR} -I${MPI_DIR} -I${HDF5_DIR}/include
mpifort -ffree-line-length-0 -g -c -cpp -fbacktrace mr_io_test_impact_mri.f90 -I.. -I${IMPACT_DIR} -I${MPI_DIR} -I${HDF5_DIR}/include
mpifort -ffree-line-length-0 -g -c -cpp -fbacktrace  *.f90 -I.. -I${IMPACT_DIR}  -I${MPI_DIR}  -I${HDF5_DIR}/include

mpifort mr_io_example_spatial.o ../mr_io.o ../mr_io_parallel.o ../mr_io_protocol.o -o mr_io_example_spatial -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_example_space_time.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_protocol.o -o mr_io_example_space_time -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran

mpifort mr_io_test_reader.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_reader -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_reader_space_time.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_reader_space_time -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_reader_flow.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_reader_flow -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_reader_segmented_flow.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_reader_segmented_flow -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran

mpifort mr_io_test_parallel_reader.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_parallel_reader -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_parallel_reader_space_time.o ../mr_io.o ../mr_io_parallel.o ../mr_io_parallel_spacetime.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_parallel_reader_space_time -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_parallel_reader_flow.o ../mr_io.o ../mr_io_parallel.o ../mr_io_parallel_spacetime.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_parallel_reader_flow -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_parallel_reader_segmented_flow.o ../mr_io.o ../mr_io_parallel.o ../mr_io_parallel_spacetime.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_parallel_reader_segmented_flow -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_parallel_reader_flow_padded.o ../mr_io.o ../mr_io_parallel.o ../mr_io_parallel_spacetime.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_parallel_reader_flow_padded -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_parallel_reader_segmented_flow_padded.o ../mr_io.o ../mr_io_parallel.o ../mr_io_parallel_spacetime.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_parallel_reader_segmented_flow_padded -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran

mpifort mr_io_test_reader_writer.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_reader_writer -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_reader_writer_space_time.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_reader_writer_space_time -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_reader_writer_flow.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_reader_writer_flow -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_reader_writer_segmented_flow.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_reader_writer_segmented_flow -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran

mpifort mr_io_test_parallel_reader_writer.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_parallel_reader_writer -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_parallel_reader_writer_space_time.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_parallel_spacetime.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_parallel_reader_writer_space_time -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_parallel_reader_writer_flow.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_parallel_spacetime.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_parallel_reader_writer_flow -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_parallel_reader_writer_segmented_flow.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_parallel_spacetime.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_parallel_reader_writer_segmented_flow -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_parallel_reader_writer_flow_padded.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_parallel_spacetime.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_parallel_reader_writer_flow_padded -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_parallel_reader_writer_flow_padded_to_space_time.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_parallel_spacetime.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_parallel_reader_writer_flow_padded_to_space_time -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_parallel_reader_writer_segmented_flow_padded.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_parallel_spacetime.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_parallel_reader_writer_segmented_flow_padded -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran
mpifort mr_io_test_parallel_reader_writer_segmented_flow_padded_to_space_time.o ../mr_io.o ../mr_io_parallel.o  ../mr_io_parallel_spacetime.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_parallel_reader_writer_segmented_flow_padded_to_space_time -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -lhdf5_fortran

mpifort mr_io_test_impact_input.o ../mr_io.o ../mr_io_parallel.o ../mr_io_parallel_spacetime.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_impact_input -L${IMPACT_DIR} -Wl,-rpath=${IMPACT_DIR} -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -limpact -lm -ldl -lz -lblas -llapack -lhdf5_fortran
mpifort mr_io_test_impact_mri.o ../mr_io.o ../mr_io_parallel.o ../mr_io_parallel_spacetime.o  ../mr_io_protocol.o mr_io_test_arg_parser.o -o mr_io_test_impact_mri -L${IMPACT_DIR} -Wl,-rpath=${IMPACT_DIR} -L${HDF5_DIR}/lib -Wl,-rpath=${HDF5_DIR}/lib -limpact -lm -ldl -lz -lblas -llapack -lhdf5_fortran
