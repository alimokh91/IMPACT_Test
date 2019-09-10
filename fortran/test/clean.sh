#!/bin/bash

set -uxo pipefail

rm  *.mod \
    *.o \
    mr_io_example_spatial \
    mr_io_test_reader \
    mr_io_test_reader_space_time \
    mr_io_test_reader_flow \
    mr_io_test_parallel_reader \
    mr_io_test_parallel_reader_space_time \
    mr_io_test_parallel_reader_flow \
    mr_io_test_parallel_reader_spacetime \
    mr_io_test_reader_writer \
    mr_io_test_parallel_reader_writer \
    mr_io_test_parallel_reader_writer_space_time \
    mr_io_test_parallel_reader_writer_flow \
    mr_io_test_reader_writer_space_time \
    mr_io_test_reader_writer_flow \
    mr_io_example_space_time \
    mr_io_test_reader_segmented_flow \
    mr_io_test_reader_writer_segmented_flow \
    mr_io_test_parallel_reader_segmented_flow \
    mr_io_test_parallel_reader_flow_padded \
    mr_io_test_parallel_reader_flow_with_padding \
    mr_io_test_parallel_reader_writer_segmented_flow \
    mr_io_test_parallel_reader_writer_flow_padded \
    mr_io_test_impact_input \
