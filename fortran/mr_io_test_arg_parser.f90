module mr_io_test_arg_parser

    ! parallel tests
    integer, dimension(3) :: mpi_cart_dims
    ! reader tests
    character(len=200) :: path
    ! reader writer tests
    character(len=200) :: in_path
    character(len=200) :: out_path
    ! padded tests
    integer, dimension(3) :: domain_padding_lhs
    integer, dimension(3) :: domain_padding_rhs

contains

! Sequential tests

subroutine mr_io_test_parse_args_reader()

    implicit none

    call get_command_argument(1,path)

    !    write(*,*) path
    !    call flush()

end subroutine mr_io_test_parse_args_reader

subroutine mr_io_test_parse_args_reader_writer()

    implicit none

    call get_command_argument(1,in_path)
    call get_command_argument(2,out_path)

    !    write(*,*) in_path
    !    write(*,*) out_path
    !    call flush()

end subroutine mr_io_test_parse_args_reader_writer

! Parallel tests

! MPI topology
subroutine mr_io_test_parse_args_parallel_topology()

    implicit none

    character(len=100) :: NB1char
    character(len=100) :: NB2char
    character(len=100) :: NB3char
    real :: NB1
    real :: NB2
    real :: NB3

    ! Read domain decomposition
    call get_command_argument(1,NB1char)
    call get_command_argument(2,NB2char)
    call get_command_argument(3,NB3char)

    read(NB1char,*)NB1
    read(NB2char,*)NB2
    read(NB3char,*)NB3

    mpi_cart_dims = (/NB1, NB2, NB3/)

    call get_command_argument(4,path)

    !    write(*,*) NB1char
    !    write(*,*) NB2char
    !    write(*,*) NB3char
    !    call flush()

end subroutine mr_io_test_parse_args_parallel_topology

subroutine mr_io_test_parse_args_parallel_reader()

    implicit none

    call mr_io_test_parse_args_parallel_topology()

    call get_command_argument(4,path)

    !    write(*,*) path
    !    call flush()

end subroutine mr_io_test_parse_args_parallel_reader

subroutine mr_io_test_parse_args_parallel_reader_writer()

    implicit none

    call mr_io_test_parse_args_parallel_topology()

    call get_command_argument(4,in_path)
    call get_command_argument(5,out_path)

    !    write(*,*) in_path
    !    write(*,*) out_path
    !    call flush()

end subroutine mr_io_test_parse_args_parallel_reader_writer

subroutine mr_io_test_parse_args_parallel_reader_padded

    character(len=100) :: pad_lhs_1_c
    character(len=100) :: pad_lhs_2_c
    character(len=100) :: pad_lhs_3_c

    character(len=100) :: pad_rhs_1_c
    character(len=100) :: pad_rhs_2_c
    character(len=100) :: pad_rhs_3_c

    call mr_io_test_parse_args_parallel_reader()

    call get_command_argument(5,pad_lhs_1_c)
    call get_command_argument(6,pad_lhs_2_c)
    call get_command_argument(7,pad_lhs_3_c)

    call get_command_argument(8, pad_rhs_1_c)
    call get_command_argument(9, pad_rhs_2_c)
    call get_command_argument(10,pad_rhs_3_c)

    read(pad_lhs_1_c,*)pad_lhs_1
    read(pad_lhs_2_c,*)pad_lhs_2
    read(pad_lhs_3_c,*)pad_lhs_3

    read(pad_rhs_1_c,*)pad_rhs_1
    read(pad_rhs_2_c,*)pad_rhs_2
    read(pad_rhs_3_c,*)pad_rhs_3

    domain_padding_lhs = (/pad_lhs_1, pad_lhs_2, pad_lhs_3/)
    domain_padding_rhs = (/pad_rhs_1, pad_rhs_2, pad_rhs_3/)

end subroutine mr_io_test_parse_args_parallel_reader_padded

subroutine mr_io_test_parse_args_parallel_reader_writer_padded

    character(len=100) :: pad_lhs_1_c
    character(len=100) :: pad_lhs_2_c
    character(len=100) :: pad_lhs_3_c

    character(len=100) :: pad_rhs_1_c
    character(len=100) :: pad_rhs_2_c
    character(len=100) :: pad_rhs_3_c

    call mr_io_test_parse_args_parallel_reader_writer()

    call get_command_argument(6,pad_lhs_1_c)
    call get_command_argument(7,pad_lhs_2_c)
    call get_command_argument(8,pad_lhs_3_c)

    call get_command_argument(9, pad_rhs_1_c)
    call get_command_argument(10,pad_rhs_2_c)
    call get_command_argument(11,pad_rhs_3_c)

    read(pad_lhs_1_c,*)pad_lhs_1
    read(pad_lhs_2_c,*)pad_lhs_2
    read(pad_lhs_3_c,*)pad_lhs_3

    read(pad_rhs_1_c,*)pad_rhs_1
    read(pad_rhs_2_c,*)pad_rhs_2
    read(pad_rhs_3_c,*)pad_rhs_3

    domain_padding_lhs = (/pad_lhs_1, pad_lhs_2, pad_lhs_3/)
    domain_padding_rhs = (/pad_rhs_1, pad_rhs_2, pad_rhs_3/)

end subroutine mr_io_test_parse_args_parallel_reader_writer_padded


end module mr_io_test_arg_parser
