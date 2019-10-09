module mr_io_test_arg_parser

    ! parallel tests
    integer, dimension(3) :: mr_io_test_mpi_cart_dims
    ! reader tests
    character(len=200) :: path
    ! reader writer tests
    character(len=200) :: in_path
    character(len=200) :: out_path
    ! padded tests
    integer, dimension(3) :: domain_padding_lhs
    integer, dimension(3) :: domain_padding_rhs

    integer :: simulation_time_refinement
    integer, dimension(3) :: simulation_spatial_refinement

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

    mr_io_test_mpi_cart_dims = (/NB1, NB2, NB3/)

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

    implicit none

    character(len=100) :: pad_lhs_1_c
    character(len=100) :: pad_lhs_2_c
    character(len=100) :: pad_lhs_3_c

    character(len=100) :: pad_rhs_1_c
    character(len=100) :: pad_rhs_2_c
    character(len=100) :: pad_rhs_3_c

    integer :: pad_lhs_1
    integer :: pad_lhs_2
    integer :: pad_lhs_3

    integer :: pad_rhs_1
    integer :: pad_rhs_2
    integer :: pad_rhs_3

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

    implicit none

    character(len=100) :: pad_lhs_1_c
    character(len=100) :: pad_lhs_2_c
    character(len=100) :: pad_lhs_3_c

    character(len=100) :: pad_rhs_1_c
    character(len=100) :: pad_rhs_2_c
    character(len=100) :: pad_rhs_3_c

    integer :: pad_lhs_1
    integer :: pad_lhs_2
    integer :: pad_lhs_3

    integer :: pad_rhs_1
    integer :: pad_rhs_2
    integer :: pad_rhs_3

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


subroutine mr_io_test_parse_args_parallel_reader_writer_padded_to_st

    implicit none

    character(len=100) :: t_refinement_c
    character(len=100) :: x_refinement_c
    character(len=100) :: y_refinement_c
    character(len=100) :: z_refinement_c

    integer :: t_refinement
    integer :: x_refinement
    integer :: y_refinement
    integer :: z_refinement

    call mr_io_test_parse_args_parallel_reader_writer_padded()

    call get_command_argument(12, t_refinement_c)
    call get_command_argument(13, x_refinement_c)
    call get_command_argument(14, y_refinement_c)
    call get_command_argument(15, z_refinement_c)

    read(t_refinement_c,*)t_refinement
    read(x_refinement_c,*)x_refinement
    read(y_refinement_c,*)y_refinement
    read(z_refinement_c,*)z_refinement

    simulation_time_refinement = t_refinement
    simulation_spatial_refinement = (/ x_refinement, y_refinement, z_refinement /)

end subroutine mr_io_test_parse_args_parallel_reader_writer_padded_to_st

end module mr_io_test_arg_parser
