program main
#ifdef USE_MPI
  use mpi
#endif
  use mod_maprof
  implicit none

  real(8):: x, y
  integer:: n = 10000000

  integer, parameter :: ALL  = 0
  integer, parameter :: SEC1 = 1
  integer, parameter :: SEC2 = 2

#ifdef USE_MPI
  integer:: ierr, myrank
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
#endif

  call maprof_time_start(ALL)

  call maprof_time_start(SEC1)
  call sub1(x, n)
  call maprof_time_stop(SEC1)

  call maprof_time_start(SEC2)
  call sub2(y, n)
  call maprof_time_stop(SEC2)

  write(*,*) 'x+y = ', x + y

  call maprof_time_stop(ALL)

  ! operation counts are integer(4)
  call maprof_set_ops(SEC1, 10, 20, 30, 40, 50, n)

  ! operation counts are integer(8)
 !call maprof_set_ops(SEC1, 10_8, 20_8, 30_8, 40_8, 50_8, n)

  ! set operation counts separately
  call maprof_set_fp_ops(SEC2, 11, n)
  call maprof_set_ld_ops(SEC2, 22_8, n)
  call maprof_set_st_ops(SEC2, 33*n)  ! 3rd argument is optional(defalut=1)
  call maprof_set_ld_min_ops(SEC2, 44_8*n)
  call maprof_set_st_min_ops(SEC2, 55, n)

call maprof_print(SEC1, 'sec1')
call maprof_print(SEC2, 'sec2')
#ifdef USE_MPI
  call maprof_print_time_mpi(SEC1, 'sec1:')  ! all processes must call this
  call maprof_print_time_mpi(SEC2, 'sec2:')  ! all processes must call this
  call maprof_print_time_mpi_full(ALL, 'all:') 
#else
  call maprof_print_time(SEC1, 'sec1 time =')
  call maprof_print_time(SEC2, 'sec2 time =')
#endif

  call maprof_setup("f_test", "1.0.0")
  call maprof_app_add_str("note", "test for fortran programs")
  call maprof_add_section("sec1", SEC1)
  call maprof_add_section("sec2", SEC2)
  call maprof_profile_add_problem_size("n", n)
  call maprof_output()



#ifdef USE_MPI
  call MPI_Finalize(ierr)
#endif

contains


subroutine sub1(x, n)
  real(8):: x
  integer:: n

  integer:: i

  x = 0.0
  do i = 1, n
    x = x + exp(-sqrt(i*0.1)**(3.14))
  end do
end subroutine sub1


subroutine sub2(x, n)
  real(8):: x
  integer:: n

  call sub1(x, n)

end subroutine sub2


end program main

