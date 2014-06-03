! Copyright (C) 2014 RIKEN AICS
! This library is released under the terms of the MIT license.
! http://fiber-miniapp.mit-license.org/

module mod_maprof
  use iso_c_binding
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  private
  ! measuring
  public maprof_time_start, maprof_time_stop
  public maprof_add_fp_ops, maprof_add_ld_ops, maprof_add_st_ops
  public maprof_add_ld_min_ops, maprof_add_st_min_ops
  public maprof_set_ops
  public maprof_set_fp_ops, maprof_set_ld_ops, maprof_set_st_ops
  public maprof_set_ld_min_ops, maprof_set_st_min_ops
  public maprof_get_time, maprof_get_flops
  public maprof_get_throughput, maprof_get_effective_throughput
  public maprof_print, maprof_print_time
#ifdef USE_MPI
  public maprof_print_time_mpi, maprof_print_time_mpi_full
#endif

  ! reporting
  public maprof_setup, maprof_output
  public maprof_app_add_str, maprof_app_add_int, maprof_app_add_float
  public maprof_add_section, maprof_profile_add_problem_size
  public maprof_profile_add_str, maprof_profile_add_int, maprof_profile_add_float

  ! maprof_stat_type
  enum, bind(c)
    enumerator :: MAPROF_ROOT, MAPROF_AVE, MAPROF_MIN, MAPROF_MAX, MAPROF_SD
  end enum

  interface

    subroutine maprof_time_start(id) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
    end subroutine maprof_time_start

    subroutine maprof_time_stop(id) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
    end subroutine maprof_time_stop

    subroutine maprof_add_fp_ops(id, ops) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      real(c_double), value :: ops
    end subroutine maprof_add_fp_ops

    subroutine maprof_add_ld_ops(id, ops) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      real(c_double), value :: ops
    end subroutine maprof_add_ld_ops

    subroutine maprof_add_st_ops(id, ops) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      real(c_double), value :: ops
    end subroutine maprof_add_st_ops

    subroutine maprof_add_ld_min_ops(id, ops) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      real(c_double), value :: ops
    end subroutine maprof_add_ld_min_ops

    subroutine maprof_add_st_min_ops(id, ops) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      real(c_double), value :: ops
    end subroutine maprof_add_st_min_ops

    real(c_double) function maprof_get_time(id, type) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      integer(c_int), value :: type
    end function maprof_get_time

    real(c_double) function maprof_get_flops(id, type) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      integer(c_int), value :: type
    end function maprof_get_flops

    real(c_double) function maprof_get_throughput(id, type) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      integer(c_int), value :: type
    end function maprof_get_throughput

    real(c_double) function maprof_get_effective_throughput(id, type) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      integer(c_int), value :: type
    end function maprof_get_effective_throughput

    subroutine maprof_flush_stdout() bind(c)
      use iso_c_binding
    end subroutine maprof_flush_stdout

  end interface

  interface maprof_set_fp_ops
    module procedure maprof_set_fp_ops_i4, maprof_set_fp_ops_i8, &
         maprof_set_fp_ops_r8
  end interface maprof_set_fp_ops

  interface maprof_set_ld_ops
    module procedure maprof_set_ld_ops_i4, maprof_set_ld_ops_i8, &
         maprof_set_ld_ops_r8
  end interface maprof_set_ld_ops

  interface maprof_set_st_ops
    module procedure maprof_set_st_ops_i4, maprof_set_st_ops_i8, &
         maprof_set_st_ops_r8
  end interface maprof_set_st_ops

  interface maprof_set_ld_min_ops
    module procedure maprof_set_ld_min_ops_i4, &
         maprof_set_ld_min_ops_i8, &
         maprof_set_ld_min_ops_r8
  end interface maprof_set_ld_min_ops

  interface maprof_set_st_min_ops
    module procedure maprof_set_st_min_ops_i4, &
         maprof_set_st_min_ops_i8, &
         maprof_set_st_min_ops_r8
  end interface maprof_set_st_min_ops

  interface maprof_set_ops
    module procedure maprof_set_ops_i4, maprof_set_ops_i8, &
         maprof_set_ops_r8
  end interface maprof_set_ops

  interface

    subroutine maprof_output() bind(c)
      use iso_c_binding
    end subroutine maprof_output

  end interface

contains

function c_string(f_str) result(c_str)
  use iso_c_binding
  character(len=*), intent(in) :: f_str
!!character(len=1, kind=c_char) :: c_str(len_trim(f_str)+1)
  character(len=1, kind=c_char) :: c_str(len(f_str)+1)
  integer :: i, n

!!n = len_trim(f_str)
  n = len(f_str)
  do i = 1, n
    c_str(i) = f_str(i:i)
  end do
  c_str(n+1) = c_null_char

end function c_string


subroutine maprof_set_fp_ops_r8(id, fp, cnt)
  integer, intent(in) :: id
  real(8), value :: fp
  integer, optional, intent(in) :: cnt

  if (present(cnt)) fp = fp * cnt;
  call maprof_add_fp_ops(id, fp)

end subroutine maprof_set_fp_ops_r8


subroutine maprof_set_fp_ops_i8(id, fp, cnt)
  integer, intent(in) :: id
  integer(8), intent(in) :: fp
  integer, optional, intent(in) :: cnt

  call maprof_set_fp_ops_r8(id, dble(fp), cnt)
  
end subroutine maprof_set_fp_ops_i8


subroutine maprof_set_fp_ops_i4(id, fp, cnt)
  integer, intent(in) :: id
  integer, intent(in) :: fp
  integer, optional, intent(in) :: cnt

  call maprof_set_fp_ops_r8(id, dble(fp), cnt)

end subroutine maprof_set_fp_ops_i4


subroutine maprof_set_ld_ops_r8(id, ld, cnt)
  integer, intent(in) :: id
  real(8), value :: ld
  integer, optional, intent(in) :: cnt

  if (present(cnt)) ld = ld * cnt
  call maprof_add_ld_ops(id, ld)

end subroutine maprof_set_ld_ops_r8


subroutine maprof_set_ld_ops_i8(id, ld, cnt)
  integer, intent(in) :: id
  integer(8), intent(in) :: ld
  integer, optional, intent(in) :: cnt

  call maprof_set_ld_ops_r8(id, dble(ld), cnt)

end subroutine maprof_set_ld_ops_i8


subroutine maprof_set_ld_ops_i4(id, ld, cnt)
  integer, intent(in) :: id
  integer, intent(in) :: ld
  integer, optional, intent(in) :: cnt

  call maprof_set_ld_ops_r8(id, dble(ld), cnt)

end subroutine maprof_set_ld_ops_i4


subroutine maprof_set_st_ops_r8(id, st, cnt)
  integer, intent(in) :: id
  real(8), value :: st
  integer, optional, intent(in) :: cnt

  if (present(cnt)) st = st * cnt
  call maprof_add_st_ops(id, st)

end subroutine maprof_set_st_ops_r8


subroutine maprof_set_st_ops_i8(id, st, cnt)
  integer, intent(in) :: id
  integer(8), intent(in) :: st
  integer, optional, intent(in) :: cnt

  call maprof_set_st_ops_r8(id, dble(st), cnt)
  
end subroutine maprof_set_st_ops_i8

subroutine maprof_set_st_ops_i4(id, st, cnt)
  integer, intent(in) :: id
  integer, intent(in) :: st
  integer, optional, intent(in) :: cnt

  call maprof_set_st_ops_r8(id, dble(st), cnt)

end subroutine maprof_set_st_ops_i4


subroutine maprof_set_ld_min_ops_r8(id, ld_min, cnt)
  integer, intent(in) :: id
  real(8), value :: ld_min
  integer, optional, intent(in) :: cnt

  if (present(cnt)) ld_min = ld_min * cnt
  call maprof_add_ld_min_ops(id, ld_min)

end subroutine maprof_set_ld_min_ops_r8


subroutine maprof_set_ld_min_ops_i8(id, ld_min, cnt)
  integer, intent(in) :: id
  integer(8), intent(in) :: ld_min
  integer, optional, intent(in) :: cnt

  call maprof_set_ld_min_ops_r8(id, dble(ld_min), cnt)

end subroutine maprof_set_ld_min_ops_i8


subroutine maprof_set_ld_min_ops_i4(id, ld_min, cnt)
  integer, intent(in) :: id
  integer, intent(in) :: ld_min
  integer, optional, intent(in) :: cnt

  call maprof_set_ld_min_ops_r8(id, dble(ld_min), cnt)

end subroutine maprof_set_ld_min_ops_i4


subroutine maprof_set_st_min_ops_r8(id, st_min, cnt)
  integer, intent(in) :: id
  real(8), value :: st_min
  integer, optional, intent(in) :: cnt

  if (present(cnt)) st_min = st_min * cnt
  call maprof_add_st_min_ops(id, st_min)

end subroutine maprof_set_st_min_ops_r8


subroutine maprof_set_st_min_ops_i8(id, st_min, cnt)
  integer, intent(in) :: id
  integer(8), intent(in) :: st_min
  integer, optional, intent(in) :: cnt

  call maprof_set_st_min_ops_r8(id, dble(st_min), cnt)

end subroutine maprof_set_st_min_ops_i8


subroutine maprof_set_st_min_ops_i4(id, st_min, cnt)
  integer, intent(in) :: id
  integer, intent(in) :: st_min
  integer, optional, intent(in) :: cnt

  call maprof_set_st_min_ops_r8(id, dble(st_min), cnt)

end subroutine maprof_set_st_min_ops_i4


subroutine maprof_set_ops_r8(sec, fp, ld, st, ld_min, st_min, cnt)
  integer, intent(in) :: sec
  real(8), intent(in) :: fp, ld, st
  real(8), intent(in) :: ld_min, st_min
  integer, optional, intent(in) :: cnt

  call maprof_set_fp_ops_r8(sec, fp, cnt)
  call maprof_set_ld_ops_r8(sec, ld, cnt)
  call maprof_set_st_ops_r8(sec, st, cnt)
  call maprof_set_ld_min_ops_r8(sec, ld_min, cnt)
  call maprof_set_st_min_ops_r8(sec, st_min, cnt)

end subroutine maprof_set_ops_r8


subroutine maprof_set_ops_i8(sec, fp, ld, st, ld_min, st_min, cnt)
  integer, intent(in) :: sec
  integer(8), intent(in) :: fp, ld, st
  integer(8), intent(in) :: ld_min, st_min
  integer, optional, intent(in) :: cnt

  call maprof_set_ops_r8(sec, dble(fp), dble(ld), dble(st), &
       dble(ld_min), dble(st_min), cnt)
  
end subroutine maprof_set_ops_i8


subroutine maprof_set_ops_i4(sec, fp, ld, st, ld_min, st_min, cnt)
  integer, intent(in) :: sec
  integer, intent(in) :: fp, ld, st
  integer, intent(in) :: ld_min, st_min
  integer, optional, intent(in) :: cnt

  call maprof_set_ops_r8(sec, dble(fp), dble(ld), dble(st), &
                      dble(ld_min), dble(st_min), cnt)

end subroutine maprof_set_ops_i4


subroutine maprof_print(id, name)
  integer, intent(in) :: id
  character(*), intent(in) :: name

  interface
    subroutine c_print(id, str) bind(c, name="maprof_print")
      use iso_c_binding
      integer(c_int), value :: id 
      character(kind=c_char), intent(in) :: str(*)
    end subroutine c_print
  end interface

  call c_print(id, c_string(name))
  call maprof_flush_stdout()

end subroutine maprof_print


subroutine maprof_print_time(id, name)
  integer, intent(in) :: id
  character(*), intent(in) :: name

  interface
    subroutine c_print_time(id, str) bind(c, name="maprof_print_time")
      use iso_c_binding
      integer(c_int), value :: id 
      character(kind=c_char), intent(in) :: str(*)
    end subroutine c_print_time
  end interface

  call c_print_time(id, c_string(name))
  call maprof_flush_stdout()

end subroutine maprof_print_time


#ifdef USE_MPI
subroutine maprof_print_time_mpi(id, name)
  integer, intent(in) :: id
  character(*), intent(in) :: name

  interface
    subroutine c_print_time_mpi(id, str) bind(c, name="maprof_print_time_mpi")
      use iso_c_binding
      integer(c_int), value :: id 
      character(kind=c_char), intent(in) :: str(*)
    end subroutine c_print_time_mpi
  end interface

  call c_print_time_mpi(id, c_string(name))
  call maprof_flush_stdout()

end subroutine maprof_print_time_mpi


subroutine maprof_print_time_mpi_full(id, name)
  integer, intent(in) :: id
  character(*), intent(in) :: name

  interface
    subroutine c_print_time_mpi_full(id, str) &
                       bind(c, name="maprof_print_time_mpi_full")
      use iso_c_binding
      integer(c_int), value :: id 
      character(kind=c_char), intent(in) :: str(*)
    end subroutine c_print_time_mpi_full
  end interface

  call c_print_time_mpi_full(id, c_string(name))
  call maprof_flush_stdout()

end subroutine maprof_print_time_mpi_full

#endif

subroutine maprof_setup(app_name, app_version)
!$ use omp_lib
  character(*), intent(in) :: app_name, app_version

  integer :: nt = 1

  interface
    subroutine c_setup(app_name, app_version) bind(c, name="maprof_setup")
      use iso_c_binding
      character(kind=c_char), intent(in) :: app_name(*), app_version(*)
    end subroutine c_setup

    subroutine maprof_set_num_threads(n) bind(c)
      use iso_c_binding
      integer(c_int), value :: n
    end subroutine maprof_set_num_threads
  end interface

!$ nt = omp_get_max_threads()
  call maprof_set_num_threads(nt)

  call c_setup(c_string(app_name), c_string(app_version))

end subroutine maprof_setup


subroutine maprof_app_add_str(key, str)
  character(*), intent(in) :: key, str

  interface
    subroutine c_app_add_str(key, str) bind(c, name="maprof_app_add_str")
      use iso_c_binding
      character(kind=c_char), intent(in) :: key(*), str(*)
    end subroutine c_app_add_str
  end interface

  call c_app_add_str(c_string(key), c_string(str))

end subroutine maprof_app_add_str


subroutine maprof_app_add_int(key, n)
  character(*), intent(in) :: key
  integer, intent(in) :: n

  interface
    subroutine c_app_add_int(key, n) bind(c, name="maprof_app_add_int")
      use iso_c_binding
      character(kind=c_char), intent(in) :: key(*)
      integer(c_int), value :: n
    end subroutine c_app_add_int
  end interface

  call c_app_add_int(c_string(key), n)

end subroutine maprof_app_add_int


subroutine maprof_app_add_float(key, r)
  character(*), intent(in) :: key
  real(8), intent(in) :: r

  interface
    subroutine c_app_add_float(key, r) bind(c, name="maprof_app_add_float")
      use iso_c_binding
      character(kind=c_char), intent(in) :: key(*)
      real(c_double), value :: r
    end subroutine c_app_add_float
  end interface

  call c_app_add_float(c_string(key), r)

end subroutine maprof_app_add_float


subroutine maprof_add_section(name, id)
  character(*), intent(in) :: name
  integer, intent(in) :: id

  interface
    subroutine c_add_section(name, id) bind(c, name="maprof_add_section")
      use iso_c_binding
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int), value :: id
    end subroutine c_add_section
  end interface

  call c_add_section(c_string(name), id)

end subroutine maprof_add_section


subroutine maprof_profile_add_problem_size(key, n)
  character(*), intent(in) :: key
  integer, intent(in) :: n

  interface
    subroutine c_profile_add_problem_size(key, n)  &
                   bind(c, name="maprof_profile_add_problem_size")
      use iso_c_binding
      character(kind=c_char), intent(in) :: key(*)
      integer(c_int), value :: n
    end subroutine c_profile_add_problem_size
  end interface

  call c_profile_add_problem_size(c_string(key), n)

end subroutine maprof_profile_add_problem_size


subroutine maprof_profile_add_str(key, str)
  character(*), intent(in) :: key, str

  interface
    subroutine c_profile_add_str(key, str) bind(c, name="maprof_profile_add_str")
      use iso_c_binding
      character(kind=c_char), intent(in) :: key(*), str(*)
    end subroutine c_profile_add_str
  end interface

  call c_profile_add_str(c_string(key), c_string(str))

end subroutine maprof_profile_add_str


subroutine maprof_profile_add_int(key, n)
  character(*), intent(in) :: key
  integer, intent(in) :: n

  interface
    subroutine c_profile_add_int(key, n) bind(c, name="maprof_profile_add_int")
      use iso_c_binding
      character(kind=c_char), intent(in) :: key(*)
      integer(c_int), value :: n
    end subroutine c_profile_add_int
  end interface

  call c_profile_add_int(c_string(key), n)

end subroutine maprof_profile_add_int


subroutine maprof_profile_add_float(key, r)
  character(*), intent(in) :: key
  real(8), intent(in) :: r

  interface
    subroutine c_profile_add_float(key, r)   &
                      bind(c, name="maprof_profile_add_float")
      use iso_c_binding
      character(kind=c_char), intent(in) :: key(*)
      real(c_double), value :: r
    end subroutine c_profile_add_float
  end interface

  call c_profile_add_float(c_string(key), r)

end subroutine maprof_profile_add_float

end module mod_maprof
