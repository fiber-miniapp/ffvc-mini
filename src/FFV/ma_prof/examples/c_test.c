#include <stdio.h>
#include <math.h>
#include "maprof.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

enum {
  ALL,
  FUNC1,
  FUNC2,
};

double func1(int n);
double func2(int n);

int main(int argc, char **argv)
{
  int n = 1000000;
  double a, b;

#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif

  maprof_time_start(ALL);

  maprof_time_start(FUNC1);
  a = func1(n);
  maprof_time_stop(FUNC1);

  maprof_time_start(FUNC2);
  b = func2(n);
  maprof_time_stop(FUNC2);

  printf("a + b = %f\n", a + b);

  maprof_time_stop(ALL);

  maprof_print(ALL, "all:");
#ifdef USE_MPI
  maprof_print_time_mpi(FUNC1, "func1: ");
  maprof_print_time_mpi(FUNC2, "func2: ");
  maprof_print_time_mpi_full(ALL, "all: ");
#else
  maprof_print_time(FUNC1, "func1: ");
  maprof_print_time(FUNC2, "func2: ");
  maprof_print_time(ALL, "all:   ");
#endif

  maprof_setup("c_test", "1.0.0");
  maprof_add_section("fucn1", FUNC1);
  maprof_add_section("fucn2", FUNC2);
  maprof_profile_add_problem_size("n", n);
  maprof_output();

#ifdef USE_MPI
  MPI_Finalize();
#endif

  return 0;
}


double func1(int n)
{
  int i;
  double x;
  x = 0.0;
  for (i = 0; i < n; i++) {
    x += exp(-pow(sqrt(fabs(sin(i*3.14))), 1.234));
  }
  return x;
}


double func2(int n)
{
  return func1(n);
}

