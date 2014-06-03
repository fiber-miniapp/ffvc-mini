//
// Copyright (c) 2012-2014 Institute of Industrial Science, The University of Tokyo.
// All right reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All right reserved.
//

#ifndef _TIMING_H_
#define _TIMING_H_

enum {
  tm_TOTAL,
  tm_INIT,
  tm_MAIN,
    tm_SETUP,
    tm_NS_FS_E_BINARY,
      tm_PVEC_MUSCL,
      tm_VP_ITR,
        tm_SOR_2_SMA,
          tm_PSOR2SMA_CORE,
          tm_POI_COMM,
          tm_POI_RESIDUAL,
          tm_POI_SRC_COMM,
        tm_UPDATE_VEC,
    tm_MONITOR,
    tm_OUTPUT,
};

#ifdef PROF_MAPROF

#include "maprof.h"
#include "mpi.h"
#include <cstdio>

#define TIME_START(id)       maprof_time_start(id)
#define TIME_STOP(id)        maprof_time_stop(id)
//#define TIME_PRINT(id, str) mpprof_print_time_mpi(id, str)
#define TIME_PRINT(id, str)  print_time(id, str)
#define ADD_FLOPS(id, ops)   maprof_add_fp_ops(id, ops)
#define FLOPS_PRINT(id, str)  print_flops(id, str)

inline void print_time(int id, const char *str) {
  double ave  = maprof_get_time(id, MAPROF_AVE);
  double min  = maprof_get_time(id, MAPROF_MIN);
  double max  = maprof_get_time(id, MAPROF_MAX);
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0) {
    printf("%s (ave) %12.6e / (min) %12.6e / (max) %12.6e [sec]\n", str, ave, min, max);
  }
}

inline void print_flops(int id, const char *str) {
  double ave  = maprof_get_flops(id, MAPROF_AVE);
  double min  = maprof_get_flops(id, MAPROF_MIN);
  double max  = maprof_get_flops(id, MAPROF_MAX);
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0) {
    printf("%s (ave) %12.6e / (min) %12.6e / (max) %12.6e [Gflops]\n", str, ave, min, max);
  }
}

#elif defined MAPROF_FAPP

#include "fj_tool/fapp.h"
#define TIME_START(id)       fapp_start(#id, 0, 0)
#define TIME_STOP(id)        fapp_stop(#id, 0, 0)
#define TIME_PRINT(id, str)  ((void)0)
#define ADD_FLOPS(id, ops)   ((void)0)
#define FLOPS_PRINT(id, str) ((void)0)

#else

#define TIME_START(id)       ((void)0)
#define TIME_STOP(id)        ((void)0)
#define TIME_PRINT(id, str)  ((void)0)
#define ADD_FLOPS(id, ops)   ((void)0)
#define FLOPS_PRINT(id, str) ((void)0)

#endif

#endif // _TIMING_H_
