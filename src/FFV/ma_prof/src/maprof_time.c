/*
 * Copyright (C) 2014 RIKEN AICS
 * This library is released under the terms of the MIT license.
 * http://fiber-miniapp.mit-license.org/
 */

#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#ifdef USE_MPI
#include "mpi.h"
#endif
#include "maprof.h"


typedef struct {
  int id;

  /* timing */
  double start;
  double time;
  long count;
  int started;

  /* operation counters */
  double fp_ops;      /* floating point */
  double ld_ops;      /* load */
  double st_ops;      /* store */
  double ld_min_ops;  /* load (effective) */
  double st_min_ops;  /* store (effective) */
} Section;


static Section sections[MAPROF_MAX_SECTIONS];


static FILE *output_file = NULL;
#define OUTPUT_STREAM  (output_file ? output_file : stdout)


static double get_current_time()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return (double)tv.tv_sec + (double)tv.tv_usec * 1.0e-6;
}


static void error_exit(const char *fmt, ...)
{
  const int error_code = 1;

  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);

#ifdef USE_MPI
  MPI_Abort(MPI_COMM_WORLD, error_code);
#else
  exit(error_code);
#endif 
}


static double get_value(double val, maprof_stat_type type)
{
  switch (type) {
    case MAPROF_ROOT:
      return val;
#ifdef USE_MPI
    case MAPROF_AVE:
    {
      int nprocs;
      double sum = 0.0;
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
      MPI_Reduce(&val, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      return sum / nprocs;
    }
    case MAPROF_MIN:
    {
      double min = 0.0;
      MPI_Reduce(&val, &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      return min;
    }
    case MAPROF_MAX:
    {
      double max = 0.0;
      MPI_Reduce(&val, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      return max;
    }
    case MAPROF_SD:
    {
      int nprocs;
      double sum2 = 0;
      double val2= val*val;
      double ave = get_value(val, MAPROF_AVE);
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
      MPI_Reduce(&val2, &sum2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      return sqrt(sum2/nprocs - ave*ave);
    }
#endif
    default:
      error_exit("***error: get_value: unknown Value Type<%d>\n", type);
  }
  /* NOTREACHED */
}


void maprof_time_start(int id)
{
  sections[id].started = 1;
  sections[id].start = get_current_time();
}


void maprof_time_stop(int id)
{
  if (!sections[id].started) error_exit("***error: maprof_stop: section<%d> not started\n", id);
  sections[id].time += get_current_time() - sections[id].start;
  sections[id].count++;
  sections[id].started = 0;
}


void maprof_add_fp_ops(int id, double ops)
{
  sections[id].fp_ops += ops;
}


void maprof_add_ld_ops(int id, double ops)
{
  sections[id].ld_ops += ops;
}


void maprof_add_st_ops(int id, double ops)
{
  sections[id].st_ops += ops;
}


void maprof_add_ld_min_ops(int id, double ops)
{
  sections[id].ld_min_ops += ops;
}


void maprof_add_st_min_ops(int id, double ops)
{
  sections[id].st_min_ops += ops;
}


double maprof_get_time(int id, maprof_stat_type type)
{
  return get_value(sections[id].time, type);
}


double maprof_get_flops(int id, maprof_stat_type type)
{
  double val_local;
  if (sections[id].time > 0.0) {
    val_local = sections[id].fp_ops / sections[id].time 
                * 1.0e-9;  /* in unit of GFLOPS */
  } else {
    val_local = 0.0;
  }
  return get_value(val_local, type);
}


double maprof_get_throughput(int id, maprof_stat_type type)
{
  double val_local;
  if (sections[id].time > 0.0) {
    val_local = (sections[id].ld_ops + sections[id].st_ops) / sections[id].time
                * 8 * 1.0e-9;  /* in unit of GB/s */
  } else {
    val_local = 0.0;
  }
  return get_value(val_local, type);
}


double maprof_get_effective_throughput(int id, maprof_stat_type type)
{
  double val_local;
  if (sections[id].time > 0.0) {
    val_local = (sections[id].ld_min_ops + sections[id].st_min_ops) / sections[id].time 
                 * 8 * 1.0e-9;  /* in unit of GB/s */
  } else {
    val_local = 0.0;
  }
  return get_value(val_local, type);
}


void maprof_print(int id, const char *name)
{
  FILE *out = OUTPUT_STREAM;
  double time, flops, tput, etput;
  int myrank = 0;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif
  if (myrank != 0) return;

  time  = maprof_get_time(id, MAPROF_ROOT);
  flops = maprof_get_flops(id, MAPROF_ROOT);
  tput  = maprof_get_throughput(id, MAPROF_ROOT);
  etput = maprof_get_effective_throughput(id, MAPROF_ROOT);

  if (time > 0.0) {
    fprintf(out, "%s Time:                 %g (S)\n", name, time);
  } else {
    fprintf(out, "%s Time:                 unknown\n", name);
  }
  if (time > 0.0 && flops > 0.0) {
    fprintf(out, "%s FLOPS:                %g (GFLOPS)\n", name, flops);
  } else {
    fprintf(out, "%s FLOPS:                unknown\n", name);
  }
  if (time > 0.0 && tput > 0.0) {
    fprintf(out, "%s Throughput:           %g (GB/S)\n", name, tput);
  } else {
    fprintf(out, "%s Throughput:           unknown\n", name);
  }
  if (time > 0.0 && etput > 0.0) {
    fprintf(out, "%s Effective Throughput: %g (GB/S)\n", name, etput);
  } else {
    fprintf(out, "%s Effective Throughput: unknown\n", name);
  }
}


void maprof_print_time(int id, const char *name)
{
  FILE *out = OUTPUT_STREAM;
  int myrank = 0;
  double time = maprof_get_time(id, MAPROF_ROOT);
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif
  if (myrank == 0) fprintf(out, "%s %f [sec]\n", name, time);
}


#ifdef USE_MPI
void maprof_print_time_mpi(int id, const char *name)
{
  FILE *out = OUTPUT_STREAM;
  int myrank;
  double time_root = maprof_get_time(id, MAPROF_ROOT);
  double time_min  = maprof_get_time(id, MAPROF_MIN);
  double time_max  = maprof_get_time(id, MAPROF_MAX);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0) {
    fprintf(out, "%s (rank0)%f / (min)%f / (max)%f [sec]\n",
           name, time_root, time_min, time_max);
  }
}


void maprof_print_time_mpi_full(int id, const char *name)
{
  FILE *out = OUTPUT_STREAM;
  int myrank;
  double time_root = maprof_get_time(id, MAPROF_ROOT);
  double time_min  = maprof_get_time(id, MAPROF_MIN);
  double time_max  = maprof_get_time(id, MAPROF_MAX);
  double time_ave  = maprof_get_time(id, MAPROF_AVE);
  double time_sd   = maprof_get_time(id, MAPROF_SD);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0) {
    fprintf(out, "%s\n", name);
    fprintf(out, "    rank0: %f [sec]\n", time_root);
    fprintf(out, "    min:   %f [sec]\n", time_min);
    fprintf(out, "    max:   %f [sec]\n", time_max);
    fprintf(out, "    ave:   %f [sec]\n", time_ave);
    fprintf(out, "    sd:    %f [sec]\n", time_sd);
  }
}
#endif


void maprof_print_to_file(FILE *file)
{
  output_file = file;
}
