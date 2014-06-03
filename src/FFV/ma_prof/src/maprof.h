/*
 * Copyright (C) 2014 RIKEN AICS
 * This library is released under the terms of the MIT license.
 * http://fiber-miniapp.mit-license.org/
 */

/** @file
 *  @brief C/C++ interface
 *  @defgroup c_interface C/C++ interface
 *  @{
 */
#ifndef MAPROF_H
#define MAPROF_H

#include <stdio.h>

/* ---------- measuring ---------- */

/** @defgroup measure mesuring
 *  @{
 */

/** type of statistical values */
typedef enum {
  MAPROF_ROOT,   /**< rank0 */
  MAPROF_AVE,    /**< average */
  MAPROF_MIN,    /**< minimum */
  MAPROF_MAX,    /**< maximum */
  MAPROF_SD,     /**< standard deviation */
} maprof_stat_type;

#ifdef __cplusplus
extern "C" {
#endif

/** 
 * start timer.
 * @param[in] id  section ID
 */
void maprof_time_start(int id);

/**
 * stop timer.
 * @param[in] id  section ID
 */
void maprof_time_stop(int id);

/** 
 * add flop count.
 * @param[in] id  section ID
 * @param[in] ops  flop count
 */
void maprof_add_fp_ops(int id, double ops);

/**
 * add load count.
 * @param[in] id  section ID
 * @param[in] ops  load count
 */
void maprof_add_ld_ops(int id, double ops);

/**
 * add store count.
 * @param[in] id  section ID
 * @param[in] ops  store count
 */
void maprof_add_st_ops(int id, double ops);

/**
 * add effective load count.
 * @param[in] id  section ID
 * @param[in] ops  load count
 */
void maprof_add_ld_min_ops(int id, double ops);

/**
 * add effective store count.
 * @param[in] id  section ID
 * @param[in] ops  store count
 */
void maprof_add_st_min_ops(int id, double ops);

/**
 * get time.
 * @param[in] id  section ID
 * @param[in] type  statistical type
 * @return time
 */
double maprof_get_time(int id, maprof_stat_type type);

/**
 * get flop count.
 * @param[in] id  section ID
 * @param[in] type  statistical type
 * @return flop count
 */
double maprof_get_flops(int id, maprof_stat_type type);

/**
 * get memory throughput.
 * @param[in] id  section ID
 * @param[in] type  statistical type
 * @return memory throughput in unit of GB/s
 */
double maprof_get_throughput(int id, maprof_stat_type type);

/**
 * get effective memory throughput.
 * @param[in] id  section ID
 * @param[in] type  statistical type
 * @return memory throughput in unit of GB/s
 */
double maprof_get_effective_throughput(int id, maprof_stat_type type);

/**
 * print all information.
 * @param[in] id  section ID
 * @param[in] name  description
 */
void maprof_print(int id, const char *name);

/**
 * print time.
 * @param[in] id  section ID
 * @param[in] name  description
 */
void maprof_print_time(int id, const char *name);

/**
 * print time for MPI.
 * @param[in] id  section ID
 * @param[in] name  description
 */
void maprof_print_time_mpi(int id, const char *name);

/**
 * print time for MPI (in detail).
 * @param[in] id  section ID
 * @param[in] name  description
 */
void maprof_print_time_mpi_full(int id, const char *name);

/**
 * print to file stream other than stdout.
 * @param[in] file output stream
 */
void maprof_print_to_file(FILE *file);

/** @} */

/* ---------- reporting ---------- */

/** @defgroup yaml_report YAML reporting
 *  @{
 */

/**
 * setup.
 * @param[in] app_name  applicaton name
 * @param[in] app_version  application version
 */
void maprof_setup(const char *app_name, const char *app_version);

/**
 * output YAML file.
 */
void maprof_output();

/**
 * add key-value(string) pair to 'app'.
 * @param[in] key  key
 * @param[in] str  value(string)
 */
void maprof_app_add_str(const char *key, const char *str);

/**
 * add key-value(int) pair to 'app'.
 * @param[in] key  key
 * @param[in] n  value(int)
 */
void maprof_app_add_int(const char *key, int n);

/**
 * add key-value(float) pair to 'app'.
 * @param[in] key  key
 * @param[in] r  value(float)
 */
void maprof_app_add_float(const char *key, double r);

/**
 * add measuring section.
 * @param[in] name  section name
 * @param[in] id  section id
 */
void maprof_add_section(const char *name, int id);

/**
 * add problem size to 'profiles'
 * @param[in] key  name
 * @param[in] n  size
 */
void maprof_profile_add_problem_size(const char *key, int n);

/**
 * add key-value pair(string) to 'profiles'
 * @param[in] key  key
 * @param[in] str  value(string)
 */
void maprof_profile_add_str(const char *key, const char *str);

/**
 * add key-value pair(int) to 'profiles'
 * @param[in] key  key
 * @param[in] n  value(int)
 */
void maprof_profile_add_int(const char *key, int n);

/**
 * add key-value pair(float) to 'profiles'
 * @param[in] key  key
 * @param[in] r  value(float)
 */
void maprof_profile_add_float(const char *key, double r);

/** @} */

/* ---------- for fortran interface ---------- */

/** @defgroup fortran_helper fortan interface helper
 *  @{
 */

/**
 * set number of threads.
 * @param[in] n  number of threads
 */
void maprof_set_num_threads(int n);

/**
 * flush stdout stream.
 */
void maprof_flush_stdout();

/** @} */

#ifdef __cplusplus
}
#endif

/** @} */

#endif /* MAPROF_H */
