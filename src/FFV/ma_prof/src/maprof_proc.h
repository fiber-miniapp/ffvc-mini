/*
 * Copyright (C) 2014 RIKEN AICS
 * This library is released under the terms of the MIT license.
 * http://fiber-miniapp.mit-license.org/
 */

#ifndef MAPROF_PROC_H
#define MAPROF_PROC_H

const char *maprof_get_proc_name();

double maprof_get_proc_clock_freq();

int maprof_get_num_core_node();

int maprof_get_num_proc_node();

double maprof_get_mem_node();

void maprof_read_cpuinfo(const char *cpuinfo);  /* for debug */

void maprof_read_meminfo(const char *meminfo);  /* for debug */


#endif /* MAPROF_PROC_H */
