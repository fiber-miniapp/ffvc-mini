/*
 * Copyright (C) 2014 RIKEN AICS
 * This library is released under the terms of the MIT license.
 * http://fiber-miniapp.mit-license.org/
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "maprof_proc.h"

int main(int argc, char *argv[])
{
  int ch;

  while ((ch = getopt(argc, argv, "c:m:")) != -1) {
    switch (ch) {
    case 'c':
      maprof_read_cpuinfo(optarg);
      break;
    case 'm':
      maprof_read_meminfo(optarg);
      break;
    default:
      fprintf(stderr, "usage: %s [-c cpuinfo_file] [-m meminfo_file]\n", argv[0]);
      exit(1);
    }
  }
  
  printf("proc_name = %s\n", maprof_get_proc_name());
  printf("proc_clock_freq = %f GHz\n", maprof_get_proc_clock_freq());
  printf("num_core_node = %d\n", maprof_get_num_core_node());
  printf("num_proc_node = %d\n", maprof_get_num_proc_node());
  printf("mem_node = %f GB\n", maprof_get_mem_node());

  return 0;
}
