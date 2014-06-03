/*
 * Copyright (C) 2014 RIKEN AICS
 * This library is released under the terms of the MIT license.
 * http://fiber-miniapp.mit-license.org/
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "maprof_proc.h"

#define BUF_SIZE 1024

static char *proc_name;
static double proc_clock_freq;
static int num_core_node;
static int num_proc_node;
static double mem_node;

static int cpuinfo_first = 1;
static int meminfo_first = 1;


static char *trim(char *s)
{
  int i;

  for (i = strlen(s)-1; i >= 0; i--) {
    if (isspace(s[i])) {
      s[i] = '\0';
    } else {
      break;
    }
  }

  for (i = 0; i < strlen(s); i++) {
    if (isalnum(s[i])) break;
  }

  return &s[i];
}


static int scan_line(FILE *f, char **key, char **val)
{
  static char buf[BUF_SIZE];
  char *k, *v, *p;

  do {
    if (fgets(buf, sizeof(buf), f) == NULL) return 0;
    p = strchr(buf, ':');
    if (p != NULL) {
      *p = '\0';
      p++;
      k = trim(buf);
      v = trim(p);
    } else {
      k = trim(buf);
      v = NULL;
    }
  } while (strlen(k) == 0);

  *key = k;
  *val = v;

  return 1;
}


static void cpuinfo_unknown()
{
#ifdef PROC_NAME
  proc_name = strdup(PROC_NAME);
#else
  proc_name = strdup("unknown");
#endif
#ifdef NUM_CORE_NODE
  num_core_node = NUM_CORE_NODE;
#else
  num_core_node = 0;
#endif
#ifdef NUM_PROC_NODE
  num_proc_node = NUM_PROC_NODE;
#else
  num_proc_node = 0;
#endif
#ifdef PROC_CLOCK_FREQ
  proc_clock_freq = PROC_CLOCK_FREQ;
#else
  proc_clock_freq = 0.0;
#endif
}


static void meminfo_unknown()
{
#ifdef MEM_NODE
  mem_node = MEM_NODE;
#else
  mem_node = 0.0;
#endif
}


static void cpuinfo_intel(FILE *f)
{
  char *key, *val;
  int n_processor = 0;
  int n_sibling = 1;
  int n_core_proc = 1;

  rewind(f);

  while (scan_line(f, &key, &val)) {
    if (strcmp(key, "processor") == 0) {
      n_processor++;
    }
    else if (strcmp(key, "model name") == 0) {
      if (n_processor == 1) proc_name = strdup(val);
    }
    else if (strcmp(key, "cpu MHz") == 0) {
      if (n_processor == 1) {
        sscanf(val, "%lf", &proc_clock_freq);
        proc_clock_freq /= 1000.0;  /* MHz -> GHz */
      }
    }
    else if (strcmp(key, "siblings") == 0) {
      sscanf(val, "%d",  &n_sibling);
    }
    else if (strcmp(key, "cpu cores") == 0) {
      sscanf(val, "%d",  &n_core_proc);
    }
  }

  num_proc_node = n_processor / n_sibling;
  num_core_node = n_core_proc * num_proc_node;
}


static void cpuinfo_sparc(FILE *f)
{
  char *key, *val;
  long int freq;

  rewind(f);

  while (scan_line(f, &key, &val)) {
    if (strcmp(key, "cpu") == 0) {
      proc_name = strdup(val);
    }
    else if (strcmp(key, "ncpus probed") == 0) {
      sscanf(val, "%d", &num_core_node);
    }
    else if (strcmp(key, "Cpu0ClkTck") == 0) {
      sscanf(val, "%lx", &freq);
      proc_clock_freq = freq * 1.0e-9;   /* Hz -> GHz */
    }
  }

  num_proc_node = 1;  /* K and FX10 */
}


#ifdef __APPLE__
static void cpuinfo_mac(FILE *f)
{
  char *key, *val;
  int total_cores;

  while (scan_line(f, &key, &val)) {
    if (strcmp(key, "Processor Name") == 0) {
      proc_name = strdup(val);
    }
    else if (strcmp(key, "Processor Speed") == 0) {
      sscanf(val, "%lf GHz", &proc_clock_freq);
    }
    else if (strcmp(key, "Number of Processors") == 0) {
      sscanf(val, "%d", &num_proc_node);
    }
    else if (strcmp(key, "Total Number of Cores") == 0) {
      sscanf(val, "%d", &total_cores);
    }
  }
  num_core_node = total_cores / num_proc_node;
}
#endif


#ifdef __APPLE__
static void meminfo_mac(FILE *f)
{
  char *key, *val;

  while (scan_line(f, &key, &val)) {
    if (strcmp(key, "Memory") == 0) {
      sscanf(val, "%lf GB", &mem_node);
    }
  }
}
#endif


void maprof_read_cpuinfo(const char *cpuinfo)
{
  FILE *f;
  char *key, *val;

  cpuinfo_first = 0;

#ifdef __APPLE__
  if (cpuinfo == NULL) {
    f = popen("system_profiler SPHardwareDataType", "r");
    if (f == NULL) {
      perror(NULL);
      cpuinfo_unknown();
      return;
    }
    cpuinfo_mac(f);
    pclose(f);
    return;
  }
#else
  if (cpuinfo == NULL) cpuinfo = "/proc/cpuinfo";
#endif

  if ((f = fopen(cpuinfo, "r")) == NULL) {
    perror(cpuinfo);
    cpuinfo_unknown();
    return;
  }

  if  (scan_line(f, &key, &val)) {
    if (strcmp(key, "processor") == 0) {
      cpuinfo_intel(f);
    } else if (strcmp(key, "cpu") == 0) {
      cpuinfo_sparc(f);
    } else {
      cpuinfo_unknown();
    }
  } else {
    cpuinfo_unknown();
  }

  fclose(f);
}


void maprof_read_meminfo(const char *meminfo)
{
  FILE *f;
  char *key, *val;
  long int mem;

  meminfo_first = 0;

#ifdef __APPLE__
  if (meminfo == NULL) {
    f = popen("system_profiler SPHardwareDataType", "r");
    if (f == NULL) {
      perror(NULL);
      meminfo_unknown();
      return;
    }
    meminfo_mac(f);
    pclose(f);
    return;
  }
#else
  if (meminfo == NULL) meminfo = "/proc/meminfo";
#endif

  if ((f = fopen(meminfo, "r")) == NULL) {
    perror(meminfo);
    meminfo_unknown();
    return;
  }

  if (scan_line(f, &key, &val) &&
      strcmp(key, "MemTotal") == 0) {
    sscanf(val, "%ld kB", &mem);
    mem_node = mem * 1.0e-6;   /* kB -> GB */
  } else {
    meminfo_unknown();
  }

  fclose(f);
}


const char *maprof_get_proc_name()
{
  if (cpuinfo_first) maprof_read_cpuinfo(NULL);
  return proc_name;
}


double maprof_get_proc_clock_freq()
{
  if (cpuinfo_first) maprof_read_cpuinfo(NULL);
  return proc_clock_freq;
}


int maprof_get_num_core_node()
{
  if (cpuinfo_first) maprof_read_cpuinfo(NULL);
  return num_core_node;
}


int maprof_get_num_proc_node()
{
  if (cpuinfo_first) maprof_read_cpuinfo(NULL);
  return num_proc_node;
}


double maprof_get_mem_node()
{
  if (meminfo_first) maprof_read_meminfo(NULL);
  return mem_node;
}
