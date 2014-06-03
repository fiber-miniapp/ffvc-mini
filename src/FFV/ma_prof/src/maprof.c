/*
 * Copyright (C) 2014 RIKEN AICS
 * This library is released under the terms of the MIT license.
 * http://fiber-miniapp.mit-license.org/
 */

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "maprof.h"
#include "maprof_yaml.h"
#include "maprof_proc.h"

#ifdef USE_MPI
#include "mpi.h"
#endif
#ifdef _OPENMP
#include "omp.h"
#endif

static maprof_yaml_node Node_root = 0;
static maprof_yaml_node Node_app = 0;
static maprof_yaml_node Node_system = 0;
static maprof_yaml_node Node_compilers = 0;
static maprof_yaml_node Node_sections = 0;
static maprof_yaml_node Node_profile = 0;
static maprof_yaml_node Node_profile_problem_size = 0;

typedef struct {
  const char *name;
  int id;
} Section;

static Section Sections[MAPROF_MAX_SECTIONS];

static int N_sections = 0;

static int N_threads = 0;


static const char *get_host_name()
{
  static char name[256];
  gethostname(name, sizeof(name));
  return name;
}


static const char *get_timestamp()
{
//static char buf[20];
  static char buf[26];
  time_t t = time(NULL);
  struct tm *l = localtime(&t);
#if 0
  sprintf(buf, "%4d-%02d-%02d %02d:%02d:%02d",
          l->tm_year+1900, l->tm_mon+1, l->tm_mday,
          l->tm_hour, l->tm_min, l->tm_sec);
#endif
  strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S %z", l);
  return buf;
}


static int get_num_processes()
{
#ifdef USE_MPI
  int n;
  MPI_Comm_size(MPI_COMM_WORLD, &n);
  return n;
#else
  return 1;
#endif
}


static int get_num_threads()
{
  if (N_threads > 0) return N_threads;
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}


static void app_setup(const char *name, const char *version)
{
  Node_app = maprof_yaml_map_node(MAPROF_YAML_BLOCK_STYLE);
  maprof_yaml_add_map_item(Node_root, "app", Node_app);
  maprof_yaml_add_map_item(Node_app, "name", maprof_yaml_str_node(name));
  maprof_yaml_add_map_item(Node_app, "version", maprof_yaml_str_node(version));
}


static void system_setup()
{
  Node_system = maprof_yaml_map_node(MAPROF_YAML_BLOCK_STYLE);
  maprof_yaml_add_map_item(Node_root, "system", Node_system);
  maprof_yaml_add_map_item(Node_system, "name", maprof_yaml_str_node(get_host_name()));
  maprof_yaml_add_map_item(Node_system, "proc_name", maprof_yaml_str_node(maprof_get_proc_name()));
  maprof_yaml_add_map_item(Node_system, "proc_clock_freq", maprof_yaml_float_node(maprof_get_proc_clock_freq()));
  maprof_yaml_add_map_item(Node_system, "num_core_node", maprof_yaml_int_node(maprof_get_num_core_node()));
  maprof_yaml_add_map_item(Node_system, "num_proc_node", maprof_yaml_int_node(maprof_get_num_proc_node()));
  maprof_yaml_add_map_item(Node_system, "mem_node", maprof_yaml_float_node(maprof_get_mem_node()));
}


static void compilers_setup()
{
  maprof_yaml_node node_lang;

  Node_compilers = maprof_yaml_seq_node(MAPROF_YAML_BLOCK_STYLE);
  maprof_yaml_add_map_item(Node_root, "compilers", Node_compilers);
#ifdef MAPROF_FC
  node_lang = maprof_yaml_map_node(MAPROF_YAML_BLOCK_STYLE);
  maprof_yaml_add_seq_item(Node_compilers, node_lang);
  maprof_yaml_add_map_item(node_lang, "lang", maprof_yaml_str_node("Fortran"));
  maprof_yaml_add_map_item(node_lang, "name", maprof_yaml_str_node(MAPROF_FC));
  /*
  maprof_yaml_add_map_item(node_lang, "version", maprof_yaml_str_node(MAPROF_FC_VERSION));
  */
  maprof_yaml_add_map_item(node_lang, "version", maprof_yaml_quoted_str_node(MAPROF_FC_VERSION));
  maprof_yaml_add_map_item(node_lang, "options", maprof_yaml_str_node(MAPROF_FFLAGS));
#endif
#ifdef MAPROF_CC
  node_lang = maprof_yaml_map_node(MAPROF_YAML_BLOCK_STYLE);
  maprof_yaml_add_seq_item(Node_compilers, node_lang);
  maprof_yaml_add_map_item(node_lang, "lang", maprof_yaml_str_node("C"));
  maprof_yaml_add_map_item(node_lang, "name", maprof_yaml_str_node(MAPROF_CC));
  /*
  maprof_yaml_add_map_item(node_lang, "version", maprof_yaml_str_node(MAPROF_CC_VERSION));
  */
  maprof_yaml_add_map_item(node_lang, "version", maprof_yaml_quoted_str_node(MAPROF_CC_VERSION));
  maprof_yaml_add_map_item(node_lang, "options", maprof_yaml_str_node(MAPROF_CFLAGS));
#endif
#ifdef MAPROF_CXX
  node_lang = maprof_yaml_map_node(MAPROF_YAML_BLOCK_STYLE);
  maprof_yaml_add_seq_item(Node_compilers, node_lang);
  maprof_yaml_add_map_item(node_lang, "lang", maprof_yaml_str_node("C++"));
  maprof_yaml_add_map_item(node_lang, "name", maprof_yaml_str_node(MAPROF_CXX));
  /*
  maprof_yaml_add_map_item(node_lang, "version", maprof_yaml_str_node(MAPROF_CXX_VERSION));
  */
  maprof_yaml_add_map_item(node_lang, "version", maprof_yaml_quoted_str_node(MAPROF_CXX_VERSION));
  maprof_yaml_add_map_item(node_lang, "options", maprof_yaml_str_node(MAPROF_CXXFLAGS));
#endif
}


static void sections_setup()
{
  int i;

  Node_sections = maprof_yaml_seq_node(MAPROF_YAML_FLOW_STYLE);
  maprof_yaml_add_map_item(Node_root, "sections", Node_sections);
  for (i = 0; i < N_sections; i++) {
    maprof_yaml_add_seq_item(Node_sections, maprof_yaml_str_node(Sections[i].name));
  }
}


static void profile_setup()
{
  maprof_yaml_node node_profiles = maprof_yaml_seq_node(MAPROF_YAML_BLOCK_STYLE);
  maprof_yaml_add_map_item(Node_root, "profiles", node_profiles);

  Node_profile = maprof_yaml_map_node(MAPROF_YAML_BLOCK_STYLE);
  maprof_yaml_add_seq_item(node_profiles, Node_profile);

  maprof_yaml_add_map_item(Node_profile, "date_time", maprof_yaml_str_node(get_timestamp()));
  maprof_yaml_add_map_item(Node_profile, "np", maprof_yaml_int_node(get_num_processes()));
  maprof_yaml_add_map_item(Node_profile, "nt", maprof_yaml_int_node(get_num_threads()));

  Node_profile_problem_size = maprof_yaml_map_node(MAPROF_YAML_BLOCK_STYLE);
  maprof_yaml_add_map_item(Node_profile, "problem_size", Node_profile_problem_size);
}


static void profile_sections_setup()
{
  int i;
  maprof_yaml_node node_sections, node_sec;

  node_sections = maprof_yaml_seq_node(MAPROF_YAML_BLOCK_STYLE);
  maprof_yaml_add_map_item(Node_profile, "sections", node_sections);

  for (i = 0; i < N_sections; i++) {
    int id = Sections[i].id;
    const char *name = Sections[i].name;
    node_sec= maprof_yaml_map_node(MAPROF_YAML_BLOCK_STYLE);
    maprof_yaml_add_seq_item(node_sections, node_sec);
    maprof_yaml_add_map_item(node_sec, "name", maprof_yaml_str_node(name));
    maprof_yaml_add_map_item(node_sec, "time",
                             maprof_yaml_float_node(maprof_get_time(id, MAPROF_ROOT)));
    maprof_yaml_add_map_item(node_sec, "flop",
                             maprof_yaml_float_node(maprof_get_flops(id, MAPROF_ROOT)));
    maprof_yaml_add_map_item(node_sec, "tput",
                             maprof_yaml_float_node(maprof_get_throughput(id, MAPROF_ROOT)));
  }
}


void maprof_setup(const char *app_name, const char *app_version)
{
  Node_root = maprof_yaml_map_node(MAPROF_YAML_BLOCK_STYLE);

  app_setup(app_name, app_version);
  system_setup();
  compilers_setup();
  sections_setup();
  profile_setup();
}


void maprof_output()
{
  const char *file;
  FILE *f;
#if USE_MPI
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank != 0) return;
#endif

  profile_sections_setup();

  if ((file = getenv("MAPROF_OUTPUT")) ==NULL) {
#ifdef MAPROF_OUTPUT
    file = MAPROF_OUTPUT;
#else
    file = "stdout";
#endif
  }

  if (strcmp(file, "stdout") == 0) {
    maprof_yaml_print(Node_root, stdout);
  } else {
    if ((f = fopen(file, "w")) == NULL) {
      perror(file);
      exit(EXIT_FAILURE);
    }
    maprof_yaml_print(Node_root, f);
    fclose(f);
  }
}


void maprof_add_section(const char *name, int id)
{
  Sections[N_sections].name = strdup(name);
  Sections[N_sections].id = id;
  N_sections++;

  if (Node_sections > 0) {
    maprof_yaml_add_seq_item(Node_sections, maprof_yaml_str_node(name));
  }
}


void maprof_app_add_str(const char *key, const char *str)
{
  maprof_yaml_add_map_item(Node_app, key, maprof_yaml_str_node(str));
}


void maprof_app_add_int(const char *key, int n)
{
  maprof_yaml_add_map_item(Node_app, key, maprof_yaml_int_node(n));
}


void maprof_app_add_float(const char *key, double r)
{
  maprof_yaml_add_map_item(Node_app, key, maprof_yaml_float_node(r));
}


void maprof_profile_add_problem_size(const char *key, int n)
{
  maprof_yaml_add_map_item(Node_profile_problem_size, key, maprof_yaml_int_node(n));
}


void maprof_profile_add_str(const char *key, const char *str)
{
  maprof_yaml_add_map_item(Node_profile, key, maprof_yaml_str_node(str));
}


void maprof_profile_add_int(const char *key, int n)
{
  maprof_yaml_add_map_item(Node_profile, key, maprof_yaml_int_node(n));
}


void maprof_profile_add_float(const char *key, double r)
{
  maprof_yaml_add_map_item(Node_profile, key, maprof_yaml_float_node(r));
}

void maprof_set_num_threads(int n)
{
  N_threads = n;
}

void maprof_flush_stdout()
{
  fflush(stdout);
}
