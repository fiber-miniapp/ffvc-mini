/*
 * Copyright (C) 2014 RIKEN AICS
 * This library is released under the terms of the MIT license.
 * http://fiber-miniapp.mit-license.org/
 */

#ifndef MAPROF_YAML_H
#define MAPROF_YAML_H

#include <stdio.h>

typedef enum {
  MAPROF_YAML_FLOW_STYLE,
  MAPROF_YAML_BLOCK_STYLE,
} maprof_yaml_style_type_t;

typedef struct maprof_yaml_node_s *maprof_yaml_node;

maprof_yaml_node maprof_yaml_str_node(const char *s);

maprof_yaml_node maprof_yaml_quoted_str_node(const char *s);

maprof_yaml_node maprof_yaml_int_node(int i);

maprof_yaml_node maprof_yaml_float_node(double f);

maprof_yaml_node maprof_yaml_seq_node(maprof_yaml_style_type_t style);

maprof_yaml_node maprof_yaml_map_node(maprof_yaml_style_type_t style);

void maprof_yaml_add_seq_item(maprof_yaml_node seq, maprof_yaml_node node);

void maprof_yaml_add_map_item(maprof_yaml_node map, const char *name, maprof_yaml_node node);

void maprof_yaml_print(maprof_yaml_node root, FILE *f);

#endif /* MAPROF_YAML_H */

