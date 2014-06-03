/*
 * Copyright (C) 2014 RIKEN AICS
 * This library is released under the terms of the MIT license.
 * http://fiber-miniapp.mit-license.org/
 */

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "maprof_yaml.h"

typedef enum {
  SCALAR,
  SEQUENCE,
  MAPPING,
} node_type_t;

typedef enum {
  STRING,
  INTEGER,
  FLOAT,
} scalar_type_t;

typedef maprof_yaml_style_type_t style_type_t;
typedef struct maprof_yaml_node_s node_t;

typedef struct seq_item_s seq_item_t;
struct seq_item_s {
  node_t *node;
  seq_item_t *next;
};

typedef struct map_item_s map_item_t;
struct map_item_s {
  const char *name;
  node_t *node;
  map_item_t *next;
};

struct maprof_yaml_node_s {
  node_type_t type;

  union {
    struct {
      scalar_type_t type;
      union {
        const char *str;
        int n;
        double r;
      } value;
    } scalar;

    struct {
      style_type_t style;
      seq_item_t *start;
      seq_item_t *end;
    } sequence;

    struct {
      style_type_t style;
      map_item_t *start;
      map_item_t *end;
    } mapping;
  } data;

};


maprof_yaml_node maprof_yaml_str_node(const char *s)
{
  node_t *n = (node_t *)malloc(sizeof(node_t));
  n->type = SCALAR;
  n->data.scalar.type = STRING;
  n->data.scalar.value.str = strdup(s);
  return n;
}


maprof_yaml_node maprof_yaml_quoted_str_node(const char *s)
{
  char *p = (char *)malloc(strlen(s) + 2 + 1);
  node_t *n = (node_t *)malloc(sizeof(node_t));
  n->type = SCALAR;
  n->data.scalar.type = STRING;
  sprintf(p, "\"%s\"", s);
  n->data.scalar.value.str = p;
  return n;
}


maprof_yaml_node maprof_yaml_int_node(int i)
{
  node_t *n = (node_t *)malloc(sizeof(node_t));
  n->type = SCALAR;
  n->data.scalar.type = INTEGER;
  n->data.scalar.value.n = i;
  return n;
}


maprof_yaml_node maprof_yaml_float_node(double f)
{
  node_t *n = (node_t *)malloc(sizeof(node_t));
  n->type = SCALAR;
  n->data.scalar.type = FLOAT;
  n->data.scalar.value.r = f;
  return n;
}


maprof_yaml_node maprof_yaml_seq_node(maprof_yaml_style_type_t style)
{
  node_t *n = (node_t *)malloc(sizeof(node_t));
  n->type = SEQUENCE;
  n->data.sequence.style = style;
  n->data.sequence.start = NULL;
  n->data.sequence.end = NULL;
  return n;
}


maprof_yaml_node maprof_yaml_map_node(maprof_yaml_style_type_t style)
{
  node_t *n = (node_t *)malloc(sizeof(node_t));
  n->type = MAPPING;
  n->data.mapping.style = style;
  n->data.mapping.start = NULL;
  n->data.mapping.end = NULL;
  return n;
}


void maprof_yaml_add_seq_item(maprof_yaml_node seq, maprof_yaml_node node) 
{
  seq_item_t *seq_end, *seq_new;
  assert(seq->type == SEQUENCE);

  seq_end = seq->data.sequence.end;
  seq_new = (seq_item_t *)malloc(sizeof(seq_item_t));
  seq_new->node = node;
  seq_new->next = NULL;
  if (seq_end == NULL) {
    assert(seq->data.sequence.start == NULL);
    seq->data.sequence.start = seq_new;
  } else {
    seq_end->next = seq_new;
  }
  seq->data.sequence.end = seq_new;
}


void maprof_yaml_add_map_item(maprof_yaml_node map, const char *name, maprof_yaml_node node)
{
  map_item_t *map_end, *map_new;
  assert(map->type == MAPPING);

  map_end = map->data.mapping.end;
  map_new = (map_item_t *)malloc(sizeof(map_item_t));
  map_new->name = strdup(name);
  map_new->node = node;
  map_new->next = NULL;
  if (map_end == NULL) {
    assert(map->data.mapping.start == NULL);
    map->data.mapping.start = map_new;
  } else {
    map_end->next = map_new;
  }
  map->data.mapping.end = map_new;
}


static void indent(int level, FILE *f)
{
  int i;
  for (i = 0; i < 2*level; i++) fputc(' ', f);
}


static void print_scalar_node(node_t *node, int level, style_type_t style, FILE *f)
{
  switch (node->data.scalar.type) {
  case STRING:
    fprintf(f, "%s", node->data.scalar.value.str);
    break;
  case INTEGER:
    fprintf(f, "%d", node->data.scalar.value.n);
    break;
  case FLOAT:
    fprintf(f, "%g", node->data.scalar.value.r);
    break;
  }
  if (style == MAPROF_YAML_BLOCK_STYLE) fputs("\n", f);
}


static void print_node(node_t *node, int level, style_type_t style, FILE *f);


static void print_sequence_node(node_t *node, int level, style_type_t style, FILE *f)
{
  seq_item_t *p = node->data.sequence.start;
  int flow_style_new_line = 0;

  if (style == MAPROF_YAML_BLOCK_STYLE) {
    if (node->data.sequence.style == MAPROF_YAML_FLOW_STYLE) {
      style = MAPROF_YAML_FLOW_STYLE;
      flow_style_new_line = 1;
    }
  }

  if (style == MAPROF_YAML_BLOCK_STYLE) {
    if (level > 0) fputs("\n", f);
    while (p != NULL) {
      indent(level, f);
      fputs("- ", f);
      print_node(p->node, level+1, style, f);
      p = p->next;
    }
  } else {
    /* FLOW STYLE */
    fputs("[", f);
    while (p != NULL) {
      print_node(p->node, level+1, style, f);
      p = p->next;
      if (p != NULL) fputs(", ", f);
    }
    fputs("]", f);
    if (flow_style_new_line) fputs("\n", f);
  }
}


static void print_mapping_node(node_t *node, int level, style_type_t style, FILE *f)
{
  map_item_t *p = node->data.mapping.start;
  int flow_style_new_line = 0;

  if (style == MAPROF_YAML_BLOCK_STYLE) {
    if (node->data.mapping.style == MAPROF_YAML_FLOW_STYLE) {
      style = MAPROF_YAML_FLOW_STYLE;
      flow_style_new_line = 1;
    }
  }

  if (style == MAPROF_YAML_BLOCK_STYLE) {
    if (level > 0) fputs("\n", f);
    while (p != NULL) {
      indent(level, f);
      fputs(p->name, f);
      fputs(": ", f);
      print_node(p->node, level+1, style, f);
      p = p->next;
    }
  } else {
    /* FLOW STYLE */
    fputs("{", f);
    while (p != NULL) {
      fputs(p->name, f);
      fputs(": ", f);
      print_node(p->node, level+1, style, f);
      p = p->next;
      if (p != NULL) fputs(", ", f);
    }
    fputs("}", f);
    if (flow_style_new_line) fputs("\n", f);
  }
}


static void print_node(node_t *node, int level, style_type_t style, FILE *f)
{
  switch (node->type) {
  case SCALAR:
    print_scalar_node(node, level, style, f);
    break;
  case SEQUENCE:
    print_sequence_node(node, level, style, f);
    break;
  case MAPPING:
    print_mapping_node(node, level, style, f);
    break;
  }
}


void maprof_yaml_print(maprof_yaml_node root, FILE *f)
{
  assert(root->type == MAPPING);
  print_node(root, 0, MAPROF_YAML_BLOCK_STYLE, f);
}
