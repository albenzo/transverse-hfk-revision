#ifndef TRANSVERSE_HFK_H
#define TRANSVERSE_HFK_H
#include <limits.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "states.h"

typedef int (*printf_t)(const char *format, ...);
printf_t print_ptr;

#define SILENT 0
#define QUIET 1
#define VERBOSE 2

int get_verbosity(void);
void set_verbosity(const int);

int mod(const int, const int);
int mod_up(const int, const int);
int min(const int, const int);

StateList swap_cols(const int, const int, const State, const Grid_t *const);
int null_homologous_D0Q(const State, const Grid_t *const);
int null_homologous_D1Q(const State, const Grid_t *const);
int null_homologous_lift(const LiftState, const LiftGrid_t * const);

VertexList prepend_vertex(const int, const VertexList);
EdgeList prepend_edge(const int, const int, const EdgeList);
EdgeList create_edge(const int, const int);
void free_edge_list(const EdgeList);
EdgeList add_mod_two_lists(const VertexList, const VertexList, const EdgeList *const);
EdgeList append_ordered(const int, const int, const EdgeList);
void special_homology(const int, const int, EdgeList *);
void contract(const int, const int, EdgeList *);

StateList new_rectangles_out_of(const StateList, const State, const Grid_t *const);
StateList new_rectangles_into(const StateList, const State, const Grid_t *const);
StateList fixed_wt_rectangles_out_of(const int, const State, const Grid_t *const);
LiftStateList new_lift_rectangles_out_of(const LiftStateList, const LiftState, const LiftGrid_t * const);
LiftStateList new_lift_rectangles_into(const LiftStateList, const LiftState, const LiftGrid_t * const);

void print_state(const State, const Grid_t *const);
void print_lift_state(const LiftState, const LiftGrid_t * const);
void print_state_short(const State, const Grid_t *const);
void print_lift_state_short(const LiftState, const LiftGrid_t * const);
void print_states(const StateList, const Grid_t *const);
void print_lift_states(const LiftStateList, const Grid_t * const);
void print_edges(const EdgeList);
void print_math_edges(const EdgeList);
void print_math_edges_a(const EdgeList);
void print_vertices(const VertexList);
void print_grid_perm(const Grid_t *const G);
void print_grid(const Grid_t *const G);
void print_lift_grid(const Grid_t *const G);
void print_tb_r(const Grid_t *const G);
void print_2AM(const Grid_t *const G, int plus);

#endif
