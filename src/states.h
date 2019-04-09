#ifndef STATES_H
#define STATES_H

#include <stdlib.h>
#include <string.h>

#define BLACK 0
#define RED 1

typedef char *State;

struct Grid {
  State Xs;
  State Os;
  int arc_index;
};

typedef struct Grid Grid_t;

typedef char **LiftState;
struct LiftGrid {
  State Xs;
  State Os;
  int arc_index;
  int sheets;
};

typedef struct LiftGrid LiftGrid_t;

typedef struct LiftStateRBTreeNode LiftStateRBTreeNode_t;
typedef LiftStateRBTreeNode_t * LiftStateRBTree;

struct LiftStateRBTreeNode {
  int color;
  int tag;
  LiftState data;
  LiftStateRBTree left;
  LiftStateRBTree right;
  LiftStateRBTree parent;
};

typedef struct StateRBTreeNode StateRBTreeNode_t;
typedef StateRBTreeNode_t * StateRBTree;

struct StateRBTreeNode {
  int color;
  int tag;
  State data;
  StateRBTree left;
  StateRBTree right;
  StateRBTree parent;
};

struct Vertex {
  int data;
  struct Vertex *nextVertex;
};

typedef struct Vertex Vertex_t;
typedef Vertex_t *VertexList;

struct EdgeNode {
  int start;
  int end;
  struct EdgeNode *nextEdge;
};

typedef struct EdgeNode EdgeNode_t;
typedef EdgeNode_t *EdgeList;

typedef struct StateNode StateNode_t;
typedef StateNode_t *StateList;

struct StateNode {
  State data;
  StateList nextState;
};

typedef struct LiftStateNode LiftStateNode_t;
typedef LiftStateNode_t *LiftStateList;

struct LiftStateNode {
  LiftState data;
  LiftStateList nextState;
};

void copy_state(State*, const State* const, const Grid_t * const);

void init_lift_state(LiftState*, const LiftGrid_t * const G);
void copy_lift_state(LiftState *, const LiftState * const, const LiftGrid_t * const);
void free_lift_state(LiftState*, const LiftGrid_t * const G);

int is_state(const State, const Grid_t *const);
int is_lift_state(const LiftState, const LiftGrid_t * const);
void mirror_lift_state(LiftState *, const LiftGrid_t * const);
int eq_state(const State a, const State b, const Grid_t * const G);
int eq_lift_state(const LiftState a, const LiftState b, const LiftGrid_t * const G);
int comp_state(const State a, const State b, const Grid_t * const G);
int comp_lift_state(const LiftState, const LiftState, const LiftGrid_t * const G);

int is_grid(const Grid_t *const);
int is_lift_grid(const LiftGrid_t * const);
LiftGrid_t * mirror_lift_grid(const LiftGrid_t * const);
int get_writhe(const Grid_t *const);
void cusps(int *, const Grid_t *const);

int NESW_pO(const State, const Grid_t *const);
int NESW_Op(const State, const Grid_t *const);
int NESW_pp(const State, const Grid_t *const);

int get_number(const State, const StateList, const Grid_t *const);
int get_lift_number(const LiftState, LiftStateList, const LiftGrid_t * const);
StateList remove_state(const State, const StateList, const Grid_t *const);
LiftStateList remove_lift_state(const LiftState, const LiftStateList, const LiftGrid_t * const);
void free_state_list(StateList);
void free_lift_state_list(LiftStateList, const LiftGrid_t * const);

void left_rotate(LiftStateRBTree*, LiftStateRBTree);
void right_rotate(LiftStateRBTree*, LiftStateRBTree);
void insert_data(LiftStateRBTree*, LiftState, const LiftGrid_t * const);
void insert_tagged_data(LiftStateRBTree*, LiftState, int, const LiftGrid_t * const);
void insert_node(LiftStateRBTree*, LiftStateRBTree, const LiftGrid_t * const);
void insert_fixup(LiftStateRBTree*, LiftStateRBTree);
void transplant(LiftStateRBTree*, LiftStateRBTree, LiftStateRBTree);
void delete_node(LiftStateRBTree*, LiftStateRBTree);
void delete_fixup(LiftStateRBTree*, LiftStateRBTree);
void delete_data(LiftStateRBTree*, LiftState, const LiftGrid_t * const);
LiftStateRBTree find_minimum_node(LiftStateRBTree*);
LiftStateRBTree find_maximum_node(LiftStateRBTree*);
LiftState find_minimum(LiftStateRBTree*);
LiftState find_maximum(LiftStateRBTree*);
LiftStateRBTree find_node(LiftStateRBTree*, LiftState, const LiftGrid_t * const);
int find_tag(LiftStateRBTree*, LiftState, const LiftGrid_t * const);
int is_member(LiftStateRBTree*, LiftState, const LiftGrid_t * const);
void free_lift_state_rbtree(LiftStateRBTree*, const LiftGrid_t * const);

void s_left_rotate(StateRBTree*, StateRBTree);
void s_right_rotate(StateRBTree*, StateRBTree);
void s_insert_data(StateRBTree*, State, const Grid_t * const);
void s_insert_tagged_data(StateRBTree*, State, int, const Grid_t * const);
void s_insert_node(StateRBTree*, StateRBTree, const Grid_t * const);
void s_insert_fixup(StateRBTree*, StateRBTree);
void s_transplant(StateRBTree*, StateRBTree, StateRBTree);
void s_delete_node(StateRBTree*, StateRBTree);
void s_delete_fixup(StateRBTree*, StateRBTree);
void s_delete_data(StateRBTree*, State, const Grid_t * const);
StateRBTree s_find_minimum_node(StateRBTree*);
StateRBTree s_find_maximum_node(StateRBTree*);
State s_find_minimum(StateRBTree*);
State s_find_maximum(StateRBTree*);
StateRBTree s_find_node(StateRBTree*, State, const Grid_t * const);
int s_find_tag(StateRBTree*, State, const Grid_t * const);
int s_is_member(StateRBTree*, State, const Grid_t * const);
void free_state_rbtree(StateRBTree*);

#endif
