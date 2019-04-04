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
  LiftState data;
  LiftStateRBTree left;
  LiftStateRBTree right;
  LiftStateRBTree parent;
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

void init_lift_state(LiftState*, const LiftGrid_t * const G);
void free_lift_state(LiftState*, const LiftGrid_t * const G);
int eq_lift_state(const LiftState a, const LiftState b, const LiftGrid_t * const G);
int comp_lift_state(const LiftState, const LiftState, const LiftGrid_t * const G);
void left_rotate(LiftStateRBTree, LiftStateRBTree);
void right_rotate(LiftStateRBTree, LiftStateRBTree);
void insert_data(LiftStateRBTree, LiftState, const LiftGrid_t * const);
void insert_node(LiftStateRBTree, LiftStateRBTree, const LiftGrid_t * const);
void insert_fixup(LiftStateRBTree, LiftStateRBTree);
void transplant(LiftStateRBTree, LiftStateRBTree, LiftStateRBTree);
void delete_node(LiftStateRBTree, LiftStateRBTree);
void delete_fixup(LiftStateRBTree, LiftStateRBTree);
void delete_data(LiftStateRBTree, LiftState, const LiftGrid_t * const);
LiftStateRBTree find_minimum_node(LiftStateRBTree);
LiftStateRBTree find_maximum_node(LiftStateRBTree);
LiftState find_minimum(LiftStateRBTree);
LiftState find_maximum(LiftStateRBTree);
LiftStateRBTree find_node(LiftStateRBTree, LiftState, const LiftGrid_t * const);
int is_member(LiftStateRBTree, LiftState, const LiftGrid_t * const);

void init_lift_state(LiftState *s, const LiftGrid_t * const G) {
  *s = (char**)malloc(sizeof(char*)*G->sheets);
  for(int i = 0; i < G->sheets; ++i) {
    (*s)[i] = malloc(sizeof(char)*G->arc_index);
  }
}

void free_lift_state(LiftState *s, const LiftGrid_t * const G) {
  for(int i = 0; i < G->sheets; ++i) {
    free((*s)[i]);
  }
  free(*s);  
}

int eq_lift_state(const LiftState a, const LiftState b, const LiftGrid_t *const G) {
  for(int i=0; i < G->sheets; ++i) {
    if(!strncmp(a[i],b[i], G->arc_index)) {
      return 0;
    }
  }
  return 1;
}

int comp_lift_state(const LiftState u, const LiftState v, const LiftGrid_t * const G) {
  for(int i=0; i< G->sheets; ++i) {
    int comp = strncmp(u[i],v[i], G->arc_index);
    if (0 != comp) {
      return comp;
    }
  }
  return 0;
}

void left_rotate(LiftStateRBTree root, LiftStateRBTree x) {
  LiftStateRBTree y = x->right;
  x->right = y->left;
  if (NULL != y->left) {
    y->left->parent = x;
  }
  y->parent = x->parent;
  if (NULL == x->parent) {
    *root = *y;
  }
  else if (x == x->parent->left) {
    x->parent->left = y;
  }
  else {
    x->parent->right = y;
  }
  y->left = x;
  x->parent = y;
}

void right_rotate(LiftStateRBTree root, LiftStateRBTree x) {
  LiftStateRBTree y = x->left;
  x->left = y->right;
  if (NULL != y->right) {
    y->right->parent = x;
  }
  y->parent = x->parent;
  if (NULL == x->parent) {
    *root = *y;
  }
  else if (x == x->parent->right) {
    x->parent->right = y;
  }
  else {
    x->parent->left = y;
  }
  y->right = x;
  x->parent = y;
}

void insert_data(LiftStateRBTree root, LiftState s, const LiftGrid_t * const G) {
  LiftStateRBTree node = malloc(sizeof(LiftStateRBTreeNode_t));
  node->data = s;
  insert_node(root, node, G);
}

void insert_node(LiftStateRBTree root, LiftStateRBTree z, const LiftGrid_t * const G) {
  LiftStateRBTree y = NULL;
  LiftStateRBTree x = root;
  while (NULL != x) {
    y = x;
    if (comp_lift_state(z->data, x->data, G) < 0) {
      x = x->left;
    }
    else {
      x = x->right;
    }
  }
  z->parent = y;
  if (NULL == y) {
    *root = *z;
  }
  else if (comp_lift_state(z->data, y->data, G) < 0) {
    y->left = z;
  }
  else {
    y->right = z;
  }
  z->left = NULL;
  z->right = NULL;
  z->color = RED;
  insert_fixup(root, z);
}

void insert_fixup(LiftStateRBTree root, LiftStateRBTree z) {
  while (RED == z->parent->color) {
    if (z->parent == z->parent->parent->left) {
      LiftStateRBTree y = z->parent->parent->right;
      if (RED == y->color) {
        z->parent->color = BLACK;
        y->color = BLACK;
        z->parent->parent->color = RED;
        z = z->parent->parent;
      }
      else if (z == z->parent->right) {
        z = z->parent;
        left_rotate(root, z);
      }
      else {
        z->parent->color = BLACK;
        z->parent->parent->color = RED;
        right_rotate(root, z->parent->parent);
      }
    }
    else {
      LiftStateRBTree y = z->parent->parent->left;
      if (RED == y->color) {
        z->parent->color = BLACK;
        y->color = BLACK;
        z->parent->parent->color = RED;
        z = z->parent->parent;
      }
      else if (z == z->parent->left) {
        z = z->parent;
        right_rotate(root, z);
      }
      else {
        z->parent->color = BLACK;
        z->parent->parent->color = RED;
        left_rotate(root, z->parent->parent);
      }
    }
  }
}

void transplant(LiftStateRBTree root, LiftStateRBTree u, LiftStateRBTree v) {
  if (NULL == u->parent) {
    *root = *v;
  }
  else if (u == u->parent->left) {
    u->parent->left = v;
  }
  else {
    u->parent->right = v;
  }
  v->parent = u->parent;
}

void delete_node(LiftStateRBTree root, LiftStateRBTree z) {
  if(NULL == root || NULL == z) {
    return;
  }
  
  LiftStateRBTree y = z;
  LiftStateRBTree x;
  int y_original_color = y->color;

  if (NULL == z->left) {
    x = z->right;
    transplant(root, z, z->right);
  }
  else if (NULL == z->right) {
    x = z->left;
    transplant(root, z, z->left);
  }
  else {
    y = find_minimum_node(z->right);
    y_original_color = y->color;
    x = y->right;
    if (y->parent == z) {
      x->parent = y;
    }
    else {
      transplant(root, y, y->right);
      y->right = z->right;
      y->right->parent = y;
    }
    transplant(root, z, y);
    y->left = z->left;
    y->left->parent = y;
    y->color = z->color;
  }
  if (BLACK == y_original_color) {
    delete_fixup(root, x);
  }
}

void delete_fixup(LiftStateRBTree root, LiftStateRBTree x) {
  while (x != root && x->color == BLACK) {
    LiftStateRBTree w;
    if (x == x->parent->left) {
      w = x->parent->right;
      if (RED == w->color) {
        w->color = BLACK;
        x->parent->color = RED;
        left_rotate(root, x->parent);
        w = x->parent->right;
      }
      if (BLACK == w->left->color && BLACK == w->right->color) {
        w->color = RED;
        x = x->parent;
      }
      else {
        if (BLACK == w->right->color) {
          w->left->color = BLACK;
          w->color = RED;
          right_rotate(root, w);
          w = x->parent->right;
        }
        w->color = x->parent->color;
        x->parent->color = BLACK;
        w->right->color = BLACK;
        left_rotate(root, x->parent);
        x = root;
      }
    }
    else {
      w = x->parent->left;
      if (RED == w->color) {
        w->color = BLACK;
        x->parent->color = RED;
        right_rotate(root, x->parent);
        w = x->parent->left;
      }
      if (BLACK == w->left->color && BLACK == w->right->color) {
        w->color = RED;
        x = x->parent;
      }
      else {
        if (BLACK == w->left->color) {
          w->right->color = BLACK;
          w->color = RED;
          left_rotate(root, w);
          w = x->parent->left;
        }
        w->color = x->parent->color;
        x->parent->color = BLACK;
        w->left->color = BLACK;
        right_rotate(root, x->parent);
        x = root;
      }
    }
  }
  x->color = BLACK;
}

void delete_data(LiftStateRBTree root, LiftState s, const LiftGrid_t * const G) {
  delete_node(root, find_node(root, s, G));
}

LiftStateRBTree find_minimum_node(LiftStateRBTree root) {
  if (NULL == root) {
    return NULL;
  }
  LiftStateRBTree iter = root;
  while (iter->left != NULL) {
    iter = iter->left;
  }
  return iter;
}

LiftStateRBTree find_maximum_node(LiftStateRBTree root) {
  if (NULL == root) {
    return NULL;
  }
  LiftStateRBTree iter = root;
  while (iter->right != NULL) {
    iter = iter->right;
  }
  return iter;
}

LiftState find_minimum(LiftStateRBTree root) {
  LiftStateRBTree min_node = find_minimum_node(root);
  if (NULL == min_node) {
    return NULL;
  }
  return min_node->data;
}

LiftState find_maximum(LiftStateRBTree root) {
  LiftStateRBTree max_node = find_maximum_node(root);
  if (NULL == max_node) {
    return NULL;
  }
  return max_node->data;
}

LiftStateRBTree find_node(LiftStateRBTree root, LiftState s, const LiftGrid_t * const G) {
  LiftStateRBTree iter = root;
  while (iter != NULL) {
    int comp = comp_lift_state(s, iter->data, G);
    if (0 == comp) {
      return iter;
    }
    else if (comp < 0) {
      iter = iter->left;
    }
    else {
      iter = iter->right;
    }
  }

  return NULL;
}

int is_member(LiftStateRBTree root, LiftState s, const LiftGrid_t * const G) {
  return find_node(root, s, G) != NULL;
}

void free_lift_state_rbtree(LiftStateRBTree root, const LiftGrid_t * const G) {
  if(NULL == root) {
    return;
  }
  free_lift_state_rbtree(root->left, G);
  free_lift_state_rbtree(root->right, G);
  free_lift_state(&root->data, G);
  free(root);
}
