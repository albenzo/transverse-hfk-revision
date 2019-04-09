#include "states.h"

/**
 * Allocates memory for the supplied lift state
 * @param s a pointer to a lift state
 * @param G a pointer to a lift grid
 */
void init_lift_state(LiftState *s, const LiftGrid_t * const G) {
  *s = (char**)malloc(sizeof(char*)*G->sheets);
  for(int i = 0; i < G->sheets; ++i) {
    (*s)[i] = malloc(sizeof(char)*G->arc_index);
  }
}

void copy_lift_state(LiftState* dest, const LiftState * const origin, const LiftGrid_t * const G) {
  for(int i = 0; i < G->sheets; ++i) {
    for(int j = 0; j < G->arc_index; ++j) {
      (*dest)[i][j] = (*origin)[i][j];
    }
  }
}

/**
 * Frees the passed in lift state
 * @param s a pointer to a lift state
 * @param g a pointer to a grid
 */
void free_lift_state(LiftState *s, const LiftGrid_t * const G) {
  for(int i = 0; i < G->sheets; ++i) {
    free((*s)[i]);
  }
  free(*s);  
}

/**
 * Determines whether the supplied state is a valid grid state for the supplied
 * grid
 * @param state the grid state
 * @param G the grid
 * @return 1 if the state is valid, 0 otherwise
 */
int is_state(const State state, const Grid_t *const G) {
  for (int i = 0; i < G->arc_index; ++i) {
    if (state[i] <= 0 || state[i] > G->arc_index) {
      return 0;
    }

    int found_i = 0;

    for (int j = 0; j < G->arc_index; ++j) {
      if (state[j] == i + 1) {
        found_i = 1;
        break;
      }
    }

    if (!found_i) {
      return 0;
    }
  }
  return 1;
}

/**
 * Determines whether the supplied lift state is valid for the supplied
 * lift grid
 * @param state a lift state
 * @param G a lift grid
 * @return 1 if the state is valid, 0 otherwise
 */
int is_lift_state(const LiftState state, const LiftGrid_t * const G) {
  return 1;
}


/**
 * Determines whether two states represent the same value
 * @param a a state
 * @param b a state
 * @param G a Grid
 */
int eq_state(const State a, const State b, const Grid_t *const G) {
  return (!strncmp(a, b, G->arc_index));
}


/**
 * Determines whether the two lift states represent the same state
 * @param a a lift state
 * @param b a lift state
 * @param G a lift grid
 * @return 1 if the states are equivalent, 0 otherwise.
 */
int eq_lift_state(const LiftState a, const LiftState b, const LiftGrid_t *const G) {
  for(int i=0; i < G->sheets; ++i) {
    if(0 != memcmp(a[i],b[i], G->arc_index)) {
      return 0;
    }
  }
  return 1;
}

/**
 * Compares two lift states with dictionary ordering
 * @param u a lift state
 * @param v a lit state
 * @param G a lift grid
 * @return 0 if the states are equivalent, -1 if u<v, and 1 if u>v
 */
int comp_lift_state(const LiftState u, const LiftState v, const LiftGrid_t * const G) {
  for(int i=0; i< G->sheets; ++i) {
    int comp = memcmp(u[i],v[i], G->arc_index);
    if (0 != comp) {
      return comp;
    }
  }
  return 0;
}

/**
 * Modifies the supplied lift state to be its mirror
 * @param state a lift state
 * @param G a grid
 */
void mirror_lift_state(LiftState * state, const LiftGrid_t * const G) {
  LiftState original;
  init_lift_state(&original, G);
  copy_lift_state(&original, state, G);
  for(int i = 0; i < G->sheets; ++i) {
    (*state)[i][0] = original[G->sheets-(i+1)][0];
    for(int j = 0; j < G->arc_index-1; ++j) {
      (*state)[i][j+1] = original[G->sheets-(i+1)][G->arc_index-(j+1)];
    }
  }
  free_lift_state(&original, G);
}

/**
 * Determines whether the supplied grid is valid
 * @param G a grid
 * @return 1 if the grid is valid, 0 otherwise
 */
int is_grid(const Grid_t *const G) {
  if (1 >= G->arc_index) {
    return 0;
  }

  for (int j = 0; j < G->arc_index; ++j) {
    if (G->Xs[j] == G->Os[j]) {
      return 0;
    }
    int found_x = 0;
    int found_o = 0;

    for (int k = 0; k < G->arc_index && (found_x == 0 || found_o == 0); ++k) {
      if (G->Xs[k] == j + 1) {
        found_x = 1;
      }

      if (G->Os[k] == j + 1) {
        found_o = 1;
      }
    }

    if (0 == found_x || 0 == found_o) {
      return 0;
    }
  }

  return 1;
}

/**
 * Determines whether the supplied lift grid is valid
 * @param G a lift grid
 * @return 1 if the grid is valid, 0 otherwise
 */
int is_lift_grid(const LiftGrid_t *const G) {
  if(G->sheets < 1) {
    return 0;
  }
  
  Grid_t H;
  H.arc_index = G->arc_index;
  H.Xs = G->Xs;
  H.Os = G->Os;
  return is_grid(&H);
}

/**
 * Allocates a new grid that is the mirror of the supplied grid
 * @param G a lift grid
 * @param the mirror of G
 */
LiftGrid_t* mirror_lift_grid(const LiftGrid_t * const G) {
  LiftGrid_t* G_mirror = malloc(sizeof(LiftGrid_t));
  G_mirror->Xs = malloc(sizeof(char)*G->arc_index);
  G_mirror->Os = malloc(sizeof(char)*G->arc_index);
  G_mirror->arc_index = G->arc_index;
  G_mirror->sheets = G->sheets;

  for(int i=0; i<G->arc_index; ++i) {
    G_mirror->Xs[i] = G->Xs[G->arc_index-(i+1)];
    G_mirror->Os[i] = G->Os[G->arc_index-(i+1)];
  }
  return G_mirror;
}

/**
 * Computes the writhe of the passed grid
 * @param G working grid
 * @return writhe of the grid
 */
int get_writhe(const Grid_t *const G) {
  int i = 1, j = 0, k = 0;
  int max_XO, min_XO;
  int writhe = 0;
  int temp_X, temp_O;
  int current_X, current_O;
  while (i < G->arc_index) {
    temp_X = G->Xs[i] - '0';
    temp_O = G->Os[i] - '0';
    min_XO = temp_X < temp_O ? temp_X : temp_O;
    if (min_XO == temp_X)
      max_XO = temp_O;
    else
      max_XO = temp_X;
    j = 0;
    while (j < i) {
      current_X = G->Xs[j] - '0';
      current_O = G->Os[j] - '0';
      if (min_XO < current_X && max_XO > current_X) {
        k = i + 1;
        while (k < G->arc_index) {
          if (G->Os[k] - '0' == current_X) {
            if (max_XO == temp_X)
              writhe++;
            else
              writhe--;
          }
          k++;
        }
      }
      if (min_XO < current_O && max_XO > current_O) {
        k = i + 1;
        while (k < G->arc_index) {
          if (G->Xs[k] - '0' == current_O) {
            if (max_XO == temp_O)
              writhe++;
            else
              writhe--;
          }
          k++;
        }
      }
      j++;
    }
    i++;
  }
  return writhe;
}

/**
 * Computes up_down_cusps array, which stores the number of up cusps in first
 * position and number of down cusps in second
 * @param up_down_cusps a length two int array
 * @param G working grid
 * @return number of up cusps stored in the first param and down cusps in the
 * second
 */
void cusps(int *up_down_cusps, const Grid_t *const G) {
  int i = 0, j = 0;
  while (i < G->arc_index) {
    if (G->Xs[i] < G->Os[i]) {
      j = i + 1;
      while (j < G->arc_index) {
        if (G->Os[j] == G->Xs[i])
          up_down_cusps[0]++;
        j++;
      }
    }
    if (G->Os[i] < G->Xs[i]) {
      j = i + 1;
      while (j < G->arc_index) {
        if (G->Xs[j] == G->Os[i])
          up_down_cusps[1]++;
        j++;
      }
    }

    if (G->Xs[i] > G->Os[i]) {
      j = 0;
      while (j < i) {
        if (G->Os[j] == G->Xs[i])
          up_down_cusps[1]++;
        j++;
      }
    }
    if (G->Os[i] > G->Xs[i]) {
      j = 0;
      while (j < i) {
        if (G->Xs[j] == G->Os[i])
          up_down_cusps[0]++;
        j++;
      }
    }

    i++;
  }
}


/**
 * Sum over each point in the permutation count the number of Os
 * that occur to the northeast
 * @param x a permutation
 * @param G working Grid
 * @return an int containing the quantity described above
 */
int NESW_pO(const State x, const Grid_t *const G) {
  int i = 0, j = 0;
  int ans = 0;
  while (i < G->arc_index) {
    j = i;
    while (j < G->arc_index) {
      if (x[i] <= G->Os[j]) {
        ans++;
      };
      j++;
    };
    i++;
  };
  return (ans);
}

/**
 * Sum over each O in Os count the number of points in the permutation
 * to the northeast
 * @param x a permutation
 * @param G working Grid
 * @return an int containing the quantity described above
 */
int NESW_Op(const State x, const Grid_t *const G) {
  int i = 0, j = 0;
  int ans = 0;
  while (i < G->arc_index) {
    j = i + 1;
    while (j < G->arc_index) {
      if (G->Os[i] < x[j]) {
        ans++;
      };
      j++;
    };
    i++;
  };
  return (ans);
}

/**
 * Sum over each point in the permutation count the number of points in
 * the same permutation that occur to the northeast
 * @param x a permutation
 * @param G working grid
 * @return an int containing the quantity described above
 */
int NESW_pp(const State x, const Grid_t *const G) {
  int i = 0, j = 0;
  int ans = 0;
  while (i < G->arc_index) {
    j = i;
    while (j < G->arc_index) {
      if (x[i] < x[j]) {
        ans++;
      };
      j++;
    };
    i++;
  };
  return (ans);
}

/**
 * Returns the index of a State within a StateList
 * Note: Indexed from 1
 * @param a State
 * @param b a StateList
 * @return index of a within b
 */
int get_number(const State a, const StateList b, const Grid_t *const G) {
  StateList temp;
  int count = 1;
  temp = b;
  while (temp != NULL) {
    if (eq_state(a, temp->data, G)) {
      return count;
    }
    temp = temp->nextState;
    count++;
  }
  return 0;
}

/**
 * Returns the index of a LiftState within a StateList
 * Note: Indexed from 1
 * @param a LiftState
 * @param b a LiftStateList
 * @return index of a within b
 */
int get_lift_number(const LiftState a, const LiftStateList b, const LiftGrid_t *const G) {
  LiftStateList temp;
  int count = 1;
  temp = b;
  while (temp != NULL) {
    if (eq_lift_state(a, temp->data, G)) {
      return count;
    }
    temp = temp->nextState;
    ++count;
  }
  return 0;
}

/**
 * Removes a state from a StateList. Equality checked using eq_state.
 * @param a a State
 * @param v a StateList
 * @return a StateList removing a from v
 * @see eq_state
 */
StateList remove_state(const State a, const StateList v,
                       const Grid_t *const G) {
  StateList temp, prev;
  StateList s_list = v;
  prev = v;
  if (v == NULL)
    return (NULL);
  else if (eq_state(a, v->data, G)) {
    temp = v;
    s_list = v->nextState;
    free(v->data);
    free(v);
    return (s_list);
  } else {
    temp = prev->nextState;
    while ((temp != NULL) && (!eq_state(a, temp->data, G))) {
      temp = temp->nextState;
      prev = prev->nextState;
    };
    if (temp != NULL) {
      prev->nextState = temp->nextState;
      free(temp->data);
      free(temp);
      return s_list;
    } else
      return s_list;
  }
}

/**
 * Remove a lift state from the supplied lift state list. Equality checked
 * using eq_lift_state
 * @param a a lift state
 * @param v a lift state list
 * @param G a lift grid
 */
LiftStateList remove_lift_state(const LiftState a, const LiftStateList v,
                       const LiftGrid_t *const G) {
  LiftStateList temp, prev;
  LiftStateList s_list = v;
  prev = v;
  if (v == NULL)
    return (NULL);
  else if (eq_lift_state(a, v->data, G)) {
    temp = v;
    s_list = v->nextState;
    free_lift_state(&v->data, G);
    free(v);
    return (s_list);
  } else {
    temp = prev->nextState;
    while ((temp != NULL) && (!eq_lift_state(a, temp->data, G))) {
      temp = temp->nextState;
      prev = prev->nextState;
    }
    if (temp != NULL) {
      prev->nextState = temp->nextState;
      free_lift_state(&temp->data, G);
      free(temp);
      return s_list;
    } else
      return s_list;
  }
}

/**
 * Frees the supplied StateList
 * @param states a StateList
 */
void free_state_list(StateList states) {
  StateList temp;
  temp = states;
  while (temp != NULL) {
    states = states->nextState;
    free(temp->data);
    free(temp);
    temp = states;
  }
}

/**
 * Frees the supplied lift state list and the data contained within
 * @param states a lift state list
 * @param G a lift grid
 */
void free_lift_state_list(LiftStateList states, const LiftGrid_t * const G) {
  LiftStateList temp;
  temp = states;
  while (temp != NULL) {
    states = states-> nextState;
    free_lift_state(&temp->data, G);
    free(temp);
    temp = states;
  }
}

/**
 * Left rotates the node x as a descendent of root
 * @param root the root of the RBTree
 * @param x a descendent node of root
 */
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

/**
 * Right rotates the node x as a descendant of root.
 * @param root the root of the RBTree
 * @param x a descendant node of root
 */
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

/**
 * Takes a lift state and inserts it into the supplied RBTree
 * @param root an RBTree
 * @param s a lift state
 * @param G a lift grid
 */
void insert_data(LiftStateRBTree root, LiftState s, const LiftGrid_t * const G) {
  LiftStateRBTree node = malloc(sizeof(LiftStateRBTreeNode_t));
  node->data = s;
  insert_node(root, node, G);
}

/**
 * Takes a node z and inserts it into the RBTree root
 * @param root the root of the RBTree
 * @param z the node to be inserted
 * @param G a lift grid
 */
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

/**
 * Repairs the RBTree after an insertion
 * @param root the root of the RBTree
 * @param z the node which is to be repaired about
 */
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

/**
 * Transplants two nodes within root
 * @param root the root of the RBTree
 * @param u a descendant node of root
 * @param v a descendant node of root
 */
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

/**
 * Deletes the supplied node from root.
 * @param root the root of the RBTree
 * @param z the node to be deleted
 */
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

/**
 * Repairs root after a delete
 * @param root the root of the RBTree
 * @param x the node to be repaired about
 */
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

/**
 * Deletes the supplied lift state once from root.
 * @param root the root of the RBTree
 * @param s the lift state to be deleted
 * @param G a lift grid
 */
void delete_data(LiftStateRBTree root, LiftState s, const LiftGrid_t * const G) {
  delete_node(root, find_node(root, s, G));
}

/**
 * Finds and returns the minimal node in root
 * @param root the root of the RBTree
 * @return the node containing the minimum value of root
 */
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

/**
 * Finds and returns the maximal node in root
 * @param root the root of the RBTree
 * @return the node containing the maximum value of root
 */
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

/**
 * Finds and returns the minimum state in root
 * @param root the root of the RBTree
 * @return the minimal lift state in root
 */
LiftState find_minimum(LiftStateRBTree root) {
  LiftStateRBTree min_node = find_minimum_node(root);
  if (NULL == min_node) {
    return NULL;
  }
  return min_node->data;
}

/**
 * Finds and returns the maximum lift state in root
 * @param root the root of the RBTree
 * @return the maximal lift state in root
 */
LiftState find_maximum(LiftStateRBTree root) {
  LiftStateRBTree max_node = find_maximum_node(root);
  if (NULL == max_node) {
    return NULL;
  }
  return max_node->data;
}

/**
 * Finds and returns a node corresponding to the supplied lift state in root
 * @param root the root of the RBTree
 * @param s the lift state to be found
 * @param G a lift grid
 * @return the node containing s, returns NULL if not found
 */
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

/**
 * Determines whether a lift state is contained in root
 * @param root the root of the RBTree
 * @param s the lift state to be tested
 * @param G a lift grid
 * @return 1 if s is contained in some node of root, 0 otherwise
 */
int is_member(LiftStateRBTree root, LiftState s, const LiftGrid_t * const G) {
  return find_node(root, s, G) != NULL;
}

/**
 * Recursively frees root, all descendants, and their contained data
 * @param root the root of the RBTree
 * @param G a lift grid
 */
void free_lift_state_rbtree(LiftStateRBTree root, const LiftGrid_t * const G) {
  if(NULL == root) {
    return;
  }
  free_lift_state_rbtree(root->left, G);
  free_lift_state_rbtree(root->right, G);
  free_lift_state(&root->data, G);
  free(root);
}
