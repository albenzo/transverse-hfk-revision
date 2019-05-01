/* Program for calculating the invariant for transverse knots, to
   accompany "Transverse knots distinguished by Knot Floer Homology"
   by L. Ng, P. S. Ozsvath, and D. P. Thurston. The invariant was defined
   in "Legendrian knots, transverse knots, and combinatorial Floer homology"
   by P. S. Ozsvath, Z. Szabo, and D. P. Thurston*/

/* In progress cleanup of the above listed file by
 * Lucas Meyers <lmeye22@lsu.edu>
 */

#include "TransverseHFK.h"

printf_t print_ptr = printf;
static int verbosity = SILENT;
static LiftStateRBTree new_lift_rectangles_out_internal(const LiftStateRBTree, const LiftState, const LiftGrid_t * const, int);

/**
 * Sets the print function to the passed in function pointer
 * @param print_fn a function pointer to a print function with the
 * same type as printf
 */
void set_print_fn(printf_t print_fn) {
  print_ptr = print_fn;
}

/**
 * Returns the current level of verbosity
 * @return the verbosity
 */
int get_verbosity(){
  return verbosity;
}

/**
 * Set the level of verbosity to SILENT, QUIET, or VERBOSE
 * @param val the verbosity level
 */
void set_verbosity(const int val) {
  verbosity = val;
}

/**
 * Shifts the input towards the interval [0,arc_index) by
 * a multiple of arc_index
 * @param a An integer
 * @param arc_index an int
 * @return a shifted towards the interval [0,arc_index)
 */
int mod(const int a, const int arc_index) {
  if (a >= arc_index) {
    return (a - arc_index);
  } else if (a < 0) {
    return (a + arc_index);
  } else {
    return (a);
  };
}

/**
 * Shifts the input towards the interval (0,arc_index] by
 * a multiple of arc_index
 * @param a An integer
 * @param arc_index an int
 * @return a shifted towards the interval (0,arc_index]
 */
int mod_up(const int a, const int arc_index) {
  if (a > arc_index) {
    return (a - arc_index);
  } else if (a <= 0) {
    return (a + arc_index);
  } else {
    return (a);
  };
}

/**
 * @param a an int
 * @param b an int
 * @return the minimum of a and b
 */
int min(const int a, const int b) {
  if (a < b) {
    return (a);
  } else {
    return (b);
  }
}


/**
 * Creates a new state equivalent to incoming with indices x1 and x2 swapped
 * @param x1 a non negative int below the arc index of G
 * @param x2 a non negative int below the arc index of G
 * @param incoming a permutation
 * @param G a grid
 * @return a new state equivalent to incoming with indices x1 and x2 swapped
 */
State swap_cols(const int x1, const int x2, const State incoming, const Grid_t *const G) {
  State ans = malloc(sizeof(char)*G->arc_index);

  for(int i=0; i< G->arc_index; ++i) {
    ans[i] = incoming[i];
  }
  ans[x1] = incoming[x2];
  ans[x2] = incoming[x1];

  return ans;
}

/**
 * Returns a single element statelist containing
 * the permutation incoming with x1 and x2 swapped
 * @param x1 a non negative int below the arc index of G
 * @param x2 a non negative int below the arc index of G
 * @param incoming a permutation
 * @param G working grid
 * @return a single element StateList containing a copy of
 * the permutation incoming with entries x1 and x2 swapped
 */
StateList swap_cols_list(const int x1, const int x2, const State incoming,
                    const Grid_t *const G) {
  StateList ans;
  int i;
  i = 0;
  ans = malloc(sizeof(StateNode_t));
  ans->data = malloc(sizeof(char) * G->arc_index);

  while (i < G->arc_index) {
    ans->data[i] = incoming[i];
    i++;
  };
  ans->data[x1] = (incoming)[x2];
  ans->data[x2] = (incoming)[x1];
  ans->nextState = NULL;
  return ans;
}

/**
 * Calculates whether the supplied state is nullhomologous
 * @param init a State
 * @param G working grid
 * @return nonzero if nullhomologous and zero otherwise.
 */
int null_homologous_D0Q(const State init, const Grid_t *const G) {
  StateRBTree new_ins, new_outs;
  StateRBTree prev_ins, prev_outs;
  StateRBTree potential_outs = EMPTY_TREE, potential_ins = EMPTY_TREE;
  int ans, prev_in_number, total_in, total_out;
  int edge_count = 0;
  int num_ins = 0;
  int num_outs = 0;
  int num_new_ins = 0;
  int num_new_outs = 0;
  prev_outs = EMPTY_TREE;
  prev_ins = EMPTY_TREE;
  new_ins = EMPTY_TREE;

  State s = malloc(sizeof(char) * G->arc_index);
  copy_state(&s, &init, G);

  // Create sentinal edge from A_0
  s_insert_tagged_data(&new_ins, s, 1, G);
  EdgeList edge_list = prepend_edge(0, 1, NULL);
  
  ans = 0;
  int current_pos = 1;
  while (new_ins != EMPTY_TREE && !ans) {
    num_new_outs = 0;
    total_in = 0;
    new_outs = EMPTY_TREE;
    EdgeList new_edges = NULL;
    if (get_verbosity() >= VERBOSE) {
      (*print_ptr)("Gathering A_%d:\n", current_pos);
    }

    // Build A_i by looking for states into B_(i-1) that are not in A_(i-1)
    StateTreeIter_t * present_iter;
    for(present_iter = s_create_iter(new_ins); s_has_next(present_iter);) {
      StateRBTree present_in = s_get_next(present_iter);
      free_state_rbtree(&potential_outs);
      total_in++;
      potential_outs = new_rectangles_into(prev_outs, present_in->data, G);

      StateTreeIter_t * potential_iter;
      for(potential_iter = s_create_iter(potential_outs); s_has_next(potential_iter);) {
        StateRBTree potential_out = s_get_next(potential_iter);
        StateRBTree node = s_find_node(&new_outs, potential_out->data, G);
        if(EMPTY_TREE == node) {
          State t = malloc(sizeof(char)*G->arc_index);
          copy_state(&t,&(potential_out->data),G);
          num_new_outs++;
          s_insert_tagged_data(&new_outs, t, num_new_outs, G);
          new_edges = prepend_edge(num_new_outs + num_outs, present_in->tag + num_ins, new_edges);
        }
        else {
          new_edges = prepend_edge(node->tag + num_outs, present_in->tag + num_ins, new_edges);
        }
        edge_count++;
      }
      s_free_iter(potential_iter);
    }
    s_free_iter(present_iter);
    
    if (get_verbosity() >= VERBOSE) {
      print_edges(new_edges);
      (*print_ptr)("\n");
    }
    free_state_rbtree(&prev_ins);
    prev_ins = new_ins;
    num_ins = num_ins + total_in;
    prev_in_number = num_ins;
    num_new_ins = 0;
    new_ins = EMPTY_TREE;
    total_out = 0;
    if (get_verbosity() >= VERBOSE) {
      (*print_ptr)("Gathering B_%d:\n", current_pos);
    }

    // Build B_i by finding states out of A_i that are not in B_(i-1)
    for(present_iter = s_create_iter(new_outs); s_has_next(present_iter);) {
      StateRBTree present_out = s_get_next(present_iter);
      free_state_rbtree(&potential_ins);
      total_out++;
      potential_ins = new_rectangles_out_of(prev_ins, present_out->data, G);

      StateTreeIter_t * potential_iter;
      for(potential_iter = s_create_iter(potential_ins); s_has_next(potential_iter);) {
        StateRBTree potential_in = s_get_next(potential_iter);
        StateRBTree node = s_find_node(&new_ins, potential_in->data, G);
        if(EMPTY_TREE == node) {
          State t = malloc(sizeof(char)*G->arc_index);
          copy_state(&t,&(potential_in->data),G);
          num_new_ins++;
          s_insert_tagged_data(&new_ins, t, num_new_ins, G);
          new_edges = prepend_edge(present_out->tag + num_outs, num_new_ins + num_ins, new_edges);
        }
        else {
          new_edges = prepend_edge(present_out->tag + num_outs, node->tag + num_ins, new_edges);
        }
        edge_count++;
      }
      s_free_iter(potential_iter);
    }
    s_free_iter(present_iter);
    
    if (get_verbosity() >= VERBOSE) {
      print_edges(new_edges);
      (*print_ptr)("\n");
    }
    free_state_rbtree(&prev_outs);
    prev_outs = new_outs;
    new_outs = EMPTY_TREE;

    new_edges = merge_sort_edges(new_edges);
    edge_list = merge_edges(edge_list, new_edges);

    if (get_verbosity() >= VERBOSE) {
      (*print_ptr)("Full edge list:\n");
      print_edges(edge_list);
    }

    if (get_verbosity() >= VERBOSE) {
      (*print_ptr)("Contracting edges from 0 to %d:\n", prev_in_number);
    }
    
    special_homology(0, prev_in_number, &edge_list);
    if (get_verbosity() >= VERBOSE) {
      print_edges(edge_list);
      (*print_ptr)("\n");
    }
    
    if ((edge_list == NULL) || (edge_list->start != 0)) {
      // If there are no edges out of A_0 (sentinal is gone) after contraction init is null-homologous
      ans = 1;
      if (get_verbosity() >= VERBOSE) {
        (*print_ptr)("No edges pointing out of A_0!\n");
      }
      free_state_rbtree(&new_ins);
      free_state_rbtree(&new_outs);
    } else if (edge_list->end <= prev_in_number) {
      // If edges out of A_0 cannot be removed anymore (sentinal will never vanish) init is not null-homologous
      ans = 0;
      if (get_verbosity() >= VERBOSE) {
        (*print_ptr)("There exist edges pointing from A_0 to B_%d! No future "
                     "contractions will remove this edge!\n",
                     current_pos - 1);
      }
      free_state_rbtree(&new_ins);
      free_state_rbtree(&new_outs);
    } else {
      num_outs = num_outs + total_out;
      if (get_verbosity() >= VERBOSE) {
        (*print_ptr)("Total number of states in B_i up to B_%d (before any "
                     "contraction): %d \n",
                     current_pos - 1, prev_in_number);
        (*print_ptr)("Total number of states in A_i up to A_%d (before any "
                     "contraction): %d \n",
                     current_pos, num_outs);
        (*print_ptr)("Total number of states in B_i up to B_%d (before any "
                     "contraction): %d \n",
                     current_pos, num_ins + total_in);
        (*print_ptr)("Total number of edges  up to A_%d and B_%d (before any "
                     "contraction): %d \n",
                     current_pos, current_pos, edge_count);
        (*print_ptr)("\n");
      }
    };
    current_pos++;
  };
  return (ans);
}

/**
 * Calculates if D1 of the supplied state is nullhomologous
 * @param init a State
 * @param G working grid
 * @return nonzero if nullhomologous and zero otherwise
 */
int null_homologous_D1Q(const State init, const Grid_t *const G) {
  StateRBTree new_ins, new_outs;
  StateRBTree prev_ins, prev_outs;
  StateRBTree potential_outs = EMPTY_TREE, potential_ins = EMPTY_TREE;
  int ans, prev_in_number, total_in, total_out;
  int edge_count = 0;
  int num_ins = 0;
  int num_outs = 0;
  int num_new_ins = 0;
  int num_new_outs = 0;
  EdgeList edge_list = prepend_edge(0, 1, NULL);
  prev_outs = EMPTY_TREE;
  prev_ins = EMPTY_TREE;
  new_ins = EMPTY_TREE;

  State s = malloc(sizeof(char) * G->arc_index);
  copy_state(&s, &init, G);

  // Calculate D1(init) and terminate if null. Otherwise build sentinal edges out of A_0
  StateList temp = fixed_wt_rectangles_out_of(1, init, G);
    
  if (NULL == temp) {
    free_state_list(temp);
    free(s);
    free_edge_list(edge_list);
    return 1;
  }

  if (temp != NULL) {
    int i = 1;
    edge_list = create_edge(0, 1);
    s_insert_tagged_data(&new_ins, temp->data, 1, G);
    temp = temp->nextState;
    while (temp != NULL) {
      i++;
      edge_list = append_ordered(0, i, edge_list);
      s_insert_tagged_data(&new_ins, temp->data, i, G);
      temp = temp->nextState;
    }
  }
  
  ans = 0;
  int current_pos = 1;
  
  while (new_ins != EMPTY_TREE && !ans) {
    num_new_outs = 0;
    total_in = 0;
    new_outs = EMPTY_TREE;
    EdgeList new_edges = NULL;
    if (get_verbosity() >= VERBOSE) {
      (*print_ptr)("Gathering A_%d:\n", current_pos);
    }

    // Build A_i by looking for states into B_(i-1) that are not in A_(i-1)
    StateTreeIter_t * present_iter;
    for(present_iter = s_create_iter(new_ins); s_has_next(present_iter);) {
      StateRBTree present_in = s_get_next(present_iter);
      free_state_rbtree(&potential_outs);
      total_in++;
      potential_outs = new_rectangles_into(prev_outs, present_in->data, G);

      StateTreeIter_t * potential_iter;
      for(potential_iter = s_create_iter(potential_outs); s_has_next(potential_iter);) {
        StateRBTree potential_out = s_get_next(potential_iter);
        StateRBTree node = s_find_node(&new_outs, potential_out->data, G);
        if(EMPTY_TREE == node) {
          State t = malloc(sizeof(char)*G->arc_index);
          copy_state(&t,&(potential_out->data),G);
          num_new_outs++;
          s_insert_tagged_data(&new_outs, t, num_new_outs, G);
          new_edges = prepend_edge(num_new_outs + num_outs, present_in->tag + num_ins, new_edges);
        }
        else {
          new_edges = prepend_edge(node->tag + num_outs, present_in->tag + num_ins, new_edges);
        }
        edge_count++;
      }
      s_free_iter(potential_iter);
    }
    s_free_iter(present_iter);
    
    if (get_verbosity() >= VERBOSE) {
      print_edges(edge_list);
      (*print_ptr)("\n");
    }
    free_state_rbtree(&prev_ins);
    prev_ins = new_ins;
    num_ins = num_ins + total_in;
    prev_in_number = num_ins;
    num_new_ins = 0;
    new_ins = EMPTY_TREE;
    total_out = 0;
    if (get_verbosity() >= VERBOSE) {
      (*print_ptr)("Gathering B_%d:\n", current_pos);
    }
    
    // Build B_i by finding states out of A_i that are not in B_(i-1)
    for(present_iter = s_create_iter(new_outs); s_has_next(present_iter);) {
      StateRBTree present_out = s_get_next(present_iter);
      free_state_rbtree(&potential_ins);
      total_out++;
      potential_ins = new_rectangles_out_of(prev_ins, present_out->data, G);

      StateTreeIter_t * potential_iter;
      for(potential_iter = s_create_iter(potential_ins); s_has_next(potential_iter);) {
        StateRBTree potential_in = s_get_next(potential_iter);
        StateRBTree node = s_find_node(&new_ins, potential_in->data, G);
        if(EMPTY_TREE == node) {
          State t = malloc(sizeof(char)*G->arc_index);
          copy_state(&t,&(potential_in->data),G);
          num_new_ins++;
          s_insert_tagged_data(&new_ins, t, num_new_ins, G);
          new_edges = prepend_edge(present_out->tag + num_outs, num_new_ins + num_ins, new_edges);
        }
        else {
          new_edges = prepend_edge(present_out->tag + num_outs, node->tag + num_ins, new_edges);
        }
        edge_count++;
      }
      s_free_iter(potential_iter);
    }
    s_free_iter(present_iter);
    
    if (get_verbosity() >= VERBOSE) {
      print_edges(edge_list);
      (*print_ptr)("\n");
    }
    free_state_rbtree(&prev_outs);
    prev_outs = new_outs;
    new_outs = EMPTY_TREE;
    if (get_verbosity() >= VERBOSE) {
      (*print_ptr)("Contracting edges from 0 to %d:\n", prev_in_number);
    }

    new_edges = merge_sort_edges(new_edges);
    edge_list = merge_edges(edge_list, new_edges);
    
    special_homology(0, prev_in_number, &edge_list);
    
    if (get_verbosity() >= VERBOSE) {
      print_edges(edge_list);
      (*print_ptr)("\n");
    }
    if ((edge_list == NULL) || (edge_list->start != 0)) {
      // If there are no edges out of A_0 (sentinal is gone) after contraction init is null-homologous
      ans = 1;
      if (get_verbosity() >= VERBOSE) {
        (*print_ptr)("No edges pointing out of A_0!\n");
      }
      free_state_rbtree(&new_ins);
      free_state_rbtree(&new_outs);
    } else if (edge_list->end <= prev_in_number) {
      // If edges out of A_0 cannot be removed anymore (sentinal will never vanish) init is not null-homologous
      ans = 0;
      if (get_verbosity() >= VERBOSE) {
        (*print_ptr)("There exist edges pointing from A_0 to B_%d! No future "
                     "contractions will remove this edge!\n",
                     current_pos - 1);
      }
      free_state_rbtree(&new_ins);
      free_state_rbtree(&new_outs);
    } else {
      num_outs = num_outs + total_out;
      if (get_verbosity() >= VERBOSE) {
        (*print_ptr)("Total number of states in B_i up to B_%d (before any "
                     "contraction): %d \n",
                     current_pos - 1, prev_in_number);
        (*print_ptr)("Total number of states in A_i up to A_%d (before any "
                     "contraction): %d \n",
                     current_pos, num_outs);
        (*print_ptr)("Total number of states in B_i up to B_%d (before any "
                     "contraction): %d \n",
                     current_pos, num_ins + total_in);
        (*print_ptr)("Total number of edges  up to A_%d and B_%d (before any "
                     "contraction): %d \n",
                     current_pos, current_pos, edge_count);
        (*print_ptr)("\n");
      }
    };
    current_pos++;
  };
  return (ans);
}

int null_homologous_lift(const LiftState init, const LiftGrid_t *const G) {
  LiftStateRBTree new_ins, new_outs;
  LiftStateRBTree prev_ins, prev_outs;
  LiftStateRBTree potential_outs = EMPTY_LIFT_TREE, potential_ins = EMPTY_LIFT_TREE;
  int ans, prev_in_number, total_in, total_out;
  int edge_count = 0;
  int num_ins = 0;
  int num_outs = 0;
  int num_new_ins = 0;
  int num_new_outs = 0;
  prev_outs = EMPTY_LIFT_TREE;
  prev_ins = EMPTY_LIFT_TREE;
  new_ins = EMPTY_LIFT_TREE;
  
  LiftState s;
  init_lift_state(&s, G);

  copy_lift_state(&s, &init, G);

  // Create sentinal edge from A_0
  insert_tagged_data(&new_ins, s, 1, G);
  EdgeList edge_list = prepend_edge(0, 1, NULL);

  ans = 0;
  int current_pos = 1;
  while (new_ins != EMPTY_LIFT_TREE && !ans) {
    num_new_outs = 0;
    total_in = 0;
    new_outs = EMPTY_LIFT_TREE;
    EdgeList new_edges = NULL;
    if (get_verbosity() >= VERBOSE) {
      (*print_ptr)("Gathering A_%d:\n", current_pos);
    }

    // Build A_i by looking for states into B_(i-1) that are not in A_(i-1) 
    LiftTreeIter_t * present_iter;
    for(present_iter = create_iter(new_ins); has_next(present_iter);) {
      LiftStateRBTree present_in = get_next(present_iter);
      free_lift_state_rbtree(&potential_outs, G);
      total_in++;
      potential_outs = new_lift_rectangles_into(prev_outs, present_in->data, G);

      LiftTreeIter_t * potential_iter;
      for(potential_iter = create_iter(potential_outs); has_next(potential_iter);) {
        LiftStateRBTree potential_out = get_next(potential_iter);
        LiftStateRBTree node = find_node(&new_outs, potential_out->data, G);
        if (EMPTY_LIFT_TREE == node) {
          LiftState t;
          init_lift_state(&t, G);
          copy_lift_state(&t,&(potential_out->data),G);
          num_new_outs++;
          insert_tagged_data(&new_outs, t, num_new_outs, G);
          new_edges = prepend_edge(num_new_outs + num_outs, present_in->tag + num_ins, new_edges);
        }
        else {
          new_edges = prepend_edge(node->tag + num_outs, present_in->tag + num_ins, new_edges);
        }
        edge_count++;
      }
      free_iter(potential_iter);
    }
    free_iter(present_iter);
    
    if (get_verbosity() >= VERBOSE) {
      print_edges(edge_list);
      (*print_ptr)("\n");
    }
    free_lift_state_rbtree(&prev_ins, G);
    prev_ins = new_ins;
    num_ins = num_ins + total_in;
    prev_in_number = num_ins;
    num_new_ins = 0;
    new_ins = EMPTY_LIFT_TREE;
    total_out = 0;
    if (get_verbosity() >= VERBOSE) {
      (*print_ptr)("Gathering B_%d:\n", current_pos);
    }

    // Build B_i by finding states out of A_i that are not in B_(i-1)
    for(present_iter = create_iter(new_outs); has_next(present_iter);) {
      LiftStateRBTree present_out = get_next(present_iter);
      free_lift_state_rbtree(&potential_ins, G);
      total_out++;
      potential_ins = new_lift_rectangles_out_of(prev_ins, present_out->data, G);

      LiftTreeIter_t * potential_iter;
      for(potential_iter = create_iter(potential_ins); has_next(potential_iter);) {
        LiftStateRBTree potential_in = get_next(potential_iter);
        LiftStateRBTree node = find_node(&new_ins, potential_in->data, G);
        if(EMPTY_LIFT_TREE == node) {
          LiftState t;
          init_lift_state(&t, G);
          copy_lift_state(&t,&(potential_in->data),G);
          num_new_ins++;
          insert_tagged_data(&new_ins, t, num_new_ins, G);
          new_edges = prepend_edge(present_out->tag + num_outs, num_new_ins + num_ins, new_edges);
        }
        else {
          new_edges = prepend_edge(present_out->tag + num_outs, node->tag + num_ins, new_edges);
        }
        edge_count++;
      }
      free_iter(potential_iter);
    }
    free_iter(present_iter);
    
    if (get_verbosity() >= VERBOSE) {
      print_edges(edge_list);
      (*print_ptr)("\n");
    }
    
    free_lift_state_rbtree(&prev_outs, G);
    prev_outs = new_outs;
    new_outs = EMPTY_LIFT_TREE;
    if (get_verbosity() >= VERBOSE) {
      (*print_ptr)("Contracting edges from 0 to %d:\n", prev_in_number);
    }

    new_edges = merge_sort_edges(new_edges);
    edge_list = merge_edges(edge_list, new_edges);
    
    special_homology(0, prev_in_number, &edge_list);
    if (get_verbosity() >= VERBOSE) {
      print_edges(edge_list);
      (*print_ptr)("\n");
    }
    if ((edge_list == NULL) || (edge_list->start != 0)) {
      // If there are no edges out of A_0 (sentinal is gone) after contraction init is null-homologous
      ans = 1;
      if (get_verbosity() >= VERBOSE) {
        (*print_ptr)("No edges pointing out of A_0!\n");
      }
      free_lift_state_rbtree(&new_ins, G);
      free_lift_state_rbtree(&new_outs, G);
      new_ins = EMPTY_LIFT_TREE;
    } else if (edge_list->end <= prev_in_number) {
      // If edges out of A_0 cannot be removed anymore (sentinal will never vanish) init is not null-homologous
      ans = 0;
      if (get_verbosity() >= VERBOSE) {
        (*print_ptr)("There exist edges pointing from A_0 to B_%d! No future "
                     "contractions will remove this edge!\n",
                     current_pos - 1);
      }
      free_lift_state_rbtree(&new_ins, G);
      free_lift_state_rbtree(&new_outs, G);
      new_ins = EMPTY_LIFT_TREE;
    } else {
      num_outs = num_outs + total_out;
      if (get_verbosity() >= VERBOSE) {
        (*print_ptr)("Total number of states in B_i up to B_%d (before any "
                     "contraction): %d \n",
                     current_pos - 1, prev_in_number);
        (*print_ptr)("Total number of states in A_i up to A_%d (before any "
                     "contraction): %d \n",
                     current_pos, num_outs);
        (*print_ptr)("Total number of states in B_i up to B_%d (before any "
                     "contraction): %d \n",
                     current_pos, num_ins + total_in);
        (*print_ptr)("Total number of edges  up to A_%d and B_%d (before any "
                     "contraction): %d \n",
                     current_pos, current_pos, edge_count);
        (*print_ptr)("\n");
      }
    };
    current_pos++;
  };
  return (ans);
}

/**
 * Creates a vertex and adds it to the start of the supplied
 * VertexList
 * @param a an int specifying the vertex
 * @param vertices a VertexList
 * @return a VertexList with a Vertex containing a prepended to vertices
 */
VertexList prepend_vertex(const int a, const VertexList vertices) {
  VertexList new_ptr;
  new_ptr = malloc(sizeof(Vertex_t));
  (new_ptr->data) = a;
  (new_ptr->nextVertex) = vertices;
  return new_ptr;
}

/**
 * Add a new edge to the start of a EdgeList
 * @param a an int indictating the source of the edge
 * @param b an int indictating the destination of the edge
 * @param e the EdgeList
 * @return e with the edge (a,b) at the front
 */
EdgeList prepend_edge(const int a, const int b, const EdgeList e) {
  EdgeList new_ptr;
  new_ptr = malloc(sizeof(EdgeNode_t));
  new_ptr->start = a;
  new_ptr->end = b;
  new_ptr->nextEdge = e;
  return (new_ptr);
}

/**
 * Creates a single element ShortEdge
 * @param a an int indicating the source of the edge
 * @param b an int indicating the destination of the edge
 * @return a single element EdgeList (a,b)
 */
EdgeList create_edge(const int a, const int b) {
  EdgeList new_ptr;
  new_ptr = malloc(sizeof(EdgeNode_t));
  new_ptr->start = a;
  new_ptr->end = b;
  new_ptr->nextEdge = NULL;
  return (new_ptr);
}

/**
 * Returns a number less than 0, equal to 0, or greater than zero
 * if respectively the head of e1 is less than, equal to, or greater than
 * the head of e2 in dictionary order by start then end.
 * @param e1 a nonempty edge list
 * @param e2 a nonempty edge list
 * @return a comparison of the heads of e1 and e2
 */
int compare_edge(EdgeList e1, EdgeList e2) {
  if (e1->start == e2->start) {
    return e1->end - e2->end;
  }
  return e1->start - e2->end;
}

/**
 * Merges two sorted ascending order edge lists into an ascending order sorted
 * edge list
 * @param list1 a sorted edge list
 * @param list2 a sorted edge list
 * @return a sorted edge list
 * @warning destructively modifies list1 and list2
 */
EdgeList merge_edges(EdgeList list1, EdgeList list2) {
  EdgeList new_list;
  EdgeList tail;

  if (NULL == list1) {
    return list2;
  }
  else if (NULL == list2) {
    return list1;
  }

  if(compare_edge(list1,list2) < 0) {
    new_list = list1;
    tail = list1;
    list1 = list1->nextEdge;
    tail->nextEdge = NULL;
  }
  else {
    new_list = list2;
    tail = list2;
    list2 = list2->nextEdge;
    tail->nextEdge = NULL;
  }

  while(NULL != list1 && NULL != list2) {
    if(compare_edge(list1,list2) > 0) {
      tail->nextEdge = list1;
      tail = tail->nextEdge;
      list1 = list1->nextEdge;
      tail->nextEdge = NULL;
    }
    else {
      tail->nextEdge = list2;
      tail = tail->nextEdge;
      list2 = list2->nextEdge;
      tail->nextEdge = NULL;
    }
  }

  if(NULL != list1) {
    tail->nextEdge = list1;
  }
  else {
    tail->nextEdge = list2;
  }
  
  return new_list;
}

/**
 * Sorts an edge list into ascending order in dictionary order by
 * start and then end
 * @param an edge list
 * @return the sorted version of edge_list
 * @warning destructively modifies edge_list
 */
EdgeList merge_sort_edges(EdgeList edge_list) {
  EdgeList list1, list2, tail1, tail2;

  if(NULL == edge_list || NULL == edge_list->nextEdge) {
    return edge_list;
  }

  list1 = edge_list;
  tail1 = edge_list;
  edge_list = edge_list->nextEdge;
  tail1->nextEdge = NULL;
  list2 = edge_list;
  tail2 = edge_list;
  edge_list = edge_list->nextEdge;
  tail2->nextEdge = NULL;

  int which_list = 1;
  while(NULL != edge_list) {
    if (which_list) {
      tail1->nextEdge = edge_list;
      tail1 = tail1->nextEdge;
      tail1->nextEdge = NULL;
    }
    else {
      tail2->nextEdge = edge_list;
      tail2 = tail2->nextEdge;
      tail2->nextEdge = NULL;
    }

    edge_list = edge_list->nextEdge;
    which_list = !which_list;
  }

  merge_sort_edges(list1);
  merge_sort_edges(list2);
  
  return merge_edges(list1, list2);
}

/**
 * Frees the supplied EdgeList
 * @param e a EdgeList
 */
void free_edge_list(const EdgeList e) {
  EdgeList temp, ntemp;
  temp = e;
  ntemp = temp;
  while (temp != NULL) {
    ntemp = temp;
    temp = temp->nextEdge;
    free(ntemp);
  };
}

/**
 * Takes in a list of parent vertices and a list of child vertices,
 * generates all edges between them and adds this EdgeList to passed edge_list
 * mod 2
 * @param parents a list of parent vertices
 * @param kids a list of child vertices
 * @param edge_list list of Edges between vertices
 * @return A shortEdges containing the result of adding them mod two.
 */
EdgeList add_mod_two_lists(const VertexList parents, const VertexList kids,
                           const EdgeList *const edge_list) {
  VertexList this_kid, this_parent, temp_vert;
  EdgeList this_edge, temp_edge, prev;
  EdgeList ans;
  if ((parents == NULL) || (kids == NULL)) {
    return *edge_list;
  } else {
    this_parent = parents;
    this_kid = kids;
    ans = NULL;
    this_edge = *edge_list;
    while (this_parent != NULL && this_edge != NULL &&
           (this_edge->start == this_parent->data &&
            this_edge->end == this_kid->data)) {
      temp_edge = this_edge;
      this_edge = this_edge->nextEdge;
      this_kid = this_kid->nextVertex;
      free(temp_edge);
      if (this_kid == NULL) {
        temp_vert = this_parent;
        this_parent = this_parent->nextVertex;
        free(temp_vert);
        this_kid = kids;
      };
    };
    if (this_edge != NULL &&
        (this_parent == NULL || (this_edge->start < this_parent->data ||
                                 (this_edge->start == this_parent->data &&
                                  this_edge->end < this_kid->data)))) {
      ans = this_edge;
      prev = ans;
      this_edge = this_edge->nextEdge;
    } else if (this_parent != NULL) {
      ans = create_edge(this_parent->data, this_kid->data);
      ans->nextEdge = NULL;
      this_kid = this_kid->nextVertex;
      if (this_kid == NULL) {
        temp_vert = this_parent;
        this_parent = this_parent->nextVertex;
        free(temp_vert);
        this_kid = kids;
      };
    } else {
      ans = NULL;
    }
    
    prev = ans;

    if (NULL == this_edge) {
      while (this_parent != NULL) {
        prev->nextEdge = create_edge(this_parent->data, this_kid->data);
        prev = prev->nextEdge;
        prev->nextEdge = NULL;
        this_kid = this_kid->nextVertex;
        if (this_kid == NULL) {
          temp_vert = this_parent;
          this_parent = this_parent->nextVertex;
          free(temp_vert);
          this_kid = kids;
        }
      }
    }
    
    while (this_edge != NULL && this_parent != NULL) {
      while (this_edge != NULL && ((this_edge->start < this_parent->data ||
                                    (this_edge->start == this_parent->data &&
                                     this_edge->end < this_kid->data)))) {
        prev->nextEdge = this_edge;
        prev = prev->nextEdge;
        this_edge = this_edge->nextEdge;
      };
      while (this_edge != NULL && this_parent != NULL &&
             (this_edge->start == this_parent->data &&
              this_edge->end == this_kid->data)) {
        temp_edge = this_edge;
        this_edge = temp_edge->nextEdge;
        prev->nextEdge = this_edge;
        free(temp_edge);
        this_kid = this_kid->nextVertex;
        if (this_kid == NULL) {
          temp_vert = this_parent;
          this_parent = this_parent->nextVertex;
          free(temp_vert);
          this_kid = kids;
        };
      };
      while (this_parent != NULL &&
             ((this_edge == NULL) || ((this_edge->start > this_parent->data ||
                                       (this_edge->start == this_parent->data &&
                                        this_edge->end > this_kid->data))))) {
        prev->nextEdge = create_edge(this_parent->data, this_kid->data);
        prev = prev->nextEdge;
        prev->nextEdge = NULL;
        this_kid = this_kid->nextVertex;
        if (this_kid == NULL) {
          temp_vert = this_parent;
          this_parent = this_parent->nextVertex;
          free(temp_vert);
          this_kid = kids;
        };
      };
    };
    if ((this_parent == NULL) && (prev != NULL)) {
      prev->nextEdge = this_edge;
    }
  };
  return (ans);
}

/**
 * contracts all edges such that the parents occur after init
 * and the children are before or at final.
 * @param init an int specifying the required start
 * @param final State to t
 * @param edge_list the EdgeList
 */
void special_homology(const int init, const int final, EdgeList *edge_list) {
  int i, j;
  EdgeList temp;
  i = 0;
  j = 0;
  temp = *edge_list;
  while ((*edge_list != NULL) && (temp != NULL)) {
    while ((temp != NULL) && (temp->start == init)) {
      temp = temp->nextEdge;
    };
    while ((temp != NULL) && (temp->end > final)) {
      temp = temp->nextEdge;
    };
    i++;
    j++;
    if (temp != NULL) {
      contract(temp->start, temp->end, edge_list);
      temp = *edge_list;
    };
  };
}

/**
 * Places a new edge into the supplied edgelist in order.
 * @param a an int specifying the parent vertex
 * @param b an int specifying the child vertex
 * @param edges a EdgeList where the new edge will be placed
 * @return a pointer to edges with the new edge added
 */
EdgeList append_ordered(const int a, const int b, const EdgeList edges) {
  EdgeList temp, prev, curr, ans;
  prev = edges;
  if ((edges == NULL) || (edges->start > a) ||
      (edges->start == a && edges->end > b)) {
    ans = malloc(sizeof(EdgeNode_t));
    ans->start = a;
    ans->end = b;
    ans->nextEdge = prev;
  } else {
    ans = edges;
    curr = prev->nextEdge;
    while (curr != NULL &&
           ((curr->start < a) || ((curr->start == a) && (curr->end < b)))) {
      curr = curr->nextEdge;
      prev = prev->nextEdge;
    };
    temp = malloc(sizeof(EdgeNode_t));
    temp->start = a;
    temp->end = b;
    temp->nextEdge = curr;
    prev->nextEdge = temp;
  };
  return (ans);
}

/**
 * contracts the edge specified by the input within global_edge_list
 * @param a the parent vertex of the edge
 * @param b the child vertex of the edge
 * @param edge_list the EdgeList
 */
void contract(const int a, const int b, EdgeList *edge_list) {
  // Initialization
  EdgeList temp;
  EdgeList prev;
  VertexList parents, kids, temp_kids, temp_parents;
  VertexList last_parent, last_kid;
  prev = *edge_list; // Multiple equal initializations
  parents = NULL;
  kids = NULL;
  temp_kids = NULL;
  temp_parents = NULL;
  last_parent = NULL;
  last_kid = NULL;
  // Loops through edge_list edges that end at b or start at a
  while (*edge_list != NULL &&
         ((*edge_list)->end == b || (*edge_list)->start == a)) {
    // If the edge goes from a to b we remove it from the edge_list
    if (((*edge_list)->end == b) && ((*edge_list)->start == a)) {
      temp = *edge_list;
      *edge_list = (*edge_list)->nextEdge;
      free(temp);
    } else {
      // if the edge ends at b but doesnt start at a we add the parent vertex of
      // the edge to a VertexList called last_parents
      if ((*edge_list)->end == b) {
        // Initializes last_parent if it hasnt yet
        if (last_parent == NULL) {
          parents = prepend_vertex((*edge_list)->start, NULL);
          last_parent = parents;
        } else {
          // Adds vertex to last_parent
          temp_parents = prepend_vertex((*edge_list)->start, NULL);
          last_parent->nextVertex = temp_parents;
          last_parent = last_parent->nextVertex;
        };
        // Removes the edge from edge_list
        temp = *edge_list;
        *edge_list = (*edge_list)->nextEdge;
        free(temp);
      } else if ((*edge_list)->start == a) {
        // If the edge starts at a but doesnt end at b we add the kid vertex of
        // the edge to a VertexList called last_kids

        // Initializes last_kid if it hasnt yet
        if (last_kid == NULL) {
          kids = prepend_vertex((*edge_list)->end, kids);
          last_kid = kids;
        } else {
          // Adds vertex to last_kids
          temp_kids = prepend_vertex((*edge_list)->end, NULL);
          last_kid->nextVertex = temp_kids;
          last_kid = last_kid->nextVertex;
        };
        // Removes the edge from edge_list
        temp = *edge_list;
        *edge_list = (*edge_list)->nextEdge;
        free(temp);
      };
    };
  };
  // Initializes prev and temp for next 2 loops
  prev = *edge_list;
  if (*edge_list != NULL) {
    temp = ((*edge_list)->nextEdge);
  } else {
    temp = NULL;
  };
  // Loops through edges that have a start vertex occuring before a
  while (temp != NULL && (temp)->start < a) {
    // If corresponding edge ends at b we add he parent of the edge to the
    // VertexList last_parent
    if ((temp)->end == b) {
      // Initializes last_parent if it hasnt before
      if (last_parent == NULL) {
        parents = prepend_vertex((temp)->start, NULL);
        last_parent = parents;
      } else {
        // adds vertex to last_parent
        temp_parents = prepend_vertex((temp)->start, NULL);
        last_parent->nextVertex = temp_parents;
        last_parent = last_parent->nextVertex;
      };
      // Iterates prev and temp to nextEdge and removes edge from temp edge_list
      (prev->nextEdge) = (temp->nextEdge);
      free(temp);
      temp = prev->nextEdge;
    } else {
      // Iterates prev and temp
      temp = (temp)->nextEdge;
      prev = (prev)->nextEdge;
    };
  };
  // Loops through edges in temp starting at a
  while (temp != NULL && (temp)->start == a) {
    // if the edge does not end at b we add its kid vertex to the vertex_list
    // last_kid
    if (temp->end != b) {
      if (last_kid == NULL) {
        kids = prepend_vertex(temp->end, NULL);
        last_kid = kids;
      } else {
        temp_kids = prepend_vertex(temp->end, NULL);
        last_kid->nextVertex = temp_kids;
        last_kid = last_kid->nextVertex;
      };
    };
    // Iterates prev and temp to the nextEdge and removes edge from temp
    // edge_list
    (prev)->nextEdge = temp->nextEdge;
    free(temp);
    temp = (prev)->nextEdge;
  };
  // Loops through remaining temp edges in order to get edges pointing from
  // vertex
  // occuring after a to b
  while (temp != NULL) {
    // If the edge ends at b we add its parent to VertexList parents
    if ((temp)->end == b) {
      if (last_parent == NULL) {
        parents = prepend_vertex(temp->start, NULL);
        last_parent = parents;
      } else {
        temp_parents = prepend_vertex((temp)->start, NULL);
        last_parent->nextVertex = temp_parents;
        last_parent = last_parent->nextVertex;
      };
      // Iterate and remove
      (prev)->nextEdge = (temp)->nextEdge;
      free(temp);
      temp = (prev)->nextEdge;
    } else {
      // Iterate
      temp = (temp)->nextEdge;
      prev = (prev)->nextEdge;
    }
  };
  // Modifies edge_list based on calculated parents and kids.
  // see add_mod_two_lists
  *edge_list = add_mod_two_lists(parents, kids, edge_list);
}

/**
 * Returns a StateList of states where a rectangle exists from incoming
 * that is not contained in prevs.
 * @param prevs Statelist containing previous states
 * @param incoming the source of rectangles used to generate statelist
 * @param G working grid
 * @return A statelist containing states reached from a rectangle leaving
 * incoming not contained in prevs.
 */
StateRBTree new_rectangles_out_of(const StateRBTree prevs, const State incoming,
                                const Grid_t *const G) {
  StateRBTree ans;
  State temp_state = malloc(sizeof(char) * G->arc_index);
  State temp;
  int LL;
  int w, h, i;
  ans = EMPTY_TREE;
  i = 0;
  while (i < G->arc_index) {
    temp_state[i] = incoming[i];
    i++;
  }
  LL = 0;
  while (LL < G->arc_index) {
    w = 1;
    h = min(mod(G->Os[LL] - incoming[LL], G->arc_index), mod(G->Xs[LL] - incoming[LL], G->arc_index));
    while (w < G->arc_index && h > 0) {
      if (mod(incoming[mod(LL + w, G->arc_index)] - incoming[LL], G->arc_index) <= h) {
        temp_state[LL] = incoming[mod(LL + w, G->arc_index)];
        temp_state[mod(LL + w, G->arc_index)] = incoming[LL];
        if (!s_is_member(&prevs, temp_state, G)) {
          if (s_is_member(&ans, temp_state, G)) {
            s_delete_data(&ans, temp_state, G);
          } else {
            temp = swap_cols(LL, mod(LL + w, G->arc_index), incoming, G);
            s_insert_data(&ans, temp, G);
          }
        };
        temp_state[LL] = incoming[LL];
        temp_state[mod(LL + w, G->arc_index)] = incoming[mod(LL + w, G->arc_index)];
        h = mod(incoming[mod(LL + w, G->arc_index)] - incoming[LL], G->arc_index);
      };
      h = min(h, min(mod(G->Os[mod(LL + w, G->arc_index)] - incoming[LL], G->arc_index),
                     mod(G->Xs[mod(LL + w, G->arc_index)] - incoming[LL], G->arc_index)));
      w++;
    };
    LL++;
  }
  
  free(temp_state);
  return ans;
}

/**
 * returns a StateRBTree containing those with a rectangle
 * pointing to the state incoming that do not overlap with prevs
 * @param incoming State that is the destination for generated rectangles
 * @param prevs StateRBTree of excluded states
 * @param G working grid
 * @return StateRBTree containing states with a rectangle to incoming.
 */
StateRBTree new_rectangles_into(const StateRBTree prevs, const State incoming,
                              const Grid_t *const G) {
  StateRBTree ans;
  State temp_state = malloc(sizeof(char) * G->arc_index);
  State temp;
  int LL;
  int w, h;
  int i;
  ans = EMPTY_TREE;
  i = 0;
  while (i < G->arc_index) {
    temp_state[i] = incoming[i];
    i++;
  }
  LL = 0;
  while (LL < G->arc_index) {
    w = 1;
    h = min(mod_up(incoming[LL] - G->Os[LL], G->arc_index),
            mod_up(incoming[LL] - G->Xs[LL], G->arc_index));
    while (w < G->arc_index && h > 0) {
      if (mod_up(incoming[LL] - incoming[mod(LL + w, G->arc_index)], G->arc_index) < h) {
        temp_state[LL] = incoming[mod(LL + w, G->arc_index)];
        temp_state[mod(LL + w, G->arc_index)] = incoming[LL];
        if (!s_is_member(&prevs, temp_state, G)) {
          if (s_is_member(&ans, temp_state, G)) {
            s_delete_data(&ans, temp_state, G);
          } else {
            temp = swap_cols(LL, mod(LL + w, G->arc_index), incoming, G);
            s_insert_data(&ans, temp, G);
          }
        };
        temp_state[LL] = incoming[LL];
        temp_state[mod(LL + w, G->arc_index)] = incoming[mod(LL + w, G->arc_index)];
        h = mod_up(incoming[LL] - incoming[mod(LL + w, G->arc_index)], G->arc_index);
      };
      h = min(h, min(mod_up(incoming[LL] - G->Os[mod(LL + w, G->arc_index)], G->arc_index),
                     mod_up(incoming[LL] - G->Xs[mod(LL + w, G->arc_index)], G->arc_index)));
      w++;
    };
    LL++;
  };
  return ans;
}

/**
 * Calculates a StateList containing all states reachable
 * by a rectangles of a fixed weight
 * @param wt an int specifying rectangle width
 * @param incoming origin state for the rectangles
 * @param G working grid
 * @return a StateList with states that are reached by a rectangle of width
 * wt from incoming.
 */
StateList fixed_wt_rectangles_out_of(const int wt, const State incoming,
                                     const Grid_t *const G) {
  StateList temp, ans;
  int LL;
  int w, h;
  int this_weight, i;
  ans = NULL;
  LL = 0;
  while (LL < G->arc_index) {
    w = 1;
    h = mod(G->Os[LL] - incoming[LL], G->arc_index);
    while (w < G->arc_index && h > 0) {
      if (mod(incoming[mod(LL + w, G->arc_index)] - incoming[LL], G->arc_index) <= h) {
        this_weight = 0;
        i = 0;
        while (i < w && this_weight <= wt + 1) {
          if (mod(G->Xs[mod(LL + i, G->arc_index)] - incoming[LL], G->arc_index) <
              mod(incoming[mod(LL + w, G->arc_index)] - incoming[LL], G->arc_index)) {
            this_weight++;
          };
          i++;
        };
        if (this_weight == wt) {
          temp = swap_cols_list(LL, mod(LL + w, G->arc_index), incoming, G);
          if (get_number(temp->data, ans, G) != 0) {
            ans = remove_state(temp->data, ans, G);
            free(temp->data);
            free(temp);
            temp = ans;
          } else {
            temp->nextState = ans;
            ans = temp;
          }
        };
        h = mod(incoming[mod(LL + w, G->arc_index)] - incoming[LL], G->arc_index);
      };
      h = min(h, mod(G->Os[mod(LL + w, G->arc_index)] - incoming[LL], G->arc_index));
      w++;
    };
    LL++;
  };
  return ans;
}

/**
 * Finds all lift states that are leaving the state incoming on G that are not in prevs. Accounts
 * for if the grid has been mirrored to calculate rectangles in.
 * @param prevs a lift state list containing previously encountered states
 * @param incoming the lift state that rectangles will be leaving
 * @param G a grid
 * @param is_mirrored pass 1 if the grid has been mirrored, 0 otherwise
 * @return a lift state list containing lift states that can be reached from incoming that are not in prevs
 */
static LiftStateRBTree new_lift_rectangles_out_internal(const LiftStateRBTree prevs, const LiftState incoming, const LiftGrid_t *const G, int is_mirrored) {
  LiftStateRBTree ans = EMPTY_LIFT_TREE;
  double *g_Xs, *g_Os;

  g_Xs = malloc(sizeof(double)*G->arc_index);
  g_Os = malloc(sizeof(double)*G->arc_index);

  for(int i=0; i < G->arc_index; ++i) {
    g_Xs[i] = ((double)G->Xs[i]) - .5;
    g_Os[i] = ((double)G->Os[i]) - .5;
  }

  for(int start_x=0; start_x < G->arc_index*G->sheets; ++start_x) {
    int jumped_down = 0;
    int jumped_up = 0;
    int start_y = mod(incoming[start_x/G->arc_index][start_x%G->arc_index]-1, G->arc_index);
    int start_sheet = start_x/G->arc_index;
    int step = 0;
    int check_index = mod((start_x + step) % G->arc_index, G->arc_index);
    int jump = start_sheet*G->arc_index;
    int height = mod((start_y - 1) % G->arc_index, G->arc_index);

    while (height != start_y) {
      check_index = mod((start_x + step)%G->arc_index, G->arc_index);
      int check_index_gen = mod((mod((start_x + step + 1) % G->arc_index, G->arc_index) + jump) % (G->arc_index * G->sheets), G->arc_index * G->sheets);
      int clear = 1;
      if (height > start_y) {
        if (g_Xs[check_index] < height && g_Xs[check_index] > start_y && clear) {
          clear = 0;
        }
        if (g_Os[check_index] < height && g_Os[check_index] > start_y && clear) {
          clear = 0;
        }
        if (g_Xs[check_index] > height && g_Os[check_index] < start_y && clear) {
          jump = jump + G->arc_index;
          jumped_up = 1;
          check_index_gen = mod((mod((start_x + step + 1) % G->arc_index, G->arc_index) + jump) % (G->arc_index * G->sheets), G->arc_index * G->sheets);
        }
        if (g_Os[check_index] > height && g_Xs[check_index] < start_y && clear) {
          jump = jump - G->arc_index;
          jumped_down = 1;
          check_index_gen = mod((mod((start_x + step + 1) % G->arc_index, G->arc_index) + jump) % (G->arc_index * G->sheets), G->arc_index * G->sheets);
        }
        if (mod(incoming[check_index_gen/G->arc_index][check_index_gen%G->arc_index]-1, G->arc_index) < height && mod(incoming[check_index_gen/G->arc_index][check_index_gen%G->arc_index]-1, G->arc_index) > start_y && clear) {
          if (jumped_down) {
            jumped_down = 0;
            jump = jump + G->arc_index;
          }
          if (jumped_up) {
            jumped_up = 0;
            jump = jump - G->arc_index;
          }
          clear = 0;
        }
        if (clear) {
          check_index_gen = mod((mod((start_x + step + 1) % G->arc_index, G->arc_index) + jump) % (G->arc_index * G->sheets), G->arc_index*G->sheets);
          if (mod(incoming[check_index_gen/G->arc_index][check_index_gen%G->arc_index]-1, G->arc_index) == height) {
            LiftState new_state = NULL;
            init_lift_state(&new_state, G);
            copy_lift_state(&new_state, &incoming, G);
            new_state[start_x/G->arc_index][start_x%G->arc_index] = incoming[check_index_gen/G->arc_index][check_index_gen%G->arc_index];
            new_state[check_index_gen/G->arc_index][check_index_gen%G->arc_index] = incoming[start_x/G->arc_index][start_x%G->arc_index];
            if (is_mirrored) {
              mirror_lift_state(&new_state, G);
            }

            if(!is_member(&prevs, new_state, G)) {
              if(!is_member(&ans, new_state, G)) {
                insert_data(&ans, new_state,G);
              }
              else {
                LiftStateRBTree temp = find_node(&ans, new_state, G);
                delete_node(&ans, temp);
                free_lift_state(&(temp->data), G);
                free(temp);
                free_lift_state(&new_state, G);
              }
            }
            
            height = mod((height - 1) % G->arc_index, G->arc_index);
          }
          step = step + 1;
          jumped_down = 0;
          jumped_up = 0;          
        }
        else {
          height = mod((height - 1) % G->arc_index,G->arc_index);
        }
      }
      else {
        if ((g_Xs[check_index] < height || g_Xs[check_index] > start_y) && clear) {
          clear = 0;
        }
        if ((g_Os[check_index] < height || g_Os[check_index] > start_y) && clear) {
          clear = 0;
        }
        if ((mod(incoming[check_index_gen/G->arc_index][check_index_gen%G->arc_index]-1, G->arc_index) < height || mod(incoming[check_index_gen/G->arc_index][check_index_gen%G->arc_index]-1, G->arc_index) >= start_y) && clear) {
          clear = 0;
        }
        if (clear) {
          if (mod(incoming[check_index_gen/G->arc_index][check_index_gen%G->arc_index]-1, G->arc_index) == height) {
            LiftState new_state = NULL;
            init_lift_state(&new_state, G);
            copy_lift_state(&new_state, &incoming, G);
            new_state[start_x/G->arc_index][start_x%G->arc_index] = incoming[check_index_gen/G->arc_index][check_index_gen%G->arc_index];
            new_state[check_index_gen/G->arc_index][check_index_gen%G->arc_index] = incoming[start_x/G->arc_index][start_x%G->arc_index];
            if (is_mirrored) {
              mirror_lift_state(&new_state, G);
            }

            if(!is_member(&prevs, new_state, G)) {
              if(!is_member(&ans, new_state, G)) {
                insert_data(&ans, new_state,G);
              }
              else {
                LiftStateRBTree temp = find_node(&ans, new_state, G);
                delete_node(&ans, temp);
                free_lift_state(&(temp->data), G);
                free(temp);
                free_lift_state(&new_state, G);
              }
            }

            height = mod((height -1) % G->arc_index,G->arc_index);
          }
          step = step + 1;
        }
        else {
          height = mod((height - 1) % G->arc_index, G->arc_index);
        }
      }
    }
  }

  free(g_Xs);
  free(g_Os);
  return ans;
}

/**
 * returns a LiftStateList containing those with a rectangle
 * pointing from the Liftstate incoming that do not overlap with prevs
 * @param incoming LiftState that is the destination for generated rectangles
 * @param prevs LiftStateList of excluded states
 * @param G working grid
 * @return LiftStateList containing states with a rectangle to incoming.
 */
LiftStateRBTree new_lift_rectangles_out_of(const LiftStateRBTree prevs, const LiftState incoming, const LiftGrid_t *const G) {
  return new_lift_rectangles_out_internal(prevs, incoming, G, 0);
}

/**
 * returns a LiftStateList containing those with a rectangle
 * pointing to the LiftState incoming that do not overlap with prevs
 * @param incoming LiftState that is the destination for generated rectangles
 * @param prevs LiftStateList of excluded states
 * @param G working grid
 * @return LiftStateList containing states with a rectangle to incoming.
 */
LiftStateRBTree new_lift_rectangles_into(const LiftStateRBTree prevs, const LiftState incoming, const LiftGrid_t *const G) {
  LiftGrid_t* G_mirror = mirror_lift_grid(G);
  LiftState incoming_mirror;
  init_lift_state(&incoming_mirror, G);
  copy_lift_state(&incoming_mirror, &incoming, G);
  mirror_lift_state(&incoming_mirror, G);
  LiftStateRBTree ans = new_lift_rectangles_out_internal(prevs, incoming_mirror, G_mirror, 1);

  free(G_mirror);
  free_lift_state(&incoming_mirror, G);

  return ans;
}

/**
 * Prints the permutation of state on the grid specified by
 * Xs and Os, as well as 2A=M=SL+1.
 * @param state a State
 * @param G working grid
 */
void print_state(const State state, const Grid_t *const G) {
  int i, j;
  j = G->arc_index;
  i = 0;
  (*print_ptr)("*---");
  while (i < G->arc_index - 1) {
    (*print_ptr)("----");
    i++;
  }
  (*print_ptr)("*\n");
  while (j > 0) {
    i = 0;
    while (i < G->arc_index) {
      if (G->Xs[i] == j) {
        (*print_ptr)("| X ");
      } else {
        if (G->Os[i] == j) {
          (*print_ptr)("| O ");
        } else {
          (*print_ptr)("|   ");
        };
      };
      i++;
    };
    (*print_ptr)("|\n");
    i = 0;
    while (i < G->arc_index) {
      if (state[i] == j) {
        (*print_ptr)("@---");
      } else {
        if (i == 0 && j > 1) {
          (*print_ptr)("|---");
        } else {
          if (j > 1) {
            (*print_ptr)("+---");
          } else {
            if (i == 0) {
              (*print_ptr)("*---");
            } else {
              if (i == 0) {
                (*print_ptr)("----");
              } else {
                (*print_ptr)("----");
              };
            };
          };
        };
      };
      i++;
    };
    if (j > 1) {
      (*print_ptr)("|\n");
    } else {
      (*print_ptr)("*\n");
    };
    j--;
  };
  (*print_ptr)("\n");
}

/**
 * Prints the permutation of state using one line notation
 * "{_,_,...}"
 * @param state a State
 * @param G working grid
 */
void print_state_short(const State state, const Grid_t *const G) {
  int i;
  i = 0;
  (*print_ptr)("{");
  while (i < G->arc_index - 1) {
    (*print_ptr)("%d,", state[i]);
    i++;
  };
  (*print_ptr)("%d}\n", state[G->arc_index - 1]);
}

/**
 * Calls print_state on each sheet of the lift state
 * @param state a lift state
 * @param G a lift grid
 * @see print_state
 */
void print_lift_state(const LiftState state, const LiftGrid_t * const G) {
  Grid_t H;
  H.arc_index = G->arc_index;
  H.Xs = G->Xs;
  H.Os = G->Os;

  for(int i=0; i < G->sheets; ++i) {
    print_state(state[i], &H);
  }
}

/**
 * Prints the permutations of a lift state
 * @param state a lift state
 * @param G a grid
 * @see print_state_short
 */
void print_lift_state_short(const LiftState state, const LiftGrid_t * const G) {
  Grid_t G_p;
  G_p.arc_index = G->arc_index;
  G_p.Xs = G->Xs;
  G_p.Os = G->Os;
  for(int i = 0; i < G->sheets; ++i) {
    print_state_short(state[i], &G_p);
  }
}

/**
 * Prints states in the form "{<state>,...}" up
 * to the first 500,000 states.
 * @param states a StateList
 * @param G working grid
 * @see print_state_short
 */
void print_states(const StateList states, const Grid_t *const G) {
  StateList temp;
  int c;
  temp = states;
  (*print_ptr)("{");
  c = 0;
  while ((temp != NULL) && c < 500000) {
    print_state_short(temp->data, G);
    temp = temp->nextState;
    if (temp != NULL) {
      (*print_ptr)(",");
    };
    c++;
  };
  if (c == 500000) {
    (*print_ptr)("...");
  };
  (*print_ptr)("}");
}

/**
 * Prints states in the form "{<state>,...}" up to the first 500,000 states.
 * @param states a LiftStateList
 * @param G a lift grid
 * @see print_states
 */
void print_lift_states(const LiftStateList states, const LiftGrid_t *const G) {
  LiftStateList temp;
  int c;
  temp = states;
  (*print_ptr)("{");
  c = 0;
  while ((temp != NULL) && c < 500000) {
    print_lift_state_short(temp->data, G);
    temp = temp->nextState;
    if (temp != NULL) {
      (*print_ptr)(",");
    };
    c++;
  };
  if (c == 500000) {
    (*print_ptr)("...");
  };
  (*print_ptr)("}");
}

/**
 * Prints all states in the supplied tree using print_state_short
 * @param states a StateRBTree
 * @param G a grid
 * @see print_state_short
 */
void print_states_tree(const StateRBTree states, const Grid_t * const G) {
  if(EMPTY_TREE == states) {
    return;
  }
  print_state_short(states->data, G);
  print_states_tree(states->left, G);
  print_states_tree(states->right, G);  
}

/**
 * Prints all lift states in the supplied tree using print_lift_state_short
 * @param states a LiftStateRBTree
 * @param G a lift grid
 * @see print_lift_state_short
 */
void print_states_lift_tree(const LiftStateRBTree states, const LiftGrid_t * const G) {
  if(EMPTY_LIFT_TREE == states) {
    return;
  }
  print_lift_state(states->data, G);
  print_states_lift_tree(states->left, G);
  print_states_lift_tree(states->right, G);
}

/**
 * As print_states_tree but also prints the tags attatched to each state
 * @param states a StateRBTree
 * @param G a grid
 * @see print_states_tree
 */
void print_states_tags(const StateRBTree states, const Grid_t * const G) {
  if(EMPTY_TREE == states) {
    return;
  }
  (*print_ptr)("%d, ", states->tag);
  print_state_short(states->data, G);
  print_states_tags(states->left, G);
  print_states_tags(states->right, G);
}

/**
 * As print_lift_states_tree but also prints the tags attatched to each lift state
 * @param states a LiftStateRBTree
 * @param G a lift grid
 * @see print_lift_states_tree
 */
void print_states_lift_tags(const LiftStateRBTree states, const LiftGrid_t * const G) {
  if(EMPTY_LIFT_TREE == states) {
    return;
  }
  (*print_ptr)("%d, ", states->tag);
  print_lift_state(states->data, G);
  print_states_lift_tree(states->left, G);
  print_states_lift_tree(states->right, G);
}

/**
 * Prints each edge in the passed EdgeList
 * @param edge_list an EdgeList
 */
void print_edges(const EdgeList edge_list) {
  EdgeList temp;
  temp = edge_list;
  while (temp != NULL) {
    (*print_ptr)("[%d -> %d]\n", temp->start, temp->end);
    temp = (temp->nextEdge);
  };
}

/**
 * Print the first 80 edges edge_list on the same line
 * @param edge_list an EdgeList
 */
void print_math_edges(const EdgeList edge_list) {
  EdgeList temp;
  int t;
  temp = edge_list;
  (*print_ptr)("{");
  t = 0;
  while (temp != NULL) {
    (*print_ptr)("[%d -> %d]", temp->start, temp->end);
    t++;
    if (t == 80) {
      temp = NULL;
      (*print_ptr)("...");
    } else {
      temp = (temp->nextEdge);
      if (temp != NULL)
        (*print_ptr)(",");
    }
  };
  (*print_ptr)("}\n");
}

/**
 * Prints the edges in a single line
 * @param edges an EdgeList
 */
void print_math_edges_a(const EdgeList edges) {
  EdgeList temp;
  temp = edges;
  (*print_ptr)("{");
  while (temp != NULL) {
    (*print_ptr)("[%d->%d]", temp->start, temp->end);
    temp = (temp->nextEdge);
    if (temp != NULL)
      (*print_ptr)(",");
  };
  (*print_ptr)("}");
}

/**
 * Prints the vertices in VertexList on a single line
 * @param v_list a VertexList
 */
void print_vertices(const VertexList v_list) {
  VertexList temp;
  temp = v_list;
  (*print_ptr)("{");
  while (temp != NULL) {
    (*print_ptr)("%d", (temp)->data);
    temp = (temp)->nextVertex;
    if (temp != NULL)
      (*print_ptr)(",");
  };
  (*print_ptr)("}");
}

/**
 * Prints the permutations of the Grid, ie the X O code inputed
 * by the user
 * @param G working grid
 */
void print_grid_perm(const Grid_t *const G) {
  int i = 0;
  (*print_ptr)("X = [");
  while (i < G->arc_index) {
    (*print_ptr)(" %d", G->Xs[i]);
    if (i != G->arc_index - 1)
      (*print_ptr)(",");
    i++;
  }
  (*print_ptr)(" ]\nO = [");
  i = 0;
  while (i < G->arc_index) {
    (*print_ptr)(" %d", G->Os[i]);
    if (i != G->arc_index - 1)
      (*print_ptr)(",");
    i++;
  }
  (*print_ptr)(" ]\n");
}

/**
 * Prints the grid without the grid state
 * @param G working grid
 */
void print_grid(const Grid_t *const G) {
  int i, j;
  j = G->arc_index;
  i = 0;
  (*print_ptr)("*---");
  while (i < G->arc_index - 1) {
    (*print_ptr)("----");
    i++;
  }
  (*print_ptr)("*\n");
  while (j > 0) {
    i = 0;
    while (i < G->arc_index) {
      if (G->Xs[i] == j) {
        (*print_ptr)("| X ");
      } else {
        if (G->Os[i] == j) {
          (*print_ptr)("| O ");
        } else {
          (*print_ptr)("|   ");
        };
      };
      i++;
    };
    (*print_ptr)("|\n");
    i = 0;
    while (i < G->arc_index) {
      if (i == 0 && j > 1) {
        (*print_ptr)("|---");
      } else {
        if (j > 1) {
          (*print_ptr)("+---");
        } else {
           if (i == 0) {
            (*print_ptr)("*---");
          } else {
              (*print_ptr)("----");
            };
          };
        };
      i++;
    };
    if (j > 1) {
      (*print_ptr)("|\n");
    } else {
      (*print_ptr)("*\n");
    };
    j--;
  };
  (*print_ptr)("\n");
  print_grid_perm(G);
  (*print_ptr)("\n");
}

/**
 * Prints Thurston Bennequin and rotation number of grid
 * @param G working grid
 */
void print_tb_r(const Grid_t *const G) {
  int writhe = 0;
  int up_down_cusps[2] = {0, 0};
  int tb;
  int r;
  writhe = get_writhe(G);
  cusps(up_down_cusps, G);
  tb = writhe - .5 * (up_down_cusps[0] + up_down_cusps[1]);
  r = .5 * (up_down_cusps[1] - up_down_cusps[0]);
  (*print_ptr)("tb = %d\n", tb);
  (*print_ptr)("r = %d\n", r);
} 

/**
 * Prints the Alexander and Maslov grading lines, which are calculated
 * from Thurston Bennequin and rotation number
 * @param G working grid
 * @param plus 1 if x^+ print, 0 if x^- print
 */
void print_2AM(const Grid_t *const G, int plus) {
  int writhe = 0;
  int up_down_cusps[2] = {0, 0};
  int tb;
  int r;
  writhe = get_writhe(G);
  cusps(up_down_cusps, G);
  tb = writhe - .5 * (up_down_cusps[0] + up_down_cusps[1]);
  r = .5 * (up_down_cusps[1] - up_down_cusps[0]);
  if(plus == 1) (*print_ptr)("2A(x^+) = M(x^+) = sl(x^+)+1 = %d\n\n"
                    , tb - r + 1);    
  if(plus == 0) (*print_ptr)("2A(x^-) = M(x^-) = %d\n\n", tb + r + 1);
} 
