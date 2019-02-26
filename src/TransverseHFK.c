/* Program for calculating the invariant for transverse knots, to
   accompany "Transverse knots distinguished by Knot Floer Homology"
   by L. Ng, P. S. Ozsvath, and D. P. Thurston. The invariant was defined
   in "Legendrian knots, transverse knots, and combinatorial Floer homology"
   by P. S. Ozsvath, Z. Szabo, and D. P. Thurston*/

/* In progress cleanup of the above listed file by
 * Lucas Meyers <lmeye22@lsu.edu>
 */

#include <argp.h>
#include <limits.h>
#include <math.h>
#include <signal.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Hardcoded max arcindex until dynamic allocation can be tested for
 * speed */
#define MAX_INDEX 30

const char *argp_program_version = "transverseHFK revision 0.0.1";
const char *argp_program_bug_address = "<lmeye22@lsu.edu>";
static const char doc[] =
    "A program to calculate the Legendrian/Transverse knot invariants\
 via the algorithm described in \"Transverse knots distinguished by\
 Knot Floer Homology\" by L. Ng, P. S. Ozsvath, and D. P. Thurston.";

static const char args_doc[] = "-i [ArcIndex] -X [Xs] -O [Os]";

static struct argp_option options[] = {
    {"verbose", 'v', 0, 0, "Produce verbose output", 0},
    {"quiet", 'q', 0, 0, "Don't produce extraneous output", 0},
    {"index", 'i', "ArcIndex", 0, "ArcIndex of the grid", 0},
    {"Xs", 'X', "[...]", 0, "List of Xs", 0},
    {"Os", 'O', "[...]", 0, "List of Os", 0},
    {"timeout", 't', "SECONDS", 0, "Maximum time to run in seconds", 0},
    {0}};

// Temporary location
static error_t parse_opt(int, char *, struct argp_state *);
static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};

bool verbose = false;
int arc_index = -1;
int max_time = -1;
char Xs[MAX_INDEX] = {};
char Os[MAX_INDEX] = {};

struct Grid {
  char Xs[MAX_INDEX];
  char Os[MAX_INDEX];
  int arc_index;
};

typedef struct Grid Grid_t;

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
  char data[MAX_INDEX];
  StateList nextState;
};

typedef char State[MAX_INDEX];

EdgeList add_mod_two_lists(VertexList, VertexList, EdgeList);
VertexList prepend_vertex(int, VertexList);
EdgeList prepend_edge(int, int, EdgeList);
StateList fixed_wt_rectangles_out_of(int, State, Grid_t*);
StateList swap_cols(int, int, State, Grid_t*);
int get_number(State, StateList);
void free_state_list(StateList);
int null_homologous_D0Q(State, Grid_t*, EdgeList);
int null_homologous_D1Q(State, Grid_t*, EdgeList);
void homology(EdgeList);
void special_homology(int, int, EdgeList);
void contract(int, int, EdgeList);
VertexList remove_vertex(int, VertexList, Grid_t*, EdgeList);
EdgeList create_edge(int, int);
void free_edge_list(EdgeList);
StateList remove_state(State, StateList);
int mod(int, Grid_t*);
int mod_up(int, Grid_t*);
StateList new_rectangles_out_of(StateList, State, Grid_t*);
StateList new_rectangles_into(StateList, State, Grid_t*);
EdgeList append_ordered(int, int, EdgeList);

int NESW_pO(char *, Grid_t*);
int NESW_Op(char *, Grid_t*);
int NESW_pp(char *, Grid_t*);

void print_state(State, Grid_t*);
void print_state_short(State, Grid_t*);
void print_edges(EdgeList);
void print_states(StateList, Grid_t*);
void print_math_edges(EdgeList);
void print_math_edges_a(EdgeList);
void print_vertices(VertexList);

void timeout(int);
int perm_len(const char *);
int is_grid(const Grid_t *);
int build_permutation(char *, char *);
int eq_state(State a, State b) { return (!strncmp(a, b, arc_index)); }

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  switch (key) {
  case 'v':
    verbose = true;
    break;
  case 'q':
    verbose = false;
    break;
  case 't':
    max_time = atoi(arg);
    if (max_time <= 0) {
      argp_failure(state, 0, 0, "Invalid timeout");
      exit(1);
    }
    break;
  case 'i':
    arc_index = atoi(arg);
    if (arc_index > MAX_INDEX || arc_index < 2) {
      argp_failure(state, 0, 0, "ArcIndex value out of range");
      exit(1);
    }
    break;
  case 'X':
    if (-1 == build_permutation(Xs, arg)) {
      argp_failure(state, 0, 0, "Malformatted Xs");
      exit(1);
    }
    break;
  case 'O':
    if (-1 == build_permutation(Os, arg)) {
      argp_failure(state, 0, 0, "Malformatted Os");
      exit(1);
    }
    break;
  default:
    break;
  }
  return 0;
}

/**
 * Signal handler for SIGALRM that exits when the signal is received.
 * @param sig
 */
void timeout(int sig) {
  if (SIGALRM == sig) {
    printf("Timeout reached. Terminating\n");
    exit(0);
  }
}

/**
 * Takes in a string and converts it into a permutation.
 * @param String of the form [_,_,...,_] where _ are integers
 * between 1 and MAX_INDEX
 * @param Destination for the permutation
 * @return 0 on success, -1 on failure
 */
int build_permutation(char *perm, char *str) {
  if (str[0] != '[') {
    return -1;
  }

  char *s = &str[1];
  long n = -1;
  int i = 0;

  while (i < MAX_INDEX) {
    errno = 0;
    n = strtol(s, &s, 10);

    if ((errno == ERANGE && (n == LONG_MAX || n == LONG_MIN)) ||
        (errno != 0 && n == 0) || (n < 1 || n > MAX_INDEX)) {
      return -1;
    }

    perm[i] = n;
    ++i;

    if (s[0] == ']') {
      if (s[1] != '\0') {
        return -1;
      }
      break;
    } else if (s[0] == ',') {
      ++s;
    }
  }

  return 0;
}

/**
 * Returns the number of characters in the character array before a zero is
 * found
 * @param p a character array
 * @return the number of characters before a zero is found. Returns -1 if not
 * found
 */
int perm_len(const char *p) {
  for (int i = 0;; ++i) {
    if (0 == p[i]) {
      return i;
    }
  }

  // Will segfault before reaching this point if no 0 is found
  return -1;
}

/**
 * Determines whether the grid specified by the parameters is valid
 * @param G working grid
 * @return 1 if the grid is valid, 0 otherwise
 */
int is_grid(const Grid_t* G) {
  if (perm_len(G->Xs) != G->arc_index || perm_len(G->Os) != G->arc_index || 1 >= G->arc_index || MAX_INDEX <= G->arc_index) {
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

int main(int argc, char **argv) {
  argp_parse(&argp, argc, argv, 0, 0, 0);

  char UR[arc_index];
  int i;

  Grid_t G;
  EdgeList edge_list = NULL;

  G.arc_index = arc_index;
  for(i=0; i< MAX_INDEX; ++i) {
    if(i < G.arc_index) {
      G.Xs[i] = Xs[i];
      G.Os[i] = Os[i];
    }
    else {
      G.Xs[i] = 0;
      G.Os[i] = 0;
    }
  }
  
  if (!is_grid(&G)) {
    printf("Invalid grid\n");
    exit(1);
  }

  if (max_time > 0) {
    if (signal(SIGALRM, timeout) == SIG_ERR) {
      perror("An error occured while setting the timer");
      exit(1);
    }
    alarm(max_time);
  }

  if (verbose) {
    printf("\n \nCalculating graph for LL invariant\n");
    print_state(G.Xs, &G);
  }
  if (null_homologous_D0Q(G.Xs, &G, edge_list)) {
    printf("LL is null-homologous\n");
  } else {
    printf("LL is NOT null-homologous\n");
  }

  if (verbose) {
    printf("\n \nCalculating graph for UR invariant\n");
  }
  if (G.Xs[G.arc_index - 1] == G.arc_index) {
    UR[0] = 1;
  } else {
    UR[0] = (char)G.Xs[arc_index - 1] + 1;
  };
  i = 1;
  while (i <= G.arc_index - 1) {
    if (G.Xs[i - 1] == arc_index) {
      UR[i] = 1;
    } else {
      UR[i] = G.Xs[i - 1] + 1;
    }
    i++;
  }

  if (verbose) {
    print_state(UR, &G);
  }

  if (null_homologous_D0Q(UR, &G, edge_list)) {
    printf("UR is null-homologous\n");
  } else {
    printf("UR is NOT null-homologous\n");
  };

  if (verbose) {
    printf("\n \nCalculating graph for D1[LL] invariant\n");
    print_state(G.Xs, &G);
  }

  if (null_homologous_D1Q(G.Xs, &G, edge_list)) {
    printf("D1[LL] is null-homologous\n");
  } else {
    printf("D1[LL] is NOT null-homologous\n");
  }

  if (verbose) {
    printf("\n \nCalculating graph for D1[UR] invariant\n");
  }

  if (Xs[arc_index - 1] == G.arc_index) {
    UR[0] = 1;
  } else {
    UR[0] = (char)G.Xs[G.arc_index - 1] + 1;
  };
  i = 1;
  while (i <= G.arc_index - 1) {
    if (Xs[i - 1] == G.arc_index) {
      UR[i] = 1;
    } else {
      UR[i] = G.Xs[i - 1] + 1;
    }
    i++;
  }

  if (verbose) {
    print_state(UR, &G);
  }

  if (null_homologous_D1Q(UR, &G, edge_list)) {
    printf("D1[UR] is null-homologous\n");
  } else {
    printf("D1[UR] is NOT null-homologous\n");
  };

  return 0;
}

/**
 * Shifts the input towards the interval [0,arc_index) by
 * a multiple of arc_index
 * @param a An integer
 * @param G working grid
 * @return a shifted towards the interval [0,arc_index)
 */
int mod(int a, Grid_t* G) {
  if (a >= G->arc_index) {
    return (a - G->arc_index);
  } else if (a < 0) {
    return (a + G->arc_index);
  } else {
    return (a);
  };
}

/**
 * Shifts the input towards the interval (0,arc_index] by
 * a multiple of arc_index
 * @param a An integer
 * @param G working grid
 * @return a shifted towards the interval (0,arc_index]
 */
int mod_up(int a, Grid_t* G) {
  if (a > G->arc_index) {
    return (a - G->arc_index);
  } else if (a <= 0) {
    return (a + G->arc_index);
  } else {
    return (a);
  };
}

int min(int a, int b) {
  if (a < b) {
    return (a);
  } else {
    return (b);
  }
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
StateList new_rectangles_out_of(StateList prevs, State incoming, Grid_t* G) {
  StateList temp, ans;
  State temp_state;
  int LL;
  int w, h, m, n, i;
  ans = NULL;
  i = 0;
  while (i < G->arc_index) {
    temp_state[i] = incoming[i];
    i++;
  }
  LL = 0;
  while (LL < G->arc_index) {
    w = 1;
    h = min(mod(G->Os[LL] - incoming[LL], G), mod(G->Xs[LL] - incoming[LL], G));
    while (w < G->arc_index && h > 0) {
      if (mod(incoming[mod(LL + w, G)] - incoming[LL], G) <= h) {
        temp_state[LL] = incoming[mod(LL + w, G)];
        temp_state[mod(LL + w, G)] = incoming[LL];
        m = get_number(temp_state, prevs);
        if (m == 0) {
          n = get_number(temp_state, ans);
          if (n != 0) {
            ans = remove_state(temp_state, ans);
          } else {
            temp = swap_cols(LL, mod(LL + w, G), incoming, G);
            temp->nextState = ans;
            ans = temp;
          }
        };
        temp_state[LL] = incoming[LL];
        temp_state[mod(LL + w, G)] = incoming[mod(LL + w, G)];
        h = mod(incoming[mod(LL + w, G)] - incoming[LL], G);
      };
      h = min(h, min(mod(G->Os[mod(LL + w, G)] - incoming[LL], G),
                     mod(G->Xs[mod(LL + w, G)] - incoming[LL], G)));
      w++;
    };
    LL++;
  };
  return ans;
}

/**
 * returns a StateList containing those with a rectangle
 * pointing to the state incoming that do not overlap with prevs
 * @param incoming State that is the destination for generated rectangles
 * @param prevs StateList of excluded states
 * @param G working grid
 * @return StateList containing states with a rectangle to incoming.
 */
StateList new_rectangles_into(StateList prevs, State incoming, Grid_t* G) {
  StateList ans, temp;
  State temp_state;
  int LL, m, n;
  int w, h;
  int i;
  ans = NULL;
  i = 0;
  while (i < G->arc_index) {
    temp_state[i] = incoming[i];
    i++;
  }
  LL = 0;
  while (LL < G->arc_index) {
    w = 1;
    h = min(mod_up(incoming[LL] - G->Os[LL], G), mod_up(incoming[LL] - G->Xs[LL], G));
    while (w < G->arc_index && h > 0) {
      if (mod_up(incoming[LL] - incoming[mod(LL + w, G)], G) < h) {
        temp_state[LL] = incoming[mod(LL + w, G)];
        temp_state[mod(LL + w, G)] = incoming[LL];
        m = get_number(temp_state, prevs);
        if (m == 0) {
          n = get_number(temp_state, ans);
          if (n != 0) {
            ans = remove_state(temp_state, ans);
          } else {
            temp = swap_cols(LL, mod(LL + w, G), incoming, G);
            temp->nextState = ans;
            ans = temp;
          }
        };
        temp_state[LL] = incoming[LL];
        temp_state[mod(LL + w, G)] = incoming[mod(LL + w, G)];
        h = mod_up(incoming[LL] - incoming[mod(LL + w, G)], G);
      };
      h = min(h, min(mod_up(incoming[LL] - G->Os[mod(LL + w, G)], G),
                     mod_up(incoming[LL] - G->Xs[mod(LL + w, G)], G)));
      w++;
    };
    LL++;
  };
  return ans;
}

/**
 * Returns a single element statelist containing
 * the permutation incoming with x1 and x2 swapped
 * @param x1 an integer between 0 and MAX_INDEX
 * @param x2 an integer between 0 and MAX_INDEX
 * @param incoming a permutation
 * @param G working grid
 * @return a single element StateList containing a copy of
 * the permutation incoming with entries x1 and x2 swapped
 */
StateList swap_cols(int x1, int x2, char *incoming, Grid_t* G) {
  StateList ans;
  int i;
  i = 0;
  ans = malloc(sizeof(StateNode_t));
  i = 0;
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
 * Prints states in the form "{<state>,...}" up
 * to the first 500,000 states.
 * @param states a StateList
 * @param G working grid
 * @see print_state_short
 */
void print_states(StateList states, Grid_t* G) {
  StateList temp;
  int c;
  temp = states;
  printf("{");
  c = 0;
  while ((temp != NULL) && c < 500000) {
    print_state_short(temp->data, G);
    temp = temp->nextState;
    if (temp != NULL) {
      printf(",");
    };
    c++;
  };
  if (c == 500000) {
    printf("...");
  };
  printf("}");
}

/**
 * Prints the permutation of state using one line notation
 * "{_,_,...}"
 * @param state a State
 * @param G working grid
 */
void print_state_short(State state, Grid_t* G) {
  int i;
  i = 0;
  printf("{");
  while (i < G->arc_index - 1) {
    printf("%d,", state[i]);
    i++;
  };
  printf("%d}", state[G->arc_index - 1]);
}

/**
 * Prints the permutation of state on the grid specified by
 * Xs and Os, as well as 2A=M=SL+1.
 * @param state a State
 * @param G working grid
 */
void print_state(State state, Grid_t* G) {
  int i, j;
  j = G->arc_index;
  while (j > 0) {
    i = 0;
    while (i < G->arc_index) {
      if (G->Xs[i] == j) {
        printf("  X  ");
      } else {
        if (G->Os[i] == j) {
          printf("  O  ");
        } else {
          printf("  -  ");
        };
      };
      i++;
    };
    printf("\n");
    i = 0;
    while (i < G->arc_index) {
      if (state[i] == j) {
        printf("*    ");
      } else {
        printf("     ");
      };
      i++;
    };
    printf("\n");
    j--;
  };
  printf("\n");
  printf("2A=M=SL+1=%d\n",
         NESW_pp(G->Xs, G) - NESW_pO(G->Xs, G) - NESW_Op(G->Xs, G) + NESW_pp(G->Os, G) + 1);
}

/**
 * Returns the index of a State within a StateList
 * Note: Indexed from 1
 * @param a State
 * @param b a StateList
 * @return index of a within b
 */
int get_number(State a, StateList b) {
  StateList temp;
  int count = 1;
  temp = b;
  while (temp != NULL) {
    if (eq_state(a, temp->data)) {
      return count;
    };
    temp = temp->nextState;
    count++;
  };
  return 0;
}

/**
 * Takes in a list of parent vertices and a list of child vertices,
 * generates all edges between them and adds this EdgeList to global_edge_list
 * mod 2
 * @param parents a list of parent vertices
 * @param kids a list of child vertices
 * @param edge_list list of Edges between vertices
 * @return A shortEdges containing the result of adding them mod two.
 */
EdgeList add_mod_two_lists(VertexList parents, VertexList kids, EdgeList edge_list) {
  VertexList this_kid, this_parent, temp_vert;
  EdgeList this_edge, temp_edge, prev;
  EdgeList ans;
  if ((parents == NULL) || (kids == NULL)) {
    ans = edge_list; // change to return
  } else {
    this_parent = parents;
    this_kid = kids;
    ans = NULL;
    this_edge = edge_list;
    while (this_parent != NULL && this_edge != NULL &&
           (this_edge->start == this_parent->data &&
            this_edge->end == this_kid->data)) {
      temp_edge = this_edge; // This may need to be freed
      this_edge = this_edge->nextEdge;
      this_kid = this_kid->nextVertex;
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
 * Places a new edge into the supplied edgelist in order.
 * @param a an int specifying the parent vertex
 * @param b an int specifying the child vertex
 * @param edges a EdgeList where the new edge will be placed
 * @return a pointer to edges with the new edge added
 */
EdgeList append_ordered(int a, int b, EdgeList edges) {
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
 * Add a new edge to the start of a EdgeList
 * @param a an int indictating the source of the edge
 * @param b an int indictating the destination of the edge
 * @param e the EdgeList
 * @return e with the edge (a,b) at the front
 */
EdgeList prepend_edge(int a, int b, EdgeList e) {
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
EdgeList create_edge(int a, int b) {
  EdgeList new_ptr;
  new_ptr = malloc(sizeof(EdgeNode_t));
  new_ptr->start = a;
  new_ptr->end = b;
  new_ptr->nextEdge = NULL;
  return (new_ptr);
}

/**
 * Creates a vertex and adds it to the start of the supplied
 * VertexList
 * @param a an int specifying the vertex
 * @param vertices a VertexList
 * @return a VertexList with a Vertex containing a prepended to vertices
 */
VertexList prepend_vertex(int a, VertexList vertices) {
  VertexList new_ptr;
  new_ptr = malloc(sizeof(Vertex_t));
  (new_ptr->data) = a;
  (new_ptr->nextVertex) = vertices;
  return new_ptr;
}

/**
 * contracts all edges such that the parents occur after init
 * and the children are before or at final.
 * @param init an int specifying the required start
 * @param final State to t
 * @param edge_list the EdgeList
 */
void special_homology(int init, int final, EdgeList edge_list) {
  int i, j;
  EdgeList temp;
  i = 0;
  j = 0;
  temp = edge_list;
  while ((edge_list != NULL) && (temp != NULL)) {
    while ((temp != NULL) && (temp->start == init)) {
      temp = temp->nextEdge;
    };
    while ((temp != NULL) && (temp->end > final)) {
      temp = temp->nextEdge;
    };
    if (verbose && j == 100) {
      j = 0;
      if (temp != NULL)
        printf("Iteration number %d; contracting edge starting at (%d,%d)\n", i,
               temp->start, temp->end);
    };
    i++;
    j++;
    if (temp != NULL) {
      contract(temp->start, temp->end, edge_list);
      temp = edge_list;
    };
  };
}

/**
 * contracts the edge specified by the input within global_edge_list
 * @param a the parent vertex of the edge
 * @param b the child vertex of the edge
 * @param edge_list the EdgeList
 */
void contract(int a, int b, EdgeList edge_list) {
  EdgeList temp;
  EdgeList prev;
  VertexList parents, kids, temp_kids, temp_parents;
  VertexList last_parent, last_kid;
  prev = edge_list; // Multiple equal initializations
  parents = NULL;
  kids = NULL;
  temp_kids = NULL;
  temp_parents = NULL;
  last_parent = NULL;
  last_kid = NULL;
  while (edge_list != NULL &&
         ((edge_list)->end == b || edge_list->start == a)) {
    if ((edge_list->end == b) && (edge_list->start == a)) {
      temp = edge_list;
      edge_list = edge_list->nextEdge;
      free(temp);
    } else {
      if (edge_list->end == b) {
        if (last_parent == NULL) {
          parents = prepend_vertex(edge_list->start, NULL);
          last_parent = parents;
        } else {
          temp_parents = prepend_vertex(edge_list->start, NULL);
          last_parent->nextVertex = temp_parents;
          last_parent = last_parent->nextVertex;
        };
        temp = edge_list;
        edge_list = edge_list->nextEdge;
        free(temp);
      } else if (edge_list->start == a) {
        if (last_kid == NULL) {
          kids = prepend_vertex(edge_list->end, kids);
          last_kid = kids;
        } else {
          temp_kids = prepend_vertex(edge_list->end, NULL);
          last_kid->nextVertex = temp_kids;
          last_kid = last_kid->nextVertex;
        };
        temp = edge_list;
        edge_list = edge_list->nextEdge;
        free(temp);
      };
    };
  };
  prev = edge_list;
  if (edge_list != NULL) {
    temp = (edge_list->nextEdge);
  } else {
    temp = NULL;
  };
  while (temp != NULL && (temp)->start < a) {
    if ((temp)->end == b) {
      if (last_parent == NULL) {
        parents = prepend_vertex((temp)->start, NULL);
        last_parent = parents;
      } else {
        temp_parents = prepend_vertex((temp)->start, NULL);
        last_parent->nextVertex = temp_parents;
        last_parent = last_parent->nextVertex;
      };
      (prev->nextEdge) = (temp->nextEdge);
      free(temp);
      temp = prev->nextEdge;
    } else {
      temp = (temp)->nextEdge;
      prev = (prev)->nextEdge;
    };
  };
  while (temp != NULL && (temp)->start == a) {
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
    (prev)->nextEdge = temp->nextEdge;
    free(temp);
    temp = (prev)->nextEdge;
  };
  while (temp != NULL) {
    if ((temp)->end == b) {
      if (last_parent == NULL) {
        parents = prepend_vertex(temp->start, NULL);
        last_parent = parents;
      } else {
        temp_parents = prepend_vertex((temp)->start, NULL);
        last_parent->nextVertex = temp_parents;
        last_parent = last_parent->nextVertex;
      };
      (prev)->nextEdge = (temp)->nextEdge;
      free(temp);
      temp = (prev)->nextEdge;
    } else {
      temp = (temp)->nextEdge;
      prev = (prev)->nextEdge;
    }
  };
  edge_list = add_mod_two_lists(parents, kids, edge_list);
}

/**
 * Removes a state from a StateList. Equality checked using eq_state.
 * @param a a State
 * @param v a StateList
 * @return a StateList removing a from v
 * @see eq_state
 */
StateList remove_state(State a, StateList v) {
  StateList temp, prev;
  StateList s_list = v;
  prev = v;
  if (v == NULL)
    return (NULL);
  else if (eq_state(a, v->data)) {
    temp = v;
    s_list = v->nextState;
    free(v);
    return (s_list);
  } else {
    temp = prev->nextState;
    while ((temp != NULL) && (!eq_state(a, temp->data))) {
      temp = temp->nextState;
      prev = prev->nextState;
    };
    if (temp != NULL) {
      prev->nextState = temp->nextState;
      free(temp);
      return s_list;
    } else
      return s_list;
  }
}

/**
 * Prints each edge in global_edge_list on a new line
 * @param edge_list an EdgeList
 */
void print_edges(EdgeList edge_list) {
  EdgeList temp;
  temp = edge_list;
  while (temp != NULL) {
    printf("%d %d\n", temp->start, temp->end);
    temp = (temp->nextEdge);
  };
}

/**
 * Print the first 80 edges in global_edge_list on the same line
 * @param edge_list an EdgeList
 */
void print_math_edges(EdgeList edge_list) {
  EdgeList temp;
  int t;
  temp = edge_list;
  printf("{");
  t = 0;
  while (temp != NULL) {
    printf("{%d,%d}", temp->start, temp->end);
    t++;
    if (t == 80) {
      temp = NULL;
      printf("...}\n");
    } else {
      temp = (temp->nextEdge);
      if (temp != NULL)
        printf(",");
    }
  };
  printf("}");
}

/**
 * Prints the edges in edges on a single line
 * @param edges an EdgeList
 */
void print_math_edges_a(EdgeList edges) {
  EdgeList temp;
  temp = edges;
  printf("{");
  while (temp != NULL) {
    printf("{%d,%d}", temp->start, temp->end);
    temp = (temp->nextEdge);
    if (temp != NULL)
      printf(",");
  };
  printf("}");
}

/**
 * Prints the vertices in VertexList on a single line
 * @param v_list a VertexList
 */
void print_vertices(VertexList v_list) {
  VertexList temp;
  temp = v_list;
  printf("{");
  while (temp != NULL) {
    printf("%d", (temp)->data);
    temp = (temp)->nextVertex;
    if (temp != NULL)
      printf(",");
  };
  printf("}");
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
    free(temp);
    temp = states;
  };
}

/**
 * Frees the supplied EdgeList
 * @param e a EdgeList
 */
void free_edge_list(EdgeList e) {
  EdgeList temp, ntemp;
  temp = e;
  ntemp = temp;
  while (temp != NULL) {
    ntemp = temp;
    temp = temp->nextEdge;
    free(ntemp);
  };
}

/* Higher differentials */

/**
 * Calculates a StateList containing all states reachable
 * by a rectangles of a fixed width
 * @param wt an int specifying rectangle width
 * @param incoming origin state for the rectangles
 * @param G working grid
 * @return a StateList with states that are reached by a rectangle of width
 * wt from incoming.
 */
StateList fixed_wt_rectangles_out_of(int wt, State incoming, Grid_t* G) {
  StateList temp, ans;
  int LL;
  int w, h;
  int this_weight, i;
  ans = NULL;
  LL = 0;
  while (LL < G->arc_index) {
    w = 1;
    h = mod(G->Os[LL] - incoming[LL], G);
    while (w < G->arc_index && h > 0) {
      if (mod(incoming[mod(LL + w, G)] - incoming[LL], G) <= h) {
        this_weight = 0;
        i = 0;
        while (i < w && this_weight <= wt + 1) {
          if (mod(G->Xs[mod(LL + i, G)] - incoming[LL], G) <
              mod(incoming[mod(LL + w, G)] - incoming[LL], G)) {
            this_weight++;
          };
          i++;
        };
        if (this_weight == wt) {
          temp = swap_cols(LL, mod(LL + w, G), incoming, G);
          if (get_number(temp->data, ans) != 0) {
            ans = remove_state(temp->data, ans);
            free(temp);
            temp = ans;
          } else {
            temp->nextState = ans;
            ans = temp;
          }
        };
        h = mod(incoming[mod(LL + w, G)] - incoming[LL], G);
      };
      h = min(h, mod(G->Os[mod(LL + w, G)] - incoming[LL], G));
      w++;
    };
    LL++;
  };
  return ans;
}

/**
 * Calculates whether the supplied state is nullhomologous. Uses
 the global variable global_edge_list.
 * @param init a State
 * @param G working grid
 * @param edge_list the edge_list
 * @return nonzero if nullhomologous and zero otherwise.
 */
int null_homologous_D0Q(State init, Grid_t* G, EdgeList edge_list) {
  StateList new_ins, new_outs, last_new_in, last_new_out, temp;
  StateList prev_ins, prev_outs;
  StateList really_new_outs = NULL, really_new_ins = NULL;
  int in_number, ans, prev_in_number;
  int out_number;
  int i;
  int edge_count = 0;
  int num_ins = 0;
  int num_outs = 0;
  int num_new_ins = 0;
  int num_new_outs = 0;
  StateList present_in, present_out;
  edge_list = prepend_edge(0, 1, NULL);
  prev_outs = NULL;
  prev_ins = NULL;
  new_ins = malloc(sizeof(StateNode_t));
  i = 0;
  while (i < G->arc_index) {
    new_ins->data[i] = init[i];
    i++;
  };
  new_ins->nextState = NULL;
  ans = 0;
  while (new_ins != NULL && !ans) {
    present_in = new_ins;
    in_number = 0;
    num_new_outs = 0;
    new_outs = NULL;
    while (present_in != NULL) {
      free_state_list(really_new_outs);
      in_number++;
      really_new_outs = new_rectangles_into(prev_outs, present_in->data, G);
      while (really_new_outs != NULL) {
        out_number = get_number(really_new_outs->data, new_outs);
        if (out_number == 0) {
          if (num_new_outs == 0) {
            new_outs = really_new_outs;
            really_new_outs = really_new_outs->nextState;
            new_outs->nextState = NULL;
            last_new_out = new_outs;
            num_new_outs++;
            out_number = num_new_outs;
          } else {
            last_new_out->nextState = really_new_outs;
            really_new_outs = really_new_outs->nextState;
            last_new_out = last_new_out->nextState;
            last_new_out->nextState = NULL;
            num_new_outs++;
            out_number = num_new_outs;
          };
        } else {
          temp = really_new_outs;
          really_new_outs = really_new_outs->nextState;
          free(temp);
        }
        edge_list = append_ordered(
            out_number + num_outs, in_number + num_ins, edge_list);
        edge_count++;
      }
      present_in = present_in->nextState;
    };
    free_state_list(prev_ins);
    prev_ins = new_ins;
    i = 1;
    num_ins = num_ins + in_number;
    prev_in_number = num_ins;
    num_new_ins = 0;
    new_ins = NULL;
    out_number = 0;
    present_out = new_outs;
    while (present_out != NULL) {
      out_number++;
      really_new_ins = new_rectangles_out_of(prev_ins, present_out->data, G);
      while (really_new_ins != NULL) {
        in_number = get_number(really_new_ins->data, new_ins);
        if (in_number == 0) {
          if (num_new_ins == 0) {
            new_ins = really_new_ins;
            really_new_ins = really_new_ins->nextState;
            new_ins->nextState = NULL;
            last_new_in = new_ins;
            num_new_ins++;
            in_number = num_new_ins;
          } else {
            last_new_in->nextState = really_new_ins;
            really_new_ins = really_new_ins->nextState;
            last_new_in = last_new_in->nextState;
            last_new_in->nextState = NULL;
            num_new_ins++;
            in_number = num_new_ins;
          };
        } else {
          temp = really_new_ins;
          really_new_ins = really_new_ins->nextState;
          free(temp);
        }
        edge_list = append_ordered(
            out_number + num_outs, in_number + num_ins, edge_list);
        edge_count++;
      };
      present_out = present_out->nextState;
    };
    free_state_list(prev_outs);
    prev_outs = new_outs;
    new_outs = NULL;
    special_homology(0, prev_in_number, edge_list);
    if ((edge_list == NULL) || (edge_list->start != 0)) {
      ans = 1;
      free_state_list(new_ins);
      free_state_list(new_outs);
      new_ins = NULL;
    } else if (edge_list->end <= prev_in_number) {
      ans = 0;
      free_state_list(new_ins);
      free_state_list(new_outs);
      new_ins = NULL;
    } else {
      num_outs = num_outs + out_number;
      if (verbose) {
        printf("%d %d %d\n", num_ins, num_outs, edge_count);
      }
    };
  };
  return (ans);
}

/**
 * Calculates if D1 of the supplied state is nullhomologous
 * @param init a State
 * @param G working grid
 * @param edge_list the EdgeList
 * @return nonzero if nullhomologous and zero otherwise
 */
int null_homologous_D1Q(State init, Grid_t* G, EdgeList edge_list) {
  StateList new_ins, new_outs, last_new_in, last_new_out, temp;
  StateList prev_ins, prev_outs;
  StateList really_new_outs = NULL, really_new_ins = NULL;
  EdgeList LastEdge;
  int in_number, ans, prev_in_number;
  int out_number;
  int i;
  int edge_count = 0;
  int num_ins = 0;
  int num_outs = 0;
  int num_new_ins = 0;
  int num_new_outs = 0;
  StateList present_in, present_out;
  LastEdge = edge_list;
  prev_outs = NULL;
  prev_ins = NULL;
  new_ins = fixed_wt_rectangles_out_of(1, init, G);
  edge_list = prepend_edge(0, 1, NULL);
  temp = new_ins;
  i = 1;
  LastEdge = NULL;
  if (temp != NULL) {
    i = 1;
    edge_list = create_edge(0, 1);
    LastEdge = edge_list;
    temp = temp->nextState;
    while (temp != NULL) {
      i++;
      LastEdge->nextEdge = create_edge(0, i);
      LastEdge = LastEdge->nextEdge;
      temp = temp->nextState;
    };
  };
  ans = 0;
  while (new_ins != NULL && !ans) {
    present_in = new_ins;
    in_number = 0;
    num_new_outs = 0;
    new_outs = NULL;
    while (present_in != NULL) {
      free_state_list(really_new_outs);
      in_number++;
      really_new_outs = new_rectangles_into(prev_outs, present_in->data, G);
      while (really_new_outs != NULL) {
        out_number = get_number(really_new_outs->data, new_outs);
        if (out_number == 0) {
          if (num_new_outs == 0) {
            new_outs = really_new_outs;
            really_new_outs = really_new_outs->nextState;
            new_outs->nextState = NULL;
            last_new_out = new_outs;
            num_new_outs++;
            out_number = num_new_outs;
          } else {
            last_new_out->nextState = really_new_outs;
            really_new_outs = really_new_outs->nextState;
            last_new_out = last_new_out->nextState;
            last_new_out->nextState = NULL;
            num_new_outs++;
            out_number = num_new_outs;
          };
        } else {
          temp = really_new_outs;
          really_new_outs = really_new_outs->nextState;
          free(temp);
        }
        edge_list = append_ordered(
                                   out_number + num_outs, in_number + num_ins, edge_list);
        edge_count++;
      }
      present_in = present_in->nextState;
    };
    free_state_list(prev_ins);
    prev_ins = new_ins;
    i = 1;
    num_ins = num_ins + in_number;
    prev_in_number = num_ins;
    num_new_ins = 0;
    new_ins = NULL;
    out_number = 0;
    present_out = new_outs;
    while (present_out != NULL) {
      out_number++;
      really_new_ins = new_rectangles_out_of(prev_ins, present_out->data, G);
      while (really_new_ins != NULL) {
        in_number = get_number(really_new_ins->data, new_ins);
        if (in_number == 0) {
          if (num_new_ins == 0) {
            new_ins = really_new_ins;
            really_new_ins = really_new_ins->nextState;
            new_ins->nextState = NULL;
            last_new_in = new_ins;
            num_new_ins++;
            in_number = num_new_ins;
          } else {
            last_new_in->nextState = really_new_ins;
            really_new_ins = really_new_ins->nextState;
            last_new_in = last_new_in->nextState;
            last_new_in->nextState = NULL;
            num_new_ins++;
            in_number = num_new_ins;
          };
        } else {
          temp = really_new_ins;
          really_new_ins = really_new_ins->nextState;
          free(temp);
        }
        edge_list = append_ordered(
                                   out_number + num_outs, in_number + num_ins, edge_list);
        edge_count++;
      };
      present_out = present_out->nextState;
    };
    free_state_list(prev_outs);
    prev_outs = new_outs;
    new_outs = NULL;
    special_homology(0, prev_in_number, edge_list);
    if ((edge_list == NULL) || (edge_list->start != 0)) {
      ans = 1;
      free_state_list(new_ins);
      free_state_list(new_outs);
      new_ins = NULL;
    } else if (edge_list->end <= prev_in_number) {
      ans = 0;
      free_state_list(new_ins);
      free_state_list(new_outs);
      new_ins = NULL;
    } else {
      num_outs = num_outs + out_number;
      if (verbose) {
        printf("%d %d %d\n", num_ins, num_outs, edge_count);
      }
    };
  };
  return (ans);
}

/**
 * Sum over each point in the permutation count the number of Os
 * that occur to the northeast
 * @param x a permutation
 * @param G working Grid
 * @return an int containing the quantity described above
 */
int NESW_pO(char *x, Grid_t* G) {
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
int NESW_Op(char *x, Grid_t* G) {
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
 * @param G working gridn
 * @return an int containing the quantity described above
 */
int NESW_pp(char *x, Grid_t* G) {
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
