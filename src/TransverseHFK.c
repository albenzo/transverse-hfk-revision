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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

const char *argp_program_version = "transverseHFK revision 0.0.1";
const char *argp_program_bug_address = "<lmeye22@lsu.edu>";
static const char doc[] =
    "A program to calculate the Legendrian/Transverse knot invariants\
 via the algorithm described in \"Transverse knots distinguished by\
 Knot Floer Homology\" by L. Ng, P. S. Ozsvath, and D. P. Thurston.\
 If the number of sheets is not equal to 1 it instead calculates the\
 theta invariant for the n-fold cyclic cover.";

static const char args_doc[] = "-i [ArcIndex] -X [Xs] -O [Os]";

static struct argp_option options[] = {
    {"verbose", 'v', 0, 0, "Produce verbose output", 0},
    {"quiet", 'q', 0, 0, "Produce some extraneous output", 0},
    {"silent", 's', 0, 0, "Don't produce any extraneous output", 0},
    {"index", 'i', "ArcIndex", 0, "ArcIndex of the grid", 0},
    {"Xs", 'X', "[...]", 0, "List of Xs", 0},
    {"Os", 'O', "[...]", 0, "List of Os", 0},
    {"sheets", 'n', 0, 0, "Number of sheets for cyclic branch cover. Default: 1", 0},
    {"timeout", 't', "SECONDS", 0, "Maximum time to run in seconds", 0},
    {0}};

// Temporary location
static error_t parse_opt(int, char *, struct argp_state *);
static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};
struct arguments {
  int arc_index;
  int sheets;
  char *Xs;
  char *Os;
  int max_time;
};

typedef int (*printf_t)(const char *format, ...);
#define SILENT 0
#define QUIET 1
#define VERBOSE 2

static int verbosity = SILENT;
printf_t print_ptr = printf;

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

EdgeList add_mod_two_lists(const VertexList, const VertexList,
                           const EdgeList *const);
VertexList prepend_vertex(const int, const VertexList);
EdgeList prepend_edge(const int, const int, const EdgeList);
StateList fixed_wt_rectangles_out_of(const int, const State,
                                     const Grid_t *const);
StateList swap_cols(const int, const int, const State, const Grid_t *const);
int get_number(const State, const StateList, const Grid_t *const);
int get_lift_number(const LiftState, LiftStateList, const LiftGrid_t * const);
void free_state_list(StateList);
int null_homologous_D0Q(const State, const Grid_t *const);
int null_homologous_D1Q(const State, const Grid_t *const);
int null_homologous_lift(const LiftState, const LiftGrid_t * const);
void special_homology(const int, const int, EdgeList *);
void contract(const int, const int, EdgeList *);
EdgeList create_edge(const int, const int);
void free_edge_list(const EdgeList);
StateList remove_state(const State, const StateList, const Grid_t *const);
LiftStateList remove_lift_state(const LiftState, const LiftStateList, const LiftGrid_t * const);
int mod(const int, const Grid_t *const);
int mod_up(const int, const Grid_t *const);
StateList new_rectangles_out_of(const StateList, const State,
                                const Grid_t *const);
StateList new_rectangles_into(const StateList, const State,
                              const Grid_t *const);
static LiftStateList new_lift_rectangles_out_internal(const LiftStateList, const LiftState, const LiftGrid_t * const, int);
LiftStateList new_lift_rectangles_out_of(const LiftStateList, const LiftState, const LiftGrid_t * const);
LiftStateList new_lift_rectangles_into(const LiftStateList, const LiftState, const LiftGrid_t * const);
EdgeList append_ordered(const int, const int, const EdgeList);

int NESW_pO(const State, const Grid_t *const);
int NESW_Op(const State, const Grid_t *const);
int NESW_pp(const State, const Grid_t *const);
int writhe(const Grid_t *const);
void cusps(int *, const Grid_t *const);

void print_state(const State, const Grid_t *const);
void print_state_short(const State, const Grid_t *const);
void print_edges(const EdgeList);
void print_states(const StateList, const Grid_t *const);
void print_math_edges(const EdgeList);
void print_math_edges_a(const EdgeList);
void print_vertices(const VertexList);
void print_grid_perm(const Grid_t *const G);
void print_self_link(const Grid_t *const G);

int get_verbosity(void);
void set_verbosity(const int);

void timeout(const int);
int is_grid(const Grid_t *const);
int is_lift_grid(const LiftGrid_t * const);
int is_state(const State, const Grid_t *const);
void init_lift_state(LiftState *, const LiftGrid_t * const);
void copy_lift_state(LiftState *, const LiftState * const, const LiftGrid_t * const);
void mirror_lift_state(LiftState *, const LiftGrid_t * const);

int build_permutation(State, char *, int);
int eq_state(const State a, const State b, const Grid_t *const G) {
  return (!strncmp(a, b, G->arc_index));
}
int eq_lift_state(const LiftState a, const LiftState b, const LiftGrid_t *const G) {
  return (!strncmp((char*)a, (char*)b, G->arc_index*G->sheets));
}

int get_verbosity() { return verbosity; }

void set_verbosity(const int val) { verbosity = val; }

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  struct arguments *args = state->input;

  switch (key) {
  case 'v':
    set_verbosity(VERBOSE);
    break;
  case 'q':
    set_verbosity(QUIET);
    break;
  case 's':
    set_verbosity(SILENT);
    break;
  case 't':
    args->max_time = atoi(arg);
    if (args->max_time <= 0) {
      argp_failure(state, 0, 0, "Invalid timeout.");
      exit(1);
    }
    break;
  case 'i':
    args->arc_index = atoi(arg);
    if (args->arc_index < 2) {
      argp_failure(state, 0, 0, "ArcIndex must be a non-negative integer greater than 1.");
      exit(1);
    }
    break;
  case 'n':
    args->sheets = atoi(arg);
    if (args->sheets < 1) {
      argp_failure(state, 0, 0, "The number of sheets must be atleast 1.");
      exit(1);
    }
    break;
  case 'X': {
    args->Xs = arg;
  } break;
  case 'O': {
    args->Os = arg;
  } break;
  default:
    break;
  }
  return 0;
}

/**
 * Signal handler for SIGALRM that exits when the signal is received.
 * @param sig
 */
void timeout(const int sig) {
  if (SIGALRM == sig) {
    (*print_ptr)("Timeout reached. Terminating\n");
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
int build_permutation(char *perm, char *str, int len) {
  if (str[0] != '[') {
    return -1;
  }

  char *s = &str[1];
  long n = -1;
  int i = 0;

  while (i < len) {
    errno = 0;
    n = strtol(s, &s, 10);

    if ((errno == ERANGE && (n == LONG_MAX || n == LONG_MIN)) ||
        (errno != 0 && n == 0) || (n < 1 || n > len)) {
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

// This is not correct
void init_lift_state(LiftState *s, const LiftGrid_t * const G) {
  (**s)[G->arc_index] = malloc(sizeof(char)*G->sheets*G->arc_index);
}

/**
 * Determines whether the grid specified by the parameters is valid
 * @param G working grid
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

int main(int argc, char **argv) {
  struct arguments args;
  args.arc_index = -1;
  args.sheets = 1;
  args.max_time = -1;
  args.Xs = NULL;
  args.Os = NULL;
  argp_parse(&argp, argc, argv, 0, 0, &args);

  if (args.arc_index == -1) {
    fprintf(stderr, "transverseHFK: Missing arc_index\n");
    exit(1);
  }

  if (args.Xs == NULL) {
    fprintf(stderr, "transverseHFK: Missing Xs\n");
    exit(1);
  }

  if (args.Os == NULL) {
    fprintf(stderr, "transverseHFK: Missing Os\n");
    exit(1);
  }

  if (args.sheets > 1) {
    // Do other invariant
    LiftGrid_t G;
    G.arc_index = args.arc_index;
    G.sheets = args.sheets;
    G.Xs = malloc(sizeof(char) * G.arc_index);
    G.Os = malloc(sizeof(char) * G.arc_index);

    if (-1 == build_permutation(G.Xs, args.Xs, args.arc_index)) {
      fprintf(stderr, "transverseHFK: Malformatted Xs\n");
      free(G.Xs);
      free(G.Os);
      exit(1);
    }

    if (-1 == build_permutation(G.Os, args.Os, args.arc_index)) {
      fprintf(stderr, "transverseHFK: Malformatted Os\n");
      free(G.Xs);
      free(G.Os);
      exit(1);
    }

    if(!is_lift_grid(&G)) {
      (*print_ptr)("Invalid grid\n");
      exit(1);
    }

    if (args.max_time > 0) {
      if (signal(SIGALRM, timeout) == SIG_ERR) {
        perror("An error occured while setting the timer");
        free(G.Xs);
        free(G.Os);
        exit(1);
      }
      alarm(args.max_time);
    }

    LiftState UR_lift;
    init_lift_state(&UR_lift, &G);
    
    for (int j = 0; j < G.sheets; ++j) {
      if (G.Xs[G.arc_index - 1] == G.arc_index) {
        UR_lift[j][0] = 1;
      } else {
        UR_lift[j][0] = (char)G.Xs[G.arc_index - 1] + 1;
      };
      for (int i = 1; i < G.arc_index; ++i) {
        if (G.Xs[i - 1] == G.arc_index) {
          UR_lift[j][i] = 1;
        } else {
          UR_lift[j][i] = G.Xs[i - 1] + 1;
        }
      }
    }

    if (QUIET <= get_verbosity()) {
      (*print_ptr)("Calculating graph for lifted invariant.\n");
      // These print statements are wrong
      //print_state(G.Xs, &G);
      //print_self_link(&G);
    }

    // Change the content of these print statements
    if (null_homologous_lift(UR_lift, &G)) {
      (*print_ptr)("%d: is null-homologous\n");
    }
    else {
      (*print_ptr)("%d: is NOT null-homologous\n");
    }

    free(G.Xs);
    free(G.Os);
    free(UR_lift);
    
    exit(0);
  }

  Grid_t G;
  G.arc_index = args.arc_index;
  G.Xs = malloc(sizeof(char) * G.arc_index);
  G.Os = malloc(sizeof(char) * G.arc_index);

  if (-1 == build_permutation(G.Xs, args.Xs, args.arc_index)) {
    fprintf(stderr, "transverseHFK: Malformatted Xs\n");
    free(G.Xs);
    free(G.Os);
    exit(1);
  }

  if (-1 == build_permutation(G.Os, args.Os, args.arc_index)) {
    fprintf(stderr, "transverseHFK: Malformatted Os\n");
    free(G.Xs);
    free(G.Os);
    exit(1);
  }

  State UR = malloc(sizeof(char) * G.arc_index);
  int i;

  if (!is_grid(&G)) {
    (*print_ptr)("Invalid grid\n");
    exit(1);
  }

  if (args.max_time > 0) {
    if (signal(SIGALRM, timeout) == SIG_ERR) {
      perror("An error occured while setting the timer");
      free(G.Xs);
      free(G.Os);
      exit(1);
    }
    alarm(args.max_time);
  }

  if (QUIET <= get_verbosity()) {
    (*print_ptr)("\n \nCalculating graph for LL invariant\n");
    print_state(G.Xs,&G);
    print_self_link(&G);
  }
  if (null_homologous_D0Q(G.Xs, &G)) {
    (*print_ptr)("LL is null-homologous\n");
  } else {
    (*print_ptr)("LL is NOT null-homologous\n");
  }

  if (QUIET <= get_verbosity()) {
    (*print_ptr)("\n \nCalculating graph for UR invariant\n");
  }
  if (G.Xs[G.arc_index - 1] == G.arc_index) {
    UR[0] = 1;
  } else {
    UR[0] = (char)G.Xs[G.arc_index - 1] + 1;
  };
  i = 1;
  while (i < G.arc_index) {
    if (G.Xs[i - 1] == G.arc_index) {
      UR[i] = 1;
    } else {
      UR[i] = G.Xs[i - 1] + 1;
    }
    ++i;
  }
  if (QUIET <= get_verbosity()) {
    print_state(UR, &G);
    print_self_link(&G);
  }

  if (null_homologous_D0Q(UR, &G)) {
    (*print_ptr)("UR is null-homologous\n");
  } else {
    (*print_ptr)("UR is NOT null-homologous\n");
  };

  if (QUIET <= get_verbosity()) {
    (*print_ptr)("\n \nCalculating graph for D1[LL] invariant\n");
    print_state(G.Xs, &G);
    print_self_link(&G);
  }

  if (null_homologous_D1Q(G.Xs, &G)) {
    (*print_ptr)("D1[LL] is null-homologous\n");
  } else {
    (*print_ptr)("D1[LL] is NOT null-homologous\n");
  }

  if (QUIET <= get_verbosity()) {
    (*print_ptr)("\n \nCalculating graph for D1[UR] invariant\n");
  }

  if (G.Xs[G.arc_index - 1] == G.arc_index) {
    UR[0] = 1;
  } else {
    UR[0] = (char)G.Xs[G.arc_index - 1] + 1;
  };
  i = 1;
  while (i <= G.arc_index - 1) {
    if (G.Xs[i - 1] == G.arc_index) {
      UR[i] = 1;
    } else {
      UR[i] = G.Xs[i - 1] + 1;
    }
    i++;
  }

  if (QUIET <= get_verbosity()) {
    print_state(UR, &G);
    print_self_link(&G);
  }

  if (null_homologous_D1Q(UR, &G)) {
    (*print_ptr)("D1[UR] is null-homologous\n");
  } else {
    (*print_ptr)("D1[UR] is NOT null-homologous\n");
  };

  free(G.Xs);
  free(G.Os);
  free(UR);

  return 0;
}

/**
 * Shifts the input towards the interval [0,arc_index) by
 * a multiple of arc_index
 * @param a An integer
 * @param G working grid
 * @return a shifted towards the interval [0,arc_index)
 */
int mod(const int a, const Grid_t *const G) {
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
int mod_up(const int a, const Grid_t *const G) {
  if (a > G->arc_index) {
    return (a - G->arc_index);
  } else if (a <= 0) {
    return (a + G->arc_index);
  } else {
    return (a);
  };
}

int min(const int a, const int b) {
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
StateList new_rectangles_out_of(const StateList prevs, const State incoming,
                                const Grid_t *const G) {
  StateList temp, ans;
  State temp_state = malloc(sizeof(char) * G->arc_index);
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
        m = get_number(temp_state, prevs, G);
        if (m == 0) {
          n = get_number(temp_state, ans, G);
          if (n != 0) {
            ans = remove_state(temp_state, ans, G);
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
StateList new_rectangles_into(const StateList prevs, const State incoming,
                              const Grid_t *const G) {
  StateList ans, temp;
  State temp_state = malloc(sizeof(char) * G->arc_index);
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
    h = min(mod_up(incoming[LL] - G->Os[LL], G),
            mod_up(incoming[LL] - G->Xs[LL], G));
    while (w < G->arc_index && h > 0) {
      if (mod_up(incoming[LL] - incoming[mod(LL + w, G)], G) < h) {
        temp_state[LL] = incoming[mod(LL + w, G)];
        temp_state[mod(LL + w, G)] = incoming[LL];
        m = get_number(temp_state, prevs, G);
        if (m == 0) {
          n = get_number(temp_state, ans, G);
          if (n != 0) {
            ans = remove_state(temp_state, ans, G);
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

static LiftStateList new_lift_rectangles_out_internal(const LiftStateList prevs, const LiftState incoming, const LiftGrid_t *const G, int is_mirrored) {
  LiftStateList ans = NULL;
  double *g_Xs, *g_Os;
  LiftState temp_state;

  g_Xs = malloc(sizeof(double)*G->arc_index);
  g_Os = malloc(sizeof(double)*G->arc_index);

  for(int i=0; i < G->arc_index; ++i) {
    g_Xs[i] = ((double)G->Xs[i]) - .5;
    g_Os[i] = ((double)G->Os[i]) - .5;
  }

  init_lift_state(&temp_state, G);

  for(int start_x=0; start_x < G->arc_index*G->sheets; ++start_x) {
    int jumped_down = 0;
    int jumped_up = 0;
    char start_y = incoming[start_x/G->arc_index][start_x%G->arc_index];
    int start_sheet = start_x/G->arc_index;
    int step = 0;
    int check_index = (start_x + step) % G->arc_index;
    int jump = start_sheet*G->arc_index;
    int height = (start_y - 1) % G->arc_index;

    while (height != start_y) {
      check_index = (start_x + step)%G->arc_index;
      int check_index_gen = (((start_x + step + 1) % G->arc_index) + jump) % (G->arc_index * G->sheets);
      int clear = 1;
      if (height > start_y) {
        if (g_Xs[check_index] < height && g_Xs[check_index] > start_y && clear) {
          clear = 0;
        }
        if (g_Os[check_index] < height && g_Os[check_index] > start_y && clear) {
          clear = 0;
        }
        if (g_Xs[check_index] > height && g_Os[check_index] < start_y && clear) {
          if (is_mirrored) {
            jump = jump - G->arc_index;
            jumped_down = 1;
          }
          else {
            jump = jump + G->arc_index;
            jumped_up = 1;
          }
          check_index_gen = (((start_x + step + 1) % G->arc_index) + jump) % (G->arc_index * G->sheets);
        }
        if (g_Os[check_index] > height && g_Xs[check_index] < start_y && clear) {
          if (is_mirrored) {
            jump = jump - G->arc_index;
            jumped_down = 1;
          }
          else {
            jump = jump + G->arc_index;
            jumped_up = 1;
          }
          check_index_gen = (((start_x + step + 1) % G->arc_index) + jump) % (G->arc_index * G->sheets);
        }
        if (incoming[check_index_gen/G->arc_index][check_index_gen%G->arc_index] < height && incoming[check_index_gen/G->arc_index][check_index_gen%G->arc_index] > start_y && clear) {
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
          check_index_gen = (((start_x + step + 1) % G->arc_index) + jump) % (G->arc_index * G->sheets);
          if (incoming[check_index_gen/G->arc_index][check_index_gen%G->arc_index] == height) {
            LiftState new_state = NULL;
            copy_lift_state(&new_state, &incoming, G);
            new_state[start_x/G->arc_index][start_x%G->arc_index] = incoming[check_index_gen/G->arc_index][check_index_gen%G->arc_index];
            new_state[check_index_gen/G->arc_index][check_index_gen%G->arc_index] = incoming[start_x/G->arc_index][start_x%G->arc_index];
            if (is_mirrored) {
              mirror_lift_state(&new_state, G);
            }
            if (0 == get_lift_number(new_state, prevs, G)) {
              if (0 == get_lift_number(new_state, ans, G)) {
                LiftStateNode_t* new_node = malloc(sizeof(LiftStateNode_t));
                new_node->data = new_state;
                new_node->nextState = ans;
                ans = new_node;
              }
              else {
                remove_lift_state(new_state, ans, G);
                free(new_state);
              }
              
              height = (height - 1) % G->arc_index;
          }
          step = step + 1;
          jumped_down = 0;
          jumped_up = 0;
        }
        else {
          height = (height - 1) % G->arc_index;
        }
      }
      else {
        if ((g_Xs[check_index] < height || g_Xs[check_index] > start_y) && clear) {
          clear = 0;
        }
        if ((g_Os[check_index] < height || g_Os[check_index] > start_y) && clear) {
          clear = 0;
        }
        if ((incoming[check_index_gen/G->arc_index][check_index_gen%G->arc_index] < height || incoming[check_index_gen/G->arc_index][check_index_gen%G->arc_index] >= start_y) && clear) {
          clear = 0;
        }
        if (clear) {
          if (incoming[check_index_gen/G->arc_index][check_index_gen%G->arc_index] == height) {
            LiftState new_state = NULL;
            init_lift_state(new_state, G);
            copy_lift_state(&new_state, &incoming, G);
            new_state[start_x/G->arc_index][start_x%G->arc_index] = incoming[check_index_gen/G->arc_index][check_index_gen%G->arc_index];
            new_state[check_index_gen/G->arc_index][check_index_gen%G->arc_index] = incoming[start_x/G->arc_index][start_x%G->arc_index];
            if (is_mirrored) {
              mirror_lift_state(&new_state, G);
            }
            if (0 == get_lift_number(new_state, prevs, G)) {
              if (0 == get_lift_number(new_state, ans, G)) {
                LiftStateNode_t* new_node = malloc(sizeof(LiftStateNode_t));
                new_node->data = new_state;
                new_node->nextState = ans;
                ans = new_node;
              }
              else {
                ans = remove_lift_state(new_state, ans, G);
                free(new_state);
              }
            height = (height -1) % G->arc_index;
          }
          step = step + 1;
        }
        else {
          height = (height - 1) % G->arc_index;
        }
      }
    }
  }


    }}
  free(g_Xs);
  free(g_Os);
  return ans;
}

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

void mirror_lift_state(LiftState * state, const LiftGrid_t * const G) {
  LiftState original;
  init_lift_state(&original, G);
  copy_lift_state(&original, state, G);
  for(int i = 0; i < G->sheets; ++i) {
    for(int j = 0; j < G->arc_index; ++j) {
      (*state)[j][i] = original[G->sheets-(j+1)][G->arc_index-(i+1)];
    }
  }
  free(original);
}

void copy_lift_state(LiftState* dest, const LiftState * const origin, const LiftGrid_t * const G) {
  for(int i = 0; i < G->sheets; ++i) {
    for(int j = 0; j < G->arc_index; ++j) {
      (*dest)[i][j] = (*origin)[i][j];
    }
  }
}

LiftStateList new_lift_rectangles_out_of(const LiftStateList prevs, const LiftState incoming, const LiftGrid_t *const G) {
  return new_lift_rectangles_out_internal(prevs, incoming, G, 0);
}

LiftStateList new_lift_rectangles_into(const LiftStateList prevs, const LiftState incoming, const LiftGrid_t *const G) {
  LiftGrid_t* G_mirror = mirror_lift_grid(G);
  LiftState incoming_mirror;
  init_lift_state(&incoming_mirror, G);
  copy_lift_state(&incoming_mirror, &incoming, G);
  mirror_lift_state(&incoming, G);
  LiftStateList ans = new_lift_rectangles_out_internal(prevs, incoming_mirror, G_mirror, 1);

  free(G_mirror);
  free(incoming_mirror);

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
StateList swap_cols(const int x1, const int x2, const State incoming,
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
  (*print_ptr)("%d}", state[G->arc_index - 1]);
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
  print_grid_perm(G);
  (*print_ptr)("\n");
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
    if (VERBOSE == get_verbosity() && j == 100) {
      j = 0;
      if (temp != NULL)
        (*print_ptr)(
            "Iteration number %d; contracting edge starting at (%d,%d)\n", i,
            temp->start, temp->end);
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
    free(v->data);
    free(v);
    return (s_list);
  } else {
    temp = prev->nextState;
    while ((temp != NULL) && (!eq_lift_state(a, temp->data, G))) {
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

void free_lift_state_list(LiftStateList states) {
  LiftStateList temp;
  temp = states;
  while (temp != NULL) {
    states = states-> nextState;
    free(temp->data);
    free(temp);
    temp = states;
  }
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

/* Higher differentials */

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
 * Calculates whether the supplied state is nullhomologous
 * @param init a State
 * @param G working grid
 * @return nonzero if nullhomologous and zero otherwise.
 */
int null_homologous_D0Q(const State init, const Grid_t *const G) {
  // Initialization of variables
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
  // Creates first edge for edge list pointing from A_0 to init
  EdgeList edge_list = prepend_edge(0, 1, NULL);
  // More initialization
  prev_outs = NULL;
  prev_ins = NULL;
  new_ins = malloc(sizeof(StateNode_t));
  new_ins->data = malloc(sizeof(char) * G->arc_index);
  i = 0;
  // sets new_ins to be initial state (IE initialize B_0)
  while (i < G->arc_index) {
    // sets new_ins to be init state
    new_ins->data[i] = init[i];
    i++;
  };
  new_ins->nextState = NULL;
  // This loop goes until either new_ins is empty or we have an answer(ans)
  // The ith loop creates the A_i and B_i sets from A_i-1 and B_i-1 and then
  // contracts edges
  // new_ins is empty implies that there is no more contracting to affect A_0
  // connected edges
  //
  // THIS IS THE ACTUAL ALGORITHM LOOP
  ans = 0;
  int current_pos = 1;
  while (new_ins != NULL && !ans) {
    // sets the present_in states to current working states (B_i-1)
    present_in = new_ins;
    // resets variables from last loop
    in_number = 0;
    num_new_outs = 0;
    new_outs = NULL;
    // loop until there are no in states left to look at in B_i-1
    // this is to make (A_i)
    if (get_verbosity() >= VERBOSE) {
      (*print_ptr)("Gathering A_%d:\n", current_pos);
    }
    while (present_in != NULL) {
      free_state_list(really_new_outs);
      in_number++;
      // sets really_new_outs to be the states with rectangles pointing into
      // present_in (B_i-1)
      // that is not in prev_outs (A_i-1) (this is part of A_i coming from the
      // current
      // working state, present_in)
      really_new_outs = new_rectangles_into(prev_outs, present_in->data, G);
      // loop through all really_new_outs
      while (really_new_outs != NULL) {
        // get position of really_new_outs in new_outs (ie position of current
        // A_i state in the set of all A_i if it is not in it yet it sets
        // to 0)
        out_number = get_number(really_new_outs->data, new_outs, G);
        // if really_new_outs state is not in list this adds it to
        // the statelist new_outs
        if (out_number == 0) {
          // creates new_outs if empty
          if (num_new_outs == 0) {
            new_outs = really_new_outs;
            really_new_outs = really_new_outs->nextState;
            new_outs->nextState = NULL;
            last_new_out = new_outs;
            num_new_outs++;
            out_number = num_new_outs;
          } else {
            // appends to new_outs if non-empty
            last_new_out->nextState = really_new_outs;
            really_new_outs = really_new_outs->nextState;
            last_new_out = last_new_out->nextState;
            last_new_out->nextState = NULL;
            num_new_outs++;
            out_number = num_new_outs;
          };
        } else {
          // removing data and skipping state if its already in the list
          temp = really_new_outs;
          really_new_outs = really_new_outs->nextState;
          free(temp->data);
          free(temp);
        }
        // Appends edge to edge_list and then increments edge count
        edge_list = append_ordered(out_number + num_outs, in_number + num_ins,
                                   edge_list);

        edge_count++;
      }
      // goes to the next state in B_i-1 to look for rectangles into it
      present_in = present_in->nextState;
    };
    if (get_verbosity() >= VERBOSE) {
      print_edges(edge_list);
      (*print_ptr)("\n");
    }
    // Initialize things to calculate B_i from A_i
    free_state_list(prev_ins);
    prev_ins = new_ins;
    i = 1;
    num_ins = num_ins + in_number;
    prev_in_number = num_ins;
    num_new_ins = 0;
    new_ins = NULL;
    out_number = 0;
    present_out = new_outs;
    // loop until there is no states left to look at in A_i to make B_i
    if (get_verbosity() >= VERBOSE) {
      (*print_ptr)("Gathering B_%d:\n", current_pos);
    }
    while (present_out != NULL) {
      out_number++;
      // set really_new_ins to be states with rectangles pointing into them from
      // the current state in present_out (ie really_new_ins is the part of B_i
      // that the current state in A_i has a recangle pointing to)
      really_new_ins = new_rectangles_out_of(prev_ins, present_out->data, G);
      // loops through really_new_ins
      while (really_new_ins != NULL) {
        // Checks position of current really_new_ins state in new_ins (B_i),
        // initializes to 0 if it is not yet in the set
        in_number = get_number(really_new_ins->data, new_ins, G);
        // if really_new_ins state is not in new_ins this checks to add it to
        // new_ins (B_i)
        if (in_number == 0) {
          if (num_new_ins == 0) {
            // creats new_ins if it does not yet exist
            new_ins = really_new_ins;
            really_new_ins = really_new_ins->nextState;
            new_ins->nextState = NULL;
            last_new_in = new_ins;
            num_new_ins++;
            in_number = num_new_ins;
          } else {
            // appends really_new_ins state to new_ins (B_i)
            last_new_in->nextState = really_new_ins;
            really_new_ins = really_new_ins->nextState;
            last_new_in = last_new_in->nextState;
            last_new_in->nextState = NULL;
            num_new_ins++;
            in_number = num_new_ins;
          };
        } else {
          // Frees memory and goes to next state if already in new_ins
          temp = really_new_ins;
          really_new_ins = really_new_ins->nextState;
          free(temp->data);
          free(temp);
        }
        // appends the edge found btw A_i and B_i to edge_list
        edge_list = append_ordered(out_number + num_outs, in_number + num_ins,
                                   edge_list);
        edge_count++;
      };
      // Goes to next state in A_i to check for new states that rectangles in
      // A_i
      // point to
      present_out = present_out->nextState;
    };
    if (get_verbosity() >= VERBOSE) {
      print_edges(edge_list);
      (*print_ptr)("\n");
    }
    // Clears memory of unneccary info (A_i-1) and initializes things for next
    // iteration of the algorithm
    free_state_list(prev_outs);
    prev_outs = new_outs;
    new_outs = NULL;
    // Contract edges in edge list
    if (get_verbosity() >= VERBOSE) {
      (*print_ptr)("Contracting edges from 0 to %d:\n", prev_in_number);
    }
    special_homology(0, prev_in_number, &edge_list);
    if (get_verbosity() >= VERBOSE) {
      print_edges(edge_list);
      (*print_ptr)("\n");
    }
    // calculates answer (ie if nullhomologous or not) if it can be calculated
    if ((edge_list == NULL) || (edge_list->start != 0)) {
      // nullhomologous if edge_list is empty or A_0 no longer points to
      // anything
      ans = 1;
      if (get_verbosity() >= VERBOSE) {
        (*print_ptr)("No edges pointing out of A_0!\n");
      }
      free_state_list(new_ins);
      free_state_list(new_outs);
      new_ins = NULL;
    } else if (edge_list->end <= prev_in_number) {
      // not nullhomologous if edge still pointing from A_0 to B_i-1
      ans = 0;
      if (get_verbosity() >= VERBOSE) {
        (*print_ptr)("There exist edges pointing from A_0 to B_%d! No future "
                     "contractions will remove this edge!\n",
                     current_pos - 1);
      }
      free_state_list(new_ins);
      free_state_list(new_outs);
      new_ins = NULL;
    } else {
      // set num_out for next loop through algorithm
      num_outs = num_outs + out_number;
      if (get_verbosity() >= VERBOSE) {
        (*print_ptr)("Total number of states in B_i up to B_%d (before any "
                     "contraction): %d \n",
                     current_pos - 1, prev_in_number);
        (*print_ptr)("Total number of states in A_i up to A_%d (before any "
                     "contraction): %d \n",
                     current_pos, num_outs);
        (*print_ptr)("Total number of states in B_i up to B_%d (before any "
                     "contraction): %d \n",
                     current_pos, num_ins + in_number);
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
  // Initialization of variables
  StateList new_ins, new_outs, last_new_in, last_new_out, temp;
  StateList prev_ins, prev_outs;
  StateList really_new_outs = NULL, really_new_ins = NULL;
  EdgeList last_edge;
  int in_number, ans, prev_in_number;
  int out_number;
  int i;
  int edge_count = 0;
  int num_ins = 0;
  int num_outs = 0;
  int num_new_ins = 0;
  int num_new_outs = 0;
  StateList present_in, present_out;
  prev_outs = NULL;
  prev_ins = NULL;
  new_ins = fixed_wt_rectangles_out_of(1, init, G);
  EdgeList edge_list = prepend_edge(0, 1, NULL);
  temp = new_ins;
  i = 1;
  last_edge = NULL;
  if (temp != NULL) {
    i = 1;
    edge_list = create_edge(0, 1);
    last_edge = edge_list;
    temp = temp->nextState;
    while (temp != NULL) {
      i++;
      last_edge->nextEdge = create_edge(0, i);
      last_edge = last_edge->nextEdge;
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
        out_number = get_number(really_new_outs->data, new_outs, G);
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
          free(temp->data);
          free(temp);
        }
        edge_list = append_ordered(out_number + num_outs, in_number + num_ins,
                                   edge_list);
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
        in_number = get_number(really_new_ins->data, new_ins, G);
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
          free(temp->data);
          free(temp);
        }
        edge_list = append_ordered(out_number + num_outs, in_number + num_ins,
                                   edge_list);
        edge_count++;
      };
      present_out = present_out->nextState;
    };
    free_state_list(prev_outs);
    prev_outs = new_outs;
    new_outs = NULL;
    special_homology(0, prev_in_number, &edge_list);
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
      if (VERBOSE == get_verbosity()) {
        (*print_ptr)("%d %d %d\n", num_ins, num_outs, edge_count);
      }
    };
  };
  return (ans);
}

int null_homologous_lift(const LiftState init, const LiftGrid_t *const G) {
  // Initialization of variables
  LiftStateList new_ins, new_outs, last_new_in, last_new_out, temp;
  LiftStateList prev_ins, prev_outs;
  LiftStateList really_new_outs = NULL, really_new_ins = NULL;
  int in_number, ans, prev_in_number;
  int out_number;
  int edge_count = 0;
  int num_ins = 0;
  int num_outs = 0;
  int num_new_ins = 0;
  int num_new_outs = 0;
  LiftStateList present_in, present_out;
  EdgeList edge_list = prepend_edge(0, 1, NULL);
  prev_outs = NULL;
  prev_ins = NULL;

  new_ins = malloc(sizeof(LiftStateNode_t));
  init_lift_state(&new_ins->data, G);
  for (int i = 0; i < G->sheets; ++i) {
    for (int j = 0; j < G->arc_index; ++j) {
      new_ins->data[i][j] = init[i][j];
    }
  }
  new_ins->nextState = NULL;

  ans = 0;
  int current_pos = 1;
  int i;
  while (new_ins != NULL && !ans) {
    present_in = new_ins;
    in_number = 0;
    num_new_outs = 0;
    new_outs = NULL;
    if (get_verbosity() >= VERBOSE) {
      (*print_ptr)("Gathering A_%d:\n", current_pos);
    }
    while (present_in != NULL) {
      free_lift_state_list(really_new_outs);
      in_number++;
      really_new_outs = new_lift_rectangles_into(prev_outs, present_in->data, G);
      while (really_new_outs != NULL) {
        out_number = get_lift_number(really_new_outs->data, new_outs, G);
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
          free(temp->data);
          free(temp);
        }
        edge_list = append_ordered(out_number + num_outs, in_number + num_ins,
                                   edge_list);

        edge_count++;
      }
      present_in = present_in->nextState;
    };
    if (get_verbosity() >= VERBOSE) {
      print_edges(edge_list);
      (*print_ptr)("\n");
    }
    free_lift_state_list(prev_ins);
    prev_ins = new_ins;
    i = 1;
    num_ins = num_ins + in_number;
    prev_in_number = num_ins;
    num_new_ins = 0;
    new_ins = NULL;
    out_number = 0;
    present_out = new_outs;
    if (get_verbosity() >= VERBOSE) {
      (*print_ptr)("Gathering B_%d:\n", current_pos);
    }
    while (present_out != NULL) {
      out_number++;
      really_new_ins = new_lift_rectangles_out_of(prev_ins, present_out->data, G);
      while (really_new_ins != NULL) {
        in_number = get_lift_number(really_new_ins->data, new_ins, G);
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
          free(temp->data);
          free(temp);
        }
        edge_list = append_ordered(out_number + num_outs, in_number + num_ins,
                                   edge_list);
        edge_count++;
      };
      present_out = present_out->nextState;
    };
    if (get_verbosity() >= VERBOSE) {
      print_edges(edge_list);
      (*print_ptr)("\n");
    }
    free_lift_state_list(prev_outs);
    prev_outs = new_outs;
    new_outs = NULL;
    if (get_verbosity() >= VERBOSE) {
      (*print_ptr)("Contracting edges from 0 to %d:\n", prev_in_number);
    }
    special_homology(0, prev_in_number, &edge_list);
    if (get_verbosity() >= VERBOSE) {
      print_edges(edge_list);
      (*print_ptr)("\n");
    }
    if ((edge_list == NULL) || (edge_list->start != 0)) {
      ans = 1;
      if (get_verbosity() >= VERBOSE) {
        (*print_ptr)("No edges pointing out of A_0!\n");
      }
      free_lift_state_list(new_ins);
      free_lift_state_list(new_outs);
      new_ins = NULL;
    } else if (edge_list->end <= prev_in_number) {
      ans = 0;
      if (get_verbosity() >= VERBOSE) {
        (*print_ptr)("There exist edges pointing from A_0 to B_%d! No future "
                     "contractions will remove this edge!\n",
                     current_pos - 1);
      }
      free_lift_state_list(new_ins);
      free_lift_state_list(new_outs);
      new_ins = NULL;
    } else {
      num_outs = num_outs + out_number;
      if (get_verbosity() >= VERBOSE) {
        (*print_ptr)("Total number of states in B_i up to B_%d (before any "
                     "contraction): %d \n",
                     current_pos - 1, prev_in_number);
        (*print_ptr)("Total number of states in A_i up to A_%d (before any "
                     "contraction): %d \n",
                     current_pos, num_outs);
        (*print_ptr)("Total number of states in B_i up to B_%d (before any "
                     "contraction): %d \n",
                     current_pos, num_ins + in_number);
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
 * Prints the permutations of the Grid, ie the X O code inputed
 * by the user
 * @param G working grid
 */
void print_grid_perm(const Grid_t *const G) {
  int i = 0;
  (*print_ptr)("X's [");
  while (i < G->arc_index) {
    (*print_ptr)(" %d", G->Xs[i]);
    if (i != G->arc_index - 1)
      (*print_ptr)(",");
    i++;
  }
  (*print_ptr)(" ]\nO's [");
  i = 0;
  while (i < G->arc_index) {
    (*print_ptr)(" %d", G->Os[i]);
    if (i != G->arc_index - 1)
      (*print_ptr)(",");
    i++;
  }
  (*print_ptr)(" ]\n");
}

void print_self_link(const Grid_t *const G) {
  int Writhe = 0;
  int up_down_cusps[2] = {0, 0};
  int tb;
  int r;
  Writhe = writhe(G);
  cusps(up_down_cusps, G);
  tb = Writhe - .5 * (up_down_cusps[0] + up_down_cusps[1]);
  r = .5 * (up_down_cusps[1] - up_down_cusps[0]);
  (*print_ptr)("Writhe = %d\n", Writhe);
  (*print_ptr)("Up Cusps: %d\nDown Cusps: %d\n", up_down_cusps[0],
               up_down_cusps[1]);
  (*print_ptr)("tb(G) = %d\n", tb);
  (*print_ptr)("r(G) = %d\n", r);
  (*print_ptr)("2A(x+) = M(x+) = sl(x+)+1 = %d\n\n", tb - r + 1);    
} 

/* Computes the writhe of the passed grid
 * @param G working grid
 * @return writhe of the grid
 */
int writhe(const Grid_t *const G) {
  int i = 1, j = 0, k = 0;
  int maxXO, minXO;
  int writhe = 0;
  int temp_X, temp_O;
  int current_X, current_O;
  while (i < G->arc_index) {
    temp_X = G->Xs[i] - '0';
    temp_O = G->Os[i] - '0';
    minXO = min(temp_X, temp_O);
    if (minXO == temp_X)
      maxXO = temp_O;
    else
      maxXO = temp_X;
    j = 0;
    while (j < i) {
      current_X = G->Xs[j] - '0';
      current_O = G->Os[j] - '0';
      if (minXO < current_X && maxXO > current_X) {
        k = i + 1;
        while (k < G->arc_index) {
          if (G->Os[k] - '0' == current_X) {
            if (maxXO == temp_X)
              writhe++;
            else
              writhe--;
          }
          k++;
        }
      }
      if (minXO < current_O && maxXO > current_O) {
        k = i + 1;
        while (k < G->arc_index) {
          if (G->Xs[k] - '0' == current_O) {
            if (maxXO == temp_O)
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

/*computes up_down_cusps array, which stores the number of up cusps in first
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
