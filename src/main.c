/**
 * Copyright (C) 2019 Lucas Meyers, Brandon Roberts
 *   
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <argp.h>
#include "TransverseHFK.h"

const char *argp_program_version = "transverseHFK revision 0.0.1";
const char *argp_program_bug_address = "<lmeye22@lsu.edu>";
static const char doc[] =
    "A program to calculate the Legendrian/Transverse knot invariants\
 via the algorithm described in \"Transverse knots distinguished by\
 Knot Floer Homology\" by L. Ng, P. S. Ozsvath, and D. P. Thurston.\
 If the number of sheets is not equal to 1 it instead calculates the\
 theta invariant for the n-fold cyclic cover.";

static const char args_doc[] = "-i <ArcIndex> -n <Sheets:1> -X [<Xs>] -O [<Os>]";

static struct argp_option options[] = {
    {"verbose", 'v', 0, 0, "Produce verbose output", 0},
    {"quiet", 'q', 0, 0, "Produce some extraneous output", 0},
    {"silent", 's', 0, 0, "Don't produce any extraneous output", 0},
    {"index", 'i', "ArcIndex", 0, "ArcIndex of the grid", 0},
    {"Xs", 'X', "[...]", 0, "List of Xs", 0},
    {"Os", 'O', "[...]", 0, "List of Os", 0},
    {"sheets", 'n', "SHEETS", 0, "Number of sheets for cyclic branch cover. Default: 1", 0},
    {"timeout", 't', "SECONDS", 0, "Maximum time to run in seconds", 0},
    {0}};

static error_t parse_opt(int, char *, struct argp_state *);
void timeout(const int);
int build_permutation(State, char *, int);

static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};
struct arguments {
  int arc_index;
  int sheets;
  char *Xs;
  char *Os;
  int max_time;
};

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
      argp_failure(state, 0, 0, "Invalid timeout");
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
      printf("Invalid grid\n");
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
      printf("Calculating graph for lifted invariant.\n");
      // These print statements are wrong
      //print_state(G.Xs, &G);
      //print_self_link(&G);
    }

    // Change the content of these print statements
    if (null_homologous_lift(UR_lift, &G)) {
      printf("theta_%d is null-homologous\n", G.sheets);
    }
    else {
      printf("theta_%d is NOT null-homologous\n", G.sheets);
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
    printf("Invalid grid\n");
    free(G.Xs);
    free(G.Os);
    free(UR);
    exit(1);
  }

  if (args.max_time > 0) {
    if (signal(SIGALRM, timeout) == SIG_ERR) {
      perror("An error occured while setting the timer");
      free(G.Xs);
      free(G.Os);
      free(UR);
      exit(1);
    }
    alarm(args.max_time);
  }
  
  if(QUIET <= get_verbosity()) {
    print_grid(&G);
    print_tb_r(&G);
  }

  if (QUIET <= get_verbosity()) {
    printf("\n \nCalculating graph for LL invariant\n");
    print_state(G.Xs,&G);
    print_2AM(&G,0);
  }
  if (null_homologous_D0Q(G.Xs, &G)) {
    printf("LL is null-homologous\n");
  } else {
    printf("LL is NOT null-homologous\n");
  }

  if (QUIET <= get_verbosity()) {
    printf("\nCalculating graph for UR invariant\n");
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
    print_2AM(&G,1);
  }

  if (null_homologous_D0Q(UR, &G)) {
    printf("UR is null-homologous\n");
  } else {
    printf("UR is NOT null-homologous\n");
  };

  if (QUIET <= get_verbosity()) {
    printf("\nCalculating graph for D1[LL] invariant\n");
    print_state(G.Xs, &G);
    print_2AM(&G,0);
  }

  if (null_homologous_D1Q(G.Xs, &G)) {
    printf("D1[LL] is null-homologous\n");
  } else {
    printf("D1[LL] is NOT null-homologous\n");
  }

  if (QUIET <= get_verbosity()) {
    printf("\nCalculating graph for D1[UR] invariant\n");
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
    print_2AM(&G,1);
  }

  if (null_homologous_D1Q(UR, &G)) {
    printf("D1[UR] is null-homologous\n");
  } else {
    printf("D1[UR] is NOT null-homologous\n");
  };

  free(G.Xs);
  free(G.Os);
  free(UR);

  return 0;
}
