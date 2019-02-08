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
    {"verbose", 'v',               0, 0,          "Produce verbose output", 0},
    {  "quiet", 'q',               0, 0, "Don't produce extraneous output", 0},
    {  "index", 'i',      "ArcIndex", 0,            "ArcIndex of the grid", 0},
    {     "Xs", 'X',         "[...]", 0,                      "List of Xs", 0},
    {     "Os", 'O',         "[...]", 0,                      "List of Os", 0},
    {"timeout", 't',       "SECONDS", 0,  "Maximum time to run in seconds", 0},
    {0}};

// Temporary location
static error_t parse_opt(int, char *, struct argp_state *);
static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};

bool verbose = false;
int ArcIndex = -1;
int max_time = -1;
char Xs[MAX_INDEX] = {};
char Os[MAX_INDEX] = {};

struct vertex {
  int data;
  struct vertex *nextVertex;
};

typedef struct vertex Vertex;
typedef Vertex *VertexList;

struct shortEdgeNode {
  int start;
  int end;
  struct shortEdgeNode *nextPtr;
};

typedef struct shortEdgeNode ShortEdgeNode;
typedef ShortEdgeNode *ShortEdges;
ShortEdges EdgeList;

typedef struct stateNode StateNode;
typedef StateNode *StateList;

struct stateNode {
  char data[MAX_INDEX];
  StateList nextState;
};

typedef char State[MAX_INDEX];

StateList BigIns;
StateList BigOuts;

VertexList ins;
VertexList outs;

ShortEdges AddModTwoLists(VertexList kids, VertexList parents);
VertexList PrependVertex(int a, VertexList v);
ShortEdges PrependEdge(int a, int b, ShortEdges e);
StateList FixedWtRectanglesOutOf(int wt, State incoming);
StateList RectanglesOutOf(State incoming);
StateList RectanglesInto(State incoming);
StateList SwapCols(int x1, int x2, State incoming);
int GetNumber(State a, StateList b);
void FreeStateList(StateList States);
int NullHomologousD0Q(State init);
int NullHomologousD1Q(State init);
void Homology(void);
void SpecialHomology(int init, int final);
void Contract(int a, int b);
VertexList RemoveVertex(int a, VertexList v);
void CreateD0Graph(State init);
void CreateD1Graph(State init);
StateList CreateStateNode(State state);
ShortEdges CreateEdge(int a, int b);
void FreeVertexList(VertexList vertices);
void FreeShortEdges(ShortEdges e);
StateList RemoveState(State a, StateList v);

int NESWpO(char *x);
int NESWOp(char *x);
int NESWpp(char *x);

void PrintState(State state);
void PrintStateShort(State state);
void PrintEdges(void);
void PrintStates(StateList states);
void PrintMathEdges(void);
void PrintMathEdgesA(ShortEdges edges);
void PrintVertices(VertexList vlist);

void timeout(int);
int perm_len(const char*);
int is_grid(const int,const char*,const char*);
int buildPermutation(char*,char*);
int EqState(State a, State b) { return (!strncmp(a, b, ArcIndex)); }

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
    if(max_time <= 0) {
      argp_failure(state, 0, 0, "Invalid timeout");
      exit(1);
    }
    break;
  case 'i':
    ArcIndex = atoi(arg);
    if (ArcIndex > MAX_INDEX || ArcIndex < 2) {
      argp_failure(state, 0, 0, "ArcIndex value out of range");
      exit(1);
    }
    break;
  case 'X':
    if (-1 == buildPermutation(Xs, arg)) {
      argp_failure(state, 0, 0, "Malformatted Xs");
      exit(1);
    }
    break;
  case 'O':
    if (-1 == buildPermutation(Os, arg)) {
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
void timeout(int sig)
{
  if(SIGALRM == sig) {
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
int buildPermutation(char *perm, char *str) {
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
      break;
    } else if (s[0] == ',') {
      ++s;
    }
  }

  return 0;
}

/**
 * Returns the number of characters in the character array before a zero is found
 * @param p a character array
 * @return the number of characters before a zero is found. Returns -1 if not found
 */
int perm_len(const char* p)
{
  for(int i=0;;++i) {
    if(0 == p[i]) {
      return i;
    }
  }

  // Will segfault before reaching this point if no 0 is found
  return -1;
}

/**
 * Determines whether the grid specified by the parameters is valid
 * @param i the grid size
 * @param Xs a permutation specifying the Xs of the grid
 * @param Os a permutation specifygin the Os of the grid
 * @return 1 if the grid is valid, 0 otherwise
 */
int is_grid(const int i, const char* Xs, const char* Os)
{
  if(perm_len(Xs) != i || perm_len(Os) != i || 1 >= i || MAX_INDEX <= i) {
    return 0;
  }

  for(int j=0; j< i; ++j) {
    if(Xs[j] == Os[j]) {
      return 0;
    }
    int found_x = 0;
    int found_o = 0;

    for(int k=0;k<i && (found_x == 0 || found_o == 0); ++k) {
      if(Xs[k] == j+1) {
        found_x = 1;
      }

      if(Os[k] == j+1) {
        found_o = 1;
      }
    }

    if(0 == found_x || 0 == found_o) {
      return 0;
    }
  }
  
  return 1;
}

int main(int argc, char **argv) {
  argp_parse(&argp, argc, argv, 0, 0, 0);

  char UR[ArcIndex];
  int i;

  if(!is_grid(ArcIndex,Xs,Os)) {
    printf("Invalid input\n");
    exit(1);
  }

  if(max_time > 0) {
    if(signal(SIGALRM, timeout) == SIG_ERR) {
      perror("An error occured while setting the timer");
      exit(1);
    }
    alarm(max_time);
  }

  if(verbose) {
    printf("\n \nCalculating graph for LL invariant\n");
    PrintState(Xs);
  }
  if (NullHomologousD0Q(Xs)) {
    printf("LL is null-homologous\n");
  } else {
    printf("LL is NOT null-homologous\n");
  }

  if (verbose) {
    printf("\n \nCalculating graph for UR invariant\n");
  }
  if (Xs[ArcIndex - 1] == ArcIndex) {
    UR[0] = 1;
  } else {
    UR[0] = (char)Xs[ArcIndex - 1] + 1;
  };
  i = 1;
  while (i <= ArcIndex - 1) {
    if (Xs[i - 1] == ArcIndex) {
      UR[i] = 1;
    } else {
      UR[i] = Xs[i - 1] + 1;
    }
    i++;
  }

  if (verbose) {
    PrintState(UR);
  }

  if (NullHomologousD0Q(UR)) {
    printf("UR is null-homologous\n");
  } else {
    printf("UR is NOT null-homologous\n");
  };

  if (verbose) {
    printf("\n \nCalculating graph for D1[LL] invariant\n");
    PrintState(Xs);
  }

  if (NullHomologousD1Q(Xs)) {
    printf("D1[LL] is null-homologous\n");
  } else {
    printf("D1[LL] is NOT null-homologous\n");
  }

  if (verbose) {
    printf("\n \nCalculating graph for D1[UR] invariant\n");
  }

  if (Xs[ArcIndex - 1] == ArcIndex) {
    UR[0] = 1;
  } else {
    UR[0] = (char)Xs[ArcIndex - 1] + 1;
  };
  i = 1;
  while (i <= ArcIndex - 1) {
    if (Xs[i - 1] == ArcIndex) {
      UR[i] = 1;
    } else {
      UR[i] = Xs[i - 1] + 1;
    }
    i++;
  }

  if (verbose) {
    PrintState(UR);
  }

  if (NullHomologousD1Q(UR)) {
    printf("D1[UR] is null-homologous\n");
  } else {
    printf("D1[UR] is NOT null-homologous\n");
  };

  return 0;
}

/**
 * Shifts the input towards the interval [0,ArcIndex) by
 * a multiple of ArcIndex
 * @param a An integer
 * @return a shifted towards the interval [0,ArcIndex)
 * @see ArcIndex
 */
int Mod(int a) {
  if (a >= ArcIndex) {
    return (a - ArcIndex);
  } else if (a < 0) {
    return (a + ArcIndex);
  } else {
    return (a);
  };
}

/**
 * Shifts the input towards the interval (0,ArcIndex] by
 * a multiple of ArcIndex
 * @param a An integer
 * @return a shifted towards the interval (0,ArcIndex]
 * @see ArcIndex
 */
int ModUp(int a) {
  if (a > ArcIndex) {
    return (a - ArcIndex);
  } else if (a <= 0) {
    return (a + ArcIndex);
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
 * that is not contained in Prevs.
 * @param Prevs Statelist containing previous states
 * @param incoming the source of rectangles used to generate statelist
 * @return A statelist containing states reached from a rectangle leaving
 * incoming not contained in Prevs.
 * @see ArcIndex
 */
StateList NewRectanglesOutOf(StateList Prevs, State incoming) {
  StateList Temp, ans;
  State TempState;
  int LL;
  int w, h, m, n, i;
  ans = NULL;
  i = 0;
  while (i < ArcIndex) {
    TempState[i] = incoming[i];
    i++;
  }
  LL = 0;
  while (LL < ArcIndex) {
    w = 1;
    h = min(Mod(Os[LL] - incoming[LL]), Mod(Xs[LL] - incoming[LL]));
    while (w < ArcIndex && h > 0) {
      if (Mod(incoming[Mod(LL + w)] - incoming[LL]) <= h) {
        TempState[LL] = incoming[Mod(LL + w)];
        TempState[Mod(LL + w)] = incoming[LL];
        m = GetNumber(TempState, Prevs);
        if (m == 0) {
          n = GetNumber(TempState, ans);
          if (n != 0) {
            ans = RemoveState(TempState, ans);
          } else {
            Temp = SwapCols(LL, Mod(LL + w), incoming);
            Temp->nextState = ans;
            ans = Temp;
          }
        };
        TempState[LL] = incoming[LL];
        TempState[Mod(LL + w)] = incoming[Mod(LL + w)];
        h = Mod(incoming[Mod(LL + w)] - incoming[LL]);
      };
      h = min(h, min(Mod(Os[Mod(LL + w)] - incoming[LL]),
                     Mod(Xs[Mod(LL + w)] - incoming[LL])));
      w++;
    };
    LL++;
  };
  return ans;
}

/**
 * Returns a StateList consisting of states that are reached by rectangles
 * leaving incoming
 * @param incoming Initial state for generated rectangles
 * @return list of states reached by a rectangle from incoming.
 * @see ArcIndex
 */
StateList RectanglesOutOf(State incoming) {
  StateList Temp, ans;
  int LL;
  int w, h;
  ans = NULL;
  LL = 0;
  while (LL < ArcIndex) {
    w = 1;
    h = min(Mod(Os[LL] - incoming[LL]), Mod(Xs[LL] - incoming[LL]));
    while (w < ArcIndex && h > 0) {
      if (Mod(incoming[Mod(LL + w)] - incoming[LL]) <= h) {
        Temp = SwapCols(LL, Mod(LL + w), incoming);
        Temp->nextState = ans;
        ans = Temp;
        h = Mod(incoming[Mod(LL + w)] - incoming[LL]);
      };
      h = min(h, min(Mod(Os[Mod(LL + w)] - incoming[LL]),
                     Mod(Xs[Mod(LL + w)] - incoming[LL])));
      w++;
    };
    LL++;
  };
  return ans;
}

/**
 * returns a StateList containing those with a rectangle
 * pointing to the state incoming
 * @param incoming State that is the destination for generated rectangles
 * @return StateList containing states with a rectangle to incoming.
 * @see ArcIndex
 */
StateList RectanglesInto(State incoming) {
  StateList Temp, ans;
  int LL;
  int w, h;
  ans = NULL;
  LL = 0;
  while (LL < ArcIndex) {
    w = 1;
    h = min(ModUp(incoming[LL] - Os[LL]), ModUp(incoming[LL] - Xs[LL]));
    while (w < ArcIndex && h > 0) {
      if (ModUp(incoming[LL] - incoming[Mod(LL + w)]) < h) {
        Temp = SwapCols(LL, Mod(LL + w), incoming);
        Temp->nextState = ans;
        ans = Temp;
        h = ModUp(incoming[LL] - incoming[Mod(LL + w)]);
      };
      h = min(h, min(ModUp(incoming[LL] - Os[Mod(LL + w)]),
                     ModUp(incoming[LL] - Xs[Mod(LL + w)])));
      w++;
    };
    LL++;
  };
  return ans;
}

/**
 * returns a StateList containing those with a rectangle
 * pointing to the state incoming that do not overlap with Prevs
 * @param incoming State that is the destination for generated rectangles
 * @param Prevs StateList of excluded states
 * @return StateList containing states with a rectangle to incoming.
 * @see ArcIndex
 */
StateList NewRectanglesInto(StateList Prevs, State incoming) {
  StateList ans, Temp;
  State TempState;
  int LL, m, n;
  int w, h;
  int i;
  ans = NULL;
  i = 0;
  while (i < ArcIndex) {
    TempState[i] = incoming[i];
    i++;
  }
  LL = 0;
  while (LL < ArcIndex) {
    w = 1;
    h = min(ModUp(incoming[LL] - Os[LL]), ModUp(incoming[LL] - Xs[LL]));
    while (w < ArcIndex && h > 0) {
      if (ModUp(incoming[LL] - incoming[Mod(LL + w)]) < h) {
        TempState[LL] = incoming[Mod(LL + w)];
        TempState[Mod(LL + w)] = incoming[LL];
        m = GetNumber(TempState, Prevs);
        if (m == 0) {
          n = GetNumber(TempState, ans);
          if (n != 0) {
            ans = RemoveState(TempState, ans);
          } else {
            Temp = SwapCols(LL, Mod(LL + w), incoming);
            Temp->nextState = ans;
            ans = Temp;
          }
        };
        TempState[LL] = incoming[LL];
        TempState[Mod(LL + w)] = incoming[Mod(LL + w)];
        h = ModUp(incoming[LL] - incoming[Mod(LL + w)]);
      };
      h = min(h, min(ModUp(incoming[LL] - Os[Mod(LL + w)]),
                     ModUp(incoming[LL] - Xs[Mod(LL + w)])));
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
 * @return a single element StateList containing a copy of
 * the permutation incoming with entries x1 and x2 swapped
 * @see ArcIndex
 */
StateList SwapCols(int x1, int x2, char *incoming) {
  StateList ans;
  int i;
  i = 0;
  ans = malloc(sizeof(StateNode));
  i = 0;
  while (i < ArcIndex) {
    ans->data[i] = incoming[i];
    i++;
  };
  ans->data[x1] = (incoming)[x2];
  ans->data[x2] = (incoming)[x1];
  ans->nextState = NULL;
  return ans;
}

/**
 * Returns the length of the StateList
 * @param states a StateList
 * @return an int equal to the length of states
 */
int LengthStateList(StateList states) {
  int c;
  StateList Temp;
  Temp = states;
  c = 0;
  while (Temp != NULL) {
    c++;
    Temp = Temp->nextState;
  };
  return c;
}

/**
 * Returns the length of a VertexList
 * @param a VertexList
 * @return an int equal to the length of states
 */
int LengthVertexList(VertexList states) {
  int c;
  VertexList Temp;
  Temp = states;
  c = 0;
  while (Temp != NULL) {
    c++;
    Temp = Temp->nextVertex;
  };
  return c;
}

/**
 * Prints states in the form "{<state>,...}" up
 * to the first 500,000 states.
 * @param states a StateList
 * @see PrintStateShort
 */
void PrintStates(StateList states) {
  StateList Temp;
  int c;
  Temp = states;
  printf("{");
  c = 0;
  while ((Temp != NULL) && c < 500000) {
    PrintStateShort(Temp->data);
    Temp = Temp->nextState;
    if (Temp != NULL) {
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
 * @see ArcIndex
 */
void PrintStateShort(State state) {
  int i;
  i = 0;
  printf("{");
  while (i < ArcIndex - 1) {
    printf("%d,", state[i]);
    i++;
  };
  printf("%d}", state[ArcIndex - 1]);
}

/**
 * Prints the permutation of state on the grid specified by
 * Xs and Os, as well as 2A=M=SL+1.
 * @param state A state
 * @see Xs
 * @see Os
 */
void PrintState(State state) {
  int i, j;
  j = ArcIndex;
  while (j > 0) {
    i = 0;
    while (i < ArcIndex) {
      if (Xs[i] == j) {
        printf("  X  ");
      } else {
        if (Os[i] == j) {
          printf("  O  ");
        } else {
          printf("  -  ");
        };
      };
      i++;
    };
    printf("\n");
    i = 0;
    while (i < ArcIndex) {
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
         NESWpp(Xs) - NESWpO(Xs) - NESWOp(Xs) + NESWpp(Os) + 1);
}

/**
 * Returns the index of a State within a StateList
 * Note: Indexed from 1
 * @param a State
 * @param b a StateList
 * @return index of a within b
 */
int GetNumber(State a, StateList b) {
  StateList temp;
  int count = 1;
  temp = b;
  while (temp != NULL) {
    if (EqState(a, temp->data)) {
      return count;
    };
    temp = temp->nextState;
    count++;
  };
  return 0;
}

/**
 * Return a single element StateList with data state
 * @param state a state to initialize the list
 * @return a single element StateList containing the data from state.
 * @see ArcIndex
 */
StateList CreateStateNode(State state) {
  StateList ans;
  int i;
  ans = malloc(sizeof(StateNode));
  ans->nextState = NULL;
  i = 0;
  while (i < ArcIndex) {
    ans->data[i] = state[i];
    i++;
  };
  return (ans);
}

/**
 * Places the given state at the end of the provided StateList
 * @param state a State
 * @param rest a StateList
 * @return appends a stateNode containing state to the end of rest
 * @see ArcIndex
 */
StateList AppendToStateList(State state, StateList rest) {
  StateList NewNode, TTTemp;
  int i;
  NewNode = malloc(sizeof(StateNode));
  i = 0;
  while (i < ArcIndex) {
    NewNode->data[i] = state[i];
    i++;
  };
  NewNode->nextState = NULL;
  if (rest == NULL) {
    return NewNode;
  } else {
    TTTemp = rest;
    while (TTTemp->nextState != NULL) {
      TTTemp = TTTemp->nextState;
    };
    TTTemp->nextState = NewNode;
  };
  return rest;
}

/**
 * Takes in two lists of vertices and adds them mod two. Uses EdgeList
 * @param parents a list of parent vertices
 * @param kids a list of child vertices
 * @return A shortEdges containing the result of adding them mod two.
 * @see EdgeList
 */
ShortEdges AddModTwoLists(VertexList parents, VertexList kids) {
  VertexList thiskid, thisparent, tempvert;
  ShortEdges thisedge, tempedge, Prev;
  ShortEdges ans;
  if ((parents == NULL) || (kids == NULL)) {
    ans = EdgeList; // change to return
  } else {
    thisparent = parents;
    thiskid = kids;
    ans = NULL;
    thisedge = EdgeList;
    while (thisparent != NULL && thisedge != NULL &&
           (thisedge->start == thisparent->data &&
            thisedge->end == thiskid->data)) {
      tempedge = thisedge; // This may need to be freed
      thisedge = thisedge->nextPtr;
      thiskid = thiskid->nextVertex;
      if (thiskid == NULL) {
        tempvert = thisparent;
        thisparent = thisparent->nextVertex;
        free(tempvert);
        thiskid = kids;
      };
    };
    if (thisedge != NULL &&
        (thisparent == NULL || (thisedge->start < thisparent->data ||
                                (thisedge->start == thisparent->data &&
                                 thisedge->end < thiskid->data)))) {
      ans = thisedge;
      Prev = ans;
      thisedge = thisedge->nextPtr;
    } else if (thisparent != NULL) {
      ans = CreateEdge(thisparent->data, thiskid->data);
      ans->nextPtr = NULL;
      thiskid = thiskid->nextVertex;
      if (thiskid == NULL) {
        tempvert = thisparent;
        thisparent = thisparent->nextVertex;
        free(tempvert);
        thiskid = kids;
      };
    } else {
      ans = NULL;
    }
    Prev = ans;
    while (thisedge != NULL && thisparent != NULL) {
      while (thisedge != NULL && ((thisedge->start < thisparent->data ||
                                   (thisedge->start == thisparent->data &&
                                    thisedge->end < thiskid->data)))) {
        Prev->nextPtr = thisedge;
        Prev = Prev->nextPtr;
        thisedge = thisedge->nextPtr;
      };
      while (thisedge != NULL && thisparent != NULL &&
             (thisedge->start == thisparent->data &&
              thisedge->end == thiskid->data)) {
        tempedge = thisedge;
        thisedge = tempedge->nextPtr;
        Prev->nextPtr = thisedge;
        free(tempedge);
        thiskid = thiskid->nextVertex;
        if (thiskid == NULL) {
          tempvert = thisparent;
          thisparent = thisparent->nextVertex;
          free(tempvert);
          thiskid = kids;
        };
      };
      while (thisparent != NULL &&
             ((thisedge == NULL) || ((thisedge->start > thisparent->data ||
                                      (thisedge->start == thisparent->data &&
                                       thisedge->end > thiskid->data))))) {
        Prev->nextPtr = CreateEdge(thisparent->data, thiskid->data);
        Prev = Prev->nextPtr;
        Prev->nextPtr = NULL;
        thiskid = thiskid->nextVertex;
        if (thiskid == NULL) {
          tempvert = thisparent;
          thisparent = thisparent->nextVertex;
          free(tempvert);
          thiskid = kids;
        };
      };
    };
    if ((thisparent == NULL) && (Prev != NULL)) {
      Prev->nextPtr = thisedge;
    }
  };
  return (ans);
}

/**
 * Places a new edge into the supplied edgelist in order.
 * @param a an int specifying the parent vertex
 * @param b an int specifying the child vertex
 * @param edges a ShortEdges where the new edge will be placed
 * @return a pointer to edges with the new edge added
 */
ShortEdges AppendOrdered(int a, int b, ShortEdges edges) {
  ShortEdges Temp, Prev, curr, ans;
  Prev = edges;
  if ((edges == NULL) || (edges->start > a) ||
      (edges->start == a && edges->end > b)) {
    ans = malloc(sizeof(ShortEdgeNode));
    ans->start = a;
    ans->end = b;
    ans->nextPtr = Prev;
  } else {
    ans = edges;
    curr = Prev->nextPtr;
    while (curr != NULL &&
           ((curr->start < a) || ((curr->start == a) && (curr->end < b)))) {
      curr = curr->nextPtr;
      Prev = Prev->nextPtr;
    };
    Temp = malloc(sizeof(ShortEdgeNode));
    Temp->start = a;
    Temp->end = b;
    Temp->nextPtr = curr;
    Prev->nextPtr = Temp;
  };
  return (ans);
}

/**
 * Add a new edge to the start of a ShortEdges
 * @param a an int indictating the source of the edge
 * @param b an int indictating the destination of the edge
 * @param e a ShortEdges
 * @return e with the edge (a,b) at the front
 */
ShortEdges PrependEdge(int a, int b, ShortEdges e) {
  ShortEdges newPtr;
  newPtr = malloc(sizeof(ShortEdgeNode));
  newPtr->start = a;
  newPtr->end = b;
  newPtr->nextPtr = e;
  return (newPtr);
}

/**
 * Creates a single element ShortEdge
 * @param a an int indicating the source of the edge
 * @param b an int indicating the destination of the edge
 * @return a single element ShortEdges (a,b)
 */
ShortEdges CreateEdge(int a, int b) {
  ShortEdges newPtr;
  newPtr = malloc(sizeof(ShortEdgeNode));
  newPtr->start = a;
  newPtr->end = b;
  newPtr->nextPtr = NULL;
  return (newPtr);
}

/**
 * Creates a vertex and adds it to the start of the supplied
 * VertexList
 * @param a an int specifying the vertex
 * @param vertices a VertexList
 * @return a VertexList with a Vertex containing a prepended to vertices
 */
VertexList PrependVertex(int a, VertexList vertices) {
  VertexList newPtr;
  newPtr = malloc(sizeof(Vertex));
  (newPtr->data) = a;
  (newPtr->nextVertex) = vertices;
  return newPtr;
}

/**
 * Calculate the Homology of Edgelist and storing the result
 * in EdgeList.
 * @see EdgeList
 */
void Homology() {
  ShortEdges Temp;
  Temp = EdgeList;
  while (Temp != NULL) {
    if (Temp != NULL) {
      Contract(Temp->start, Temp->end);
      Temp = EdgeList;
    };
  };
}

/**
 * Calculates the homology of EdgeList where EdgeList must
 * have a specified initial state and a terminates at a specified
 * final state.
 * @param init an int specifying the required start
 * @param final State to t
 * @return
 * @see EdgeList
 */
void SpecialHomology(int init, int final) {
  int i, j, t;
  ShortEdges Temp;
  i = 0;
  j = 0;
  Temp = EdgeList;
  if ((EdgeList == NULL) || (EdgeList->start != init)) {
    printf("FOOO");
    scanf("%d", &t);
    FreeShortEdges(EdgeList);
    EdgeList = NULL;
  };
  while ((EdgeList != NULL) && (Temp != NULL)) {
    while ((Temp != NULL) && (Temp->start == init)) {
      Temp = Temp->nextPtr;
    };
    while ((Temp != NULL) && (Temp->end > final)) {
      Temp = Temp->nextPtr;
    };
    if (verbose && j == 100) {
      j = 0;
      if (Temp != NULL)
        printf("Iteration number %d; Contracting edge starting at (%d,%d)\n", i,
               Temp->start, Temp->end);
    };
    i++;
    j++;
    if (Temp != NULL) {
      Contract(Temp->start, Temp->end);
      Temp = EdgeList;
    };
  };
}

/**
 * Contracts the edge specified by the input within EdgeList
 * @param a the parent vertex of the edge
 * @param b the child vertex of the edge
 * @see EdgeList
 */
void Contract(int a, int b) {
  ShortEdges Temp;
  ShortEdges Prev;
  VertexList parents, kids, tempkids, tempparents;
  VertexList LastParent, LastKid;
  Prev = EdgeList; // Multiple equal initializations
  parents = NULL;
  kids = NULL;
  tempkids = NULL;
  tempparents = NULL;
  LastParent = NULL; 
  LastKid = NULL;
  while (EdgeList != NULL && ((EdgeList)->end == b || EdgeList->start == a)) {
    if ((EdgeList->end == b) && (EdgeList->start == a)) {
      Temp = EdgeList;
      EdgeList = EdgeList->nextPtr;
      free(Temp);
    } else {
      if (EdgeList->end == b) {
        if (LastParent == NULL) {
          parents = PrependVertex(EdgeList->start, NULL);
          LastParent = parents;
        } else {
          tempparents = PrependVertex(EdgeList->start, NULL);
          LastParent->nextVertex = tempparents;
          LastParent = LastParent->nextVertex;
        };
        Temp = EdgeList;
        EdgeList = EdgeList->nextPtr;
        free(Temp);
      } else if (EdgeList->start == a) {
        if (LastKid == NULL) {
          kids = PrependVertex(EdgeList->end, kids);
          LastKid = kids;
        } else {
          tempkids = PrependVertex(EdgeList->end, NULL);
          LastKid->nextVertex = tempkids;
          LastKid = LastKid->nextVertex;
        };
        Temp = EdgeList;
        EdgeList = EdgeList->nextPtr;
        free(Temp);
      };
    };
  };
  Prev = EdgeList;
  if (EdgeList != NULL) {
    Temp = (EdgeList->nextPtr);
  } else {
    Temp = NULL;
  };
  while (Temp != NULL && (Temp)->start < a) {
    if ((Temp)->end == b) {
      if (LastParent == NULL) {
        parents = PrependVertex((Temp)->start, NULL);
        LastParent = parents;
      } else {
        tempparents = PrependVertex((Temp)->start, NULL);
        LastParent->nextVertex = tempparents;
        LastParent = LastParent->nextVertex;
      };
      (Prev->nextPtr) = (Temp->nextPtr);
      free(Temp);
      Temp = Prev->nextPtr;
    } else {
      Temp = (Temp)->nextPtr;
      Prev = (Prev)->nextPtr;
    };
  };
  while (Temp != NULL && (Temp)->start == a) {
    if (Temp->end != b) {
      if (LastKid == NULL) {
        kids = PrependVertex(Temp->end, NULL);
        LastKid = kids;
      } else {
        tempkids = PrependVertex(Temp->end, NULL);
        LastKid->nextVertex = tempkids;
        LastKid = LastKid->nextVertex;
      };
    };
    (Prev)->nextPtr = Temp->nextPtr;
    free(Temp);
    Temp = (Prev)->nextPtr;
  };
  while (Temp != NULL) {
    if ((Temp)->end == b) {
      if (LastParent == NULL) {
        parents = PrependVertex(Temp->start, NULL);
        LastParent = parents;
      } else {
        tempparents = PrependVertex((Temp)->start, NULL);
        LastParent->nextVertex = tempparents;
        LastParent = LastParent->nextVertex;
      };
      (Prev)->nextPtr = (Temp)->nextPtr;
      free(Temp);
      Temp = (Prev)->nextPtr;
    } else {
      Temp = (Temp)->nextPtr;
      Prev = (Prev)->nextPtr;
    }
  };
  EdgeList = AddModTwoLists(parents, kids);
}

/**
 * Removes a state from a StateList. Equality checked using EqState.
 * @param a a State
 * @param v a StateList
 * @return a StateList removing a from v
 * @see EqState
 */
StateList RemoveState(State a, StateList v) {
  StateList Temp, Prev;
  StateList sList = v;
  Prev = v;
  if (v == NULL)
    return (NULL);
  else if (EqState(a, v->data)) {
    Temp = v;
    sList = v->nextState;
    free(v);
    return (sList);
  } else {
    Temp = Prev->nextState;
    while ((Temp != NULL) && (!EqState(a, Temp->data))) {
      Temp = Temp->nextState;
      Prev = Prev->nextState;
    };
    if (Temp != NULL) {
      Prev->nextState = Temp->nextState;
      free(Temp);
      return sList;
    } else
      return sList;
  }
}

/**
 * Removes a vertex containing the supplied int from a VertexList
 * @param a an int
 * @param v a VertexList
 * @return a VertexList with the Vertex containing a removed
 */
VertexList RemoveVertex(int a, VertexList v) {
  VertexList Temp, Prev;
  VertexList vList = v;
  Prev = v;
  if (v == NULL)
    return NULL;
  else if (v->data == a) {
    Temp = v;
    vList = v->nextVertex;
    free(Temp);
    return (vList);
  } else {
    Temp = Prev->nextVertex;
    while ((Temp != NULL) && (Temp->data) < a) {
      Temp = Temp->nextVertex;
      Prev = Prev->nextVertex;
    };
    if (Temp != NULL) {
      Prev->nextVertex = Temp->nextVertex;
      free(Temp);
      return vList;
    } else
      return vList;
  };
}

/**
 * Prints each edge in EdgeList on a new line
 * @see EdgeList
 */
void PrintEdges() {
  ShortEdges Temp;
  Temp = EdgeList;
  while (Temp != NULL) {
    printf("%d %d\n", Temp->start, Temp->end);
    Temp = (Temp->nextPtr);
  };
}

/**
 * Print the first 80 edges in EdgeList on the same line
 * @see EdgeList
 */
void PrintMathEdges() {
  ShortEdges Temp;
  int t;
  Temp = EdgeList;
  printf("{");
  t = 0;
  while (Temp != NULL) {
    printf("{%d,%d}", Temp->start, Temp->end);
    t++;
    if (t == 80) {
      Temp = NULL;
      printf("...}\n");
    } else {
      Temp = (Temp->nextPtr);
      if (Temp != NULL)
        printf(",");
    }
  };
  printf("}");
}

/**
 * Prints the edges in edges on a single line
 * @param edges
 */
void PrintMathEdgesA(ShortEdges edges) {
  ShortEdges Temp;
  Temp = edges;
  printf("{");
  while (Temp != NULL) {
    printf("{%d,%d}", Temp->start, Temp->end);
    Temp = (Temp->nextPtr);
    if (Temp != NULL)
      printf(",");
  };
  printf("}");
}

/**
 * Prints the vertices in VertexList on a single line
 * @param vlist
 */
void PrintVertices(VertexList vlist) {
  VertexList temp;
  temp = vlist;
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
 * @param states
 */
void FreeStateList(StateList states) {
  StateList Temp;
  Temp = states;
  while (Temp != NULL) {
    states = states->nextState;
    free(Temp);
    Temp = states;
  };
}

/**
 * Frees the supplied ShortEdges
 * @param e
 */
void FreeShortEdges(ShortEdges e) {
  ShortEdges Temp, nTemp;
  Temp = e;
  nTemp = Temp;
  while (Temp != NULL) {
    nTemp = Temp;
    Temp = Temp->nextPtr;
    free(nTemp);
  };
}

/**
 * Frees the supplied VertexList
 * @param vertices
 */
void FreeVertexList(VertexList vertices) {
  VertexList Temp, nTemp;
  Temp = vertices;
  nTemp = Temp;
  while (Temp != NULL) {
    nTemp = Temp;
    Temp = Temp->nextVertex;
    free(nTemp);
  };
}

/* Higher differentials */

/**
 * Calculates a StateList containing all states reachable
 * by a rectangles of a fixed width
 * @param wt an int specifying rectangle width
 * @param incoming origin state for the rectangles
 * @return a StateList with states that are reached by a rectangle of width
 * wt from incoming.
 * @see ArcIndex
 * @see Xs
 * @see Os
 */
StateList FixedWtRectanglesOutOf(int wt, State incoming) {
  StateList Temp, ans;
  int LL;
  int w, h;
  int thisweight, i;
  ans = NULL;
  LL = 0;
  while (LL < ArcIndex) {
    w = 1;
    h = Mod(Os[LL] - incoming[LL]);
    while (w < ArcIndex && h > 0) {
      if (Mod(incoming[Mod(LL + w)] - incoming[LL]) <= h) {
        thisweight = 0;
        i = 0;
        while (i < w && thisweight <= wt + 1) {
          if (Mod(Xs[Mod(LL + i)] - incoming[LL]) <
              Mod(incoming[Mod(LL + w)] - incoming[LL])) {
            thisweight++;
          };
          i++;
        };
        if (thisweight == wt) {
          Temp = SwapCols(LL, Mod(LL + w), incoming);
          if (GetNumber(Temp->data, ans) != 0) {
            ans = RemoveState(Temp->data, ans);
            free(Temp);
            Temp = ans;
          } else {
            Temp->nextState = ans;
            ans = Temp;
          }
        };
        h = Mod(incoming[Mod(LL + w)] - incoming[LL]);
      };
      h = min(h, Mod(Os[Mod(LL + w)] - incoming[LL]));
      w++;
    };
    LL++;
  };
  return ans;
}

/**
 * Calculates whether the supplied state is nullhomologous. Uses
 the global variable EdgeList.
 * @param init a State
 * @return nonzero if nullhomologous and zero otherwise.
 * @see EdgeList
 * @see ArcIndex
 */
int NullHomologousD0Q(State init) {
  StateList NewIns, NewOuts, LastNewIn, LastNewOut, Temp;
  StateList PrevIns, PrevOuts;
  StateList ReallyNewOuts = NULL, ReallyNewIns = NULL;
  int innumber, ans, previnnumber;
  int outnumber;
  int i;
  int edgecount = 0;
  int numIns = 0;
  int numOuts = 0;
  int numNewIns = 0;
  int numNewOuts = 0;
  StateList PresentIn, PresentOut;
  EdgeList = PrependEdge(0, 1, NULL);
  PrevOuts = NULL;
  PrevIns = NULL;
  NewIns = malloc(sizeof(StateNode));
  i = 0;
  while (i < ArcIndex) {
    NewIns->data[i] = init[i];
    i++;
  };
  NewIns->nextState = NULL;
  ans = 0;
  while (NewIns != NULL && !ans) {
    PresentIn = NewIns;
    innumber = 0;
    numNewOuts = 0;
    NewOuts = NULL;
    while (PresentIn != NULL) {
      FreeStateList(ReallyNewOuts);
      innumber++;
      ReallyNewOuts = NewRectanglesInto(PrevOuts, PresentIn->data);
      while (ReallyNewOuts != NULL) {
        outnumber = GetNumber(ReallyNewOuts->data, NewOuts);
        if (outnumber == 0) {
          if (numNewOuts == 0) {
            NewOuts = ReallyNewOuts;
            ReallyNewOuts = ReallyNewOuts->nextState;
            NewOuts->nextState = NULL;
            LastNewOut = NewOuts;
            numNewOuts++;
            outnumber = numNewOuts;
          } else {
            LastNewOut->nextState = ReallyNewOuts;
            ReallyNewOuts = ReallyNewOuts->nextState;
            LastNewOut = LastNewOut->nextState;
            LastNewOut->nextState = NULL;
            numNewOuts++;
            outnumber = numNewOuts;
          };
        } else {
          Temp = ReallyNewOuts;
          ReallyNewOuts = ReallyNewOuts->nextState;
          free(Temp);
        }
        EdgeList =
            AppendOrdered(outnumber + numOuts, innumber + numIns, EdgeList);
        edgecount++;
      }
      PresentIn = PresentIn->nextState;
    };
    FreeStateList(PrevIns);
    PrevIns = NewIns;
    i = 1;
    numIns = numIns + innumber;
    previnnumber = numIns;
    numNewIns = 0;
    NewIns = NULL;
    outnumber = 0;
    PresentOut = NewOuts;
    while (PresentOut != NULL) {
      outnumber++;
      ReallyNewIns = NewRectanglesOutOf(PrevIns, PresentOut->data);
      while (ReallyNewIns != NULL) {
        innumber = GetNumber(ReallyNewIns->data, NewIns);
        if (innumber == 0) {
          if (numNewIns == 0) {
            NewIns = ReallyNewIns;
            ReallyNewIns = ReallyNewIns->nextState;
            NewIns->nextState = NULL;
            LastNewIn = NewIns;
            numNewIns++;
            innumber = numNewIns;
          } else {
            LastNewIn->nextState = ReallyNewIns;
            ReallyNewIns = ReallyNewIns->nextState;
            LastNewIn = LastNewIn->nextState;
            LastNewIn->nextState = NULL;
            numNewIns++;
            innumber = numNewIns;
          };
        } else {
          Temp = ReallyNewIns;
          ReallyNewIns = ReallyNewIns->nextState;
          free(Temp);
        }
        EdgeList =
            AppendOrdered(outnumber + numOuts, innumber + numIns, EdgeList);
        edgecount++;
      };
      PresentOut = PresentOut->nextState;
    };
    FreeStateList(PrevOuts);
    PrevOuts = NewOuts;
    NewOuts = NULL;
    SpecialHomology(0, previnnumber);
    if ((EdgeList == NULL) || (EdgeList->start != 0)) {
      ans = 1;
      FreeStateList(NewIns);
      FreeStateList(NewOuts);
      NewIns = NULL;
    } else if (EdgeList->end <= previnnumber) {
      ans = 0;
      FreeStateList(NewIns);
      FreeStateList(NewOuts);
      NewIns = NULL;
    } else {
      numOuts = numOuts + outnumber;
      if (verbose) {
        printf("%d %d %d\n", numIns, numOuts, edgecount);
      }
    };
  };
  return (ans);
}

/**
 * Calculates if D1 of the supplied state is nullhomologous
 * @param init a State
 * @return nonzero if nullhomologous and zero otherwise
 * @see EdgeList
 * @see ArcIndex
 */
int NullHomologousD1Q(State init) {
  StateList NewIns, NewOuts, LastNewIn, LastNewOut, Temp;
  StateList PrevIns, PrevOuts;
  StateList ReallyNewOuts = NULL, ReallyNewIns = NULL;
  ShortEdges LastEdge;
  int innumber, ans, previnnumber;
  int outnumber;
  int i;
  int edgecount = 0;
  int numIns = 0;
  int numOuts = 0;
  int numNewIns = 0;
  int numNewOuts = 0;
  StateList PresentIn, PresentOut;
  LastEdge = EdgeList;
  PrevOuts = NULL;
  PrevIns = NULL;
  NewIns = FixedWtRectanglesOutOf(1, init);
  EdgeList = PrependEdge(0, 1, NULL);
  Temp = NewIns;
  i = 1;
  LastEdge = NULL;
  if (Temp != NULL) {
    i = 1;
    EdgeList = CreateEdge(0, 1);
    LastEdge = EdgeList;
    Temp = Temp->nextState;
    while (Temp != NULL) {
      i++;
      LastEdge->nextPtr = CreateEdge(0, i);
      LastEdge = LastEdge->nextPtr;
      Temp = Temp->nextState;
    };
  };
  ans = 0;
  while (NewIns != NULL && !ans) {
    PresentIn = NewIns;
    innumber = 0;
    numNewOuts = 0;
    NewOuts = NULL;
    while (PresentIn != NULL) {
      FreeStateList(ReallyNewOuts);
      innumber++;
      ReallyNewOuts = NewRectanglesInto(PrevOuts, PresentIn->data);
      while (ReallyNewOuts != NULL) {
        outnumber = GetNumber(ReallyNewOuts->data, NewOuts);
        if (outnumber == 0) {
          if (numNewOuts == 0) {
            NewOuts = ReallyNewOuts;
            ReallyNewOuts = ReallyNewOuts->nextState;
            NewOuts->nextState = NULL;
            LastNewOut = NewOuts;
            numNewOuts++;
            outnumber = numNewOuts;
          } else {
            LastNewOut->nextState = ReallyNewOuts;
            ReallyNewOuts = ReallyNewOuts->nextState;
            LastNewOut = LastNewOut->nextState;
            LastNewOut->nextState = NULL;
            numNewOuts++;
            outnumber = numNewOuts;
          };
        } else {
          Temp = ReallyNewOuts;
          ReallyNewOuts = ReallyNewOuts->nextState;
          free(Temp);
        }
        EdgeList =
            AppendOrdered(outnumber + numOuts, innumber + numIns, EdgeList);
        edgecount++;
      }
      PresentIn = PresentIn->nextState;
    };
    FreeStateList(PrevIns);
    PrevIns = NewIns;
    i = 1;
    numIns = numIns + innumber;
    previnnumber = numIns;
    numNewIns = 0;
    NewIns = NULL;
    outnumber = 0;
    PresentOut = NewOuts;
    while (PresentOut != NULL) {
      outnumber++;
      ReallyNewIns = NewRectanglesOutOf(PrevIns, PresentOut->data);
      while (ReallyNewIns != NULL) {
        innumber = GetNumber(ReallyNewIns->data, NewIns);
        if (innumber == 0) {
          if (numNewIns == 0) {
            NewIns = ReallyNewIns;
            ReallyNewIns = ReallyNewIns->nextState;
            NewIns->nextState = NULL;
            LastNewIn = NewIns;
            numNewIns++;
            innumber = numNewIns;
          } else {
            LastNewIn->nextState = ReallyNewIns;
            ReallyNewIns = ReallyNewIns->nextState;
            LastNewIn = LastNewIn->nextState;
            LastNewIn->nextState = NULL;
            numNewIns++;
            innumber = numNewIns;
          };
        } else {
          Temp = ReallyNewIns;
          ReallyNewIns = ReallyNewIns->nextState;
          free(Temp);
        }
        EdgeList =
            AppendOrdered(outnumber + numOuts, innumber + numIns, EdgeList);
        edgecount++;
      };
      PresentOut = PresentOut->nextState;
    };
    FreeStateList(PrevOuts);
    PrevOuts = NewOuts;
    NewOuts = NULL;
    SpecialHomology(0, previnnumber);
    if ((EdgeList == NULL) || (EdgeList->start != 0)) {
      ans = 1;
      FreeStateList(NewIns);
      FreeStateList(NewOuts);
      NewIns = NULL;
    } else if (EdgeList->end <= previnnumber) {
      ans = 0;
      FreeStateList(NewIns);
      FreeStateList(NewOuts);
      NewIns = NULL;
    } else {
      numOuts = numOuts + outnumber;
      if (verbose) {
        printf("%d %d %d\n", numIns, numOuts, edgecount);
      }
    };
  };
  return (ans);
}

/**
 * For each point in the permutation count the number of Os
 * that occur to the northeast
 * @param x a permutation
 * @return an int containing the quantity described above
 * @see Os
 */
int NESWpO(char *x) {
  int i = 0, j = 0;
  int ans = 0;
  while (i < ArcIndex) {
    j = i;
    while (j < ArcIndex) {
      if (x[i] <= Os[j]) {
        ans++;
      };
      j++;
    };
    i++;
  };
  return (ans);
}

/**
 * For each O in Os count the number of points in the permutation
 * to the northeast
 * @param x a permutation
 * @return an int containing the quantity described above
 * @see Os
 */
int NESWOp(char *x) {
  int i = 0, j = 0;
  int ans = 0;
  while (i < ArcIndex) {
    j = i + 1;
    while (j < ArcIndex) {
      if (Os[i] < x[j]) {
        ans++;
      };
      j++;
    };
    i++;
  };
  return (ans);
}

/**
 * For each point in the permutation count the number of points in
 * the same permutation that occur to the northeast
 * @param x a permutation
 * @return an int containing the quantity described above
 */
int NESWpp(char *x) {
  int i = 0, j = 0;
  int ans = 0;
  while (i < ArcIndex) {
    j = i;
    while (j < ArcIndex) {
      if (x[i] < x[j]) {
        ans++;
      };
      j++;
    };
    i++;
  };
  return (ans);
}
