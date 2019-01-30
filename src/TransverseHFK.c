/* Program for calculating the invariant for transverse knots, to
   accompany "Transverse knots distinguished by Knot Floer Homology"
   by L. Ng, P. S. Ozsvath, and D. P. Thurston. The invariant was defined
   in "Legendrian knots, transverse knots, and combinatorial Floer homology"
   by P. S. Ozsvath, Z. Szabo, and D. P. Thurston*/

/* In progress cleanup of the above listed file by
 * Lucas Meyers <lmeye22@lsu.edu>
 */

#include <argp.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
    {0, 'X', "Xs", 0, "List of Xs", 0},
    {0, 'O', "Os", 0, "List of Os", 0},
    {0}};

// Temporary location
static error_t parse_opt(int, char *, struct argp_state *);
static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};

bool verbose = false;
int ArcIndex = -1;
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

ShortEdges AddModTwo(int a, int b, ShortEdges edges);
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

int buildPermutation(char *, char *);
int EqState(State a, State b) { return (!strncmp(a, b, ArcIndex)); }

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  switch (key) {
  case 'v':
    verbose = true;
    break;
  case 'q':
    verbose = false;
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
      argp_failure(state, 0, 0, "Malformated Xs");
      exit(1);
    }
    break;
  case 'O':
    if (-1 == buildPermutation(Os, arg)) {
      argp_failure(state, 0, 0, "Malformated Os");
      exit(1);
    }
    break;
  default:
    break;
  }
  return 0;
}


/**
 * Takes in a string and converts it into a permutation.
 * @param String of the form [_,_,...,_] where _ are integers between 1 and MAX_INDEX
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

int main(int argc, char **argv) {
  argp_parse(&argp, argc, argv, 0, 0, 0);

  char UR[ArcIndex];
  int i;

  if (verbose) {
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

ShortEdges AddModTwoLists(VertexList parents, VertexList kids) {
  VertexList thiskid, thisparent, tempvert;
  ShortEdges thisedge, tempedge, Prev;
  ShortEdges ans;
  if ((parents == NULL) || (kids == NULL)) {
    ans = EdgeList;
  } else {
    thisparent = parents;
    thiskid = kids;
    ans = NULL;
    thisedge = EdgeList;
    while (thisparent != NULL && thisedge != NULL &&
           (thisedge->start == thisparent->data &&
            thisedge->end == thiskid->data)) {
      tempedge = thisedge;
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

ShortEdges AddModTwo(int a, int b, ShortEdges edges) {
  ShortEdges Temp, Prev, curr, ans;
  Prev = edges;
  if ((edges == NULL) || (edges->start > a) ||
      (edges->start == a && edges->end > b)) {
    ans = malloc(sizeof(ShortEdgeNode));
    ans->start = a;
    ans->end = b;
    ans->nextPtr = Prev;
  } else if (edges->start == a && edges->end == b) {
    Temp = edges;
    ans = edges->nextPtr;
    free(Temp);
  } else {
    ans = edges;
    curr = Prev->nextPtr;
    while (curr != NULL &&
           ((curr->start < a) || ((curr->start == a) && (curr->end < b)))) {
      curr = curr->nextPtr;
      Prev = Prev->nextPtr;
    };
    if ((curr != NULL) && (curr->start == a) && (curr->end == b)) {
      Prev->nextPtr = curr->nextPtr;
      free(curr);
    } else {
      Temp = malloc(sizeof(ShortEdgeNode));
      Temp->start = a;
      Temp->end = b;
      Temp->nextPtr = curr;
      Prev->nextPtr = Temp;
    };
  };
  return (ans);
}

ShortEdges PrependEdge(int a, int b, ShortEdges e) {
  ShortEdges newPtr;
  newPtr = malloc(sizeof(ShortEdgeNode));
  newPtr->start = a;
  newPtr->end = b;
  newPtr->nextPtr = e;
  return (newPtr);
}

ShortEdges CreateEdge(int a, int b) {
  ShortEdges newPtr;
  newPtr = malloc(sizeof(ShortEdgeNode));
  newPtr->start = a;
  newPtr->end = b;
  newPtr->nextPtr = NULL;
  return (newPtr);
}

VertexList PrependVertex(int a, VertexList vertices) {
  VertexList newPtr;
  newPtr = malloc(sizeof(Vertex));
  (newPtr->data) = a;
  (newPtr->nextVertex) = vertices;
  return newPtr;
}

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

void Contract(int a, int b) {
  ShortEdges Temp;
  ShortEdges Prev;
  VertexList parents, kids, tempkids, tempparents;
  VertexList LastParent, LastKid;
  Prev = EdgeList;
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

int OrderedQ(VertexList L) {
  int max = 0;
  int t;
  VertexList Temp;
  Temp = L;
  while (Temp != NULL) {
    if (Temp->data < max) {
      printf("XXXX");
      scanf("%d", &t);
      return 0;
    }
    max = Temp->data;
    Temp = Temp->nextVertex;
  };
  return 1; // Temp will always be one if we hit this point.
  // if (Temp==NULL) { return 1; }
}

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

void PrintEdges() {
  ShortEdges Temp;
  Temp = EdgeList;
  while (Temp != NULL) {
    printf("%d %d\n", Temp->start, Temp->end);
    Temp = (Temp->nextPtr);
  };
}

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

void FreeStateList(StateList states) {
  StateList Temp;
  Temp = states;
  while (Temp != NULL) {
    states = states->nextState;
    free(Temp);
    Temp = states;
  };
}

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

void CreateD1Graph(State init) {
  StateList NewIns, NewOuts, LastNewIn, LastNewOut, Temp;
  StateList PrevIns, PrevOuts;
  StateList ReallyNewOuts = NULL, ReallyNewIns = NULL;
  ShortEdges LastEdge;
  ShortEdges PresentEdgeList;
  int innumber;
  int outnumber;
  int i;
  int edgecount = 0;
  int numIns = 0;
  int numOuts = 0;
  int numNewIns = 0;
  int numNewOuts = 0;
  StateList PresentIn, PresentOut;
  PrevOuts = NULL;
  PrevIns = NULL;
  NewIns = FixedWtRectanglesOutOf(1, init);
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
  while (NewIns != NULL) {
    PresentIn = NewIns;
    PresentEdgeList = NULL;
    innumber = 0;
    numNewOuts = 0;
    NewOuts = NULL;
    while (PresentIn != NULL) {
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
        PresentEdgeList = AppendOrdered(outnumber + numOuts, innumber + numIns,
                                        PresentEdgeList);
        edgecount++;
      }
      PresentIn = PresentIn->nextState;
    };
    FreeStateList(PrevIns);
    PrevIns = NewIns;
    numIns = numIns + innumber;
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
        PresentEdgeList = AppendOrdered(outnumber + numOuts, innumber + numIns,
                                        PresentEdgeList);
        edgecount++;
      };
      PresentOut = PresentOut->nextState;
    };
    if (PresentEdgeList != NULL) {
      LastEdge->nextPtr = PresentEdgeList;
      LastEdge = PresentEdgeList;
      while (LastEdge->nextPtr != NULL) {
        LastEdge = LastEdge->nextPtr;
      };
      PresentEdgeList = NULL;
    };
    FreeStateList(PrevOuts);
    PrevOuts = NewOuts;
    NewOuts = NULL;
    numOuts = numOuts + outnumber;
    if (verbose) {
      printf("%d %d %d\n", numIns, numOuts, edgecount);
    }
  };
}

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
