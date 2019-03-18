#include <Python.h>
#include <string.h>
#include "TransverseHFK.c"

static PyObject* error = NULL;
static PyObject* out_stream = NULL;

int print_py(const char * format, ...) {
  int ret = -1;
  if (NULL == out_stream) {
    return ret;
  }
  char* s = malloc(201*sizeof(char));
  va_list args;
  va_start(args, format);
  vsnprintf(s, 200, format, args);

  PyObject_CallMethod(out_stream, "write", "(s)", s);

  ret = strlen(s);
  free(s);
  
  return ret;
}

static PyObject* null_homologous_D0Q_py(PyObject* self, PyObject* args, PyObject* keywds) {
  PyObject* py_Xs = NULL;
  PyObject* py_Os = NULL;
  PyObject* py_state = NULL;
  int py_verbosity = 0;
  PyObject* py_out_stream = NULL;
  
  Grid_t G;
  char state[MAX_INDEX];

  G.arc_index = 0;
  for(int i = 0; i < MAX_INDEX; ++i) {
    state[i] = 0;
    G.Xs[i] = 0;
    G.Os[i] = 0;
  }

  const char* keyword_list[] = {"state","Xs", "Os", "out_stream", "verbosity", 0};

  if(!(PyArg_ParseTupleAndKeywords(args, keywds, "OOOOi:null_homologous_D0Q", (char**)keyword_list, &py_state, &py_Xs, &py_Os, &py_out_stream, &py_verbosity))) {
    return NULL;
  }
      
  if(!(PySequence_Check(py_Xs) && PySequence_Check(py_Os) && PySequence_Check(py_state))) {
    PyErr_SetString(error, "The state, Xs, and Os must be sequences.");
    return NULL;
  }

  G.arc_index = PySequence_Length(py_Xs);

  if(PySequence_Length(py_Os) != G.arc_index || PySequence_Length(py_state) != G.arc_index) {
    PyErr_SetString(error, "The state, Xs, and Os must be the same length");
      return NULL;
  }

  if(G.arc_index > 30 || G.arc_index < 2) {
    PyErr_SetString(error, "The grid size must be between 2 and 30");
    return NULL;
  }

  int failed = 0;
  
  for(int i = 0; i < G.arc_index; ++i) {
    PyObject *elem = PySequence_GetItem(py_Xs, i);

    if(NULL == elem || !PyInt_Check(elem)) {
      failed = 1;
      break;
    }
    G.Xs[i] = (int)PyInt_AS_LONG(elem);
    Py_DECREF(elem);

    elem = PySequence_GetItem(py_Os, i);

    if(NULL == elem || !PyInt_Check(elem)) {
      failed = 1;
      break;
    }
    G.Os[i] = (int)PyInt_AS_LONG(elem);
    Py_DECREF(elem);

    elem = PySequence_GetItem(py_state, i);

    if(NULL == elem || !PyInt_Check(elem)) {
      failed = 1;
      break;
    }
    state[i] = (char)PyInt_AS_LONG(elem);
    Py_DECREF(elem);    
  }

  if(failed || !is_grid(&G) || !is_state(state,&G)) {
    PyErr_SetString(error, "state, Xs, and Os must be lists containing [1,...,N] exactly once with no matching indices between Xs and Os");
    return NULL;
  }
  else if (verbosity < 0 || 2 < verbosity) {
    PyErr_SetString(error, "verbosity must be passed an integer 0, 1, or 2.");
    return NULL;
  }
  else {
    set_verbosity(py_verbosity);
  }

  if(NULL == py_out_stream) {
    PyErr_SetString(error, "An out stream must be specified.");
    return NULL;
  }
  else if (!PyObject_HasAttrString(py_out_stream, "write")){
    PyErr_SetString(error, "The out stream must implement the write method.");
    return NULL;
  }
  else {
    out_stream = py_out_stream;
  }
  
  if(null_homologous_D0Q(state,&G)) {
    Py_RETURN_TRUE;
  }
  else {
    Py_RETURN_FALSE;
  }
}

static PyObject* null_homologous_D1Q_py(PyObject* self, PyObject* args, PyObject* keywds) {
  PyObject* py_Xs = NULL;
  PyObject* py_Os = NULL;
  PyObject* py_state = NULL;
  int py_verbosity = 0;
  PyObject* py_out_stream = NULL;
  
  Grid_t G;
  char state[MAX_INDEX];

  G.arc_index = 0;
  for(int i = 0; i < MAX_INDEX; ++i) {
    state[i] = 0;
    G.Xs[i] = 0;
    G.Os[i] = 0;
  }

  const char* keyword_list[] = {"state","Xs", "Os", "out_stream", "verbosity", 0};

  if(!(PyArg_ParseTupleAndKeywords(args, keywds, "OOOOi:null_homologous_D1Q", (char**)keyword_list, &py_state, &py_Xs, &py_Os, &py_out_stream, &py_verbosity))) {
    return NULL;
  }
      
  if(!(PySequence_Check(py_Xs) && PySequence_Check(py_Os) && PySequence_Check(py_state))) {
    PyErr_SetString(error, "The state, Xs, and Os must be sequences.");
    return NULL;
  }

  G.arc_index = PySequence_Length(py_Xs);

  if(PySequence_Length(py_Os) != G.arc_index || PySequence_Length(py_state) != G.arc_index) {
    PyErr_SetString(error, "The state, Xs, and Os must be the same length");
      return NULL;
  }

  if(G.arc_index > 30 || G.arc_index < 2) {
    PyErr_SetString(error, "The grid size must be between 2 and 30");
    return NULL;
  }

  int failed = 0;
  
  for(int i = 0; i < G.arc_index; ++i) {
    PyObject *elem = PySequence_GetItem(py_Xs, i);

    if(NULL == elem || !PyInt_Check(elem)) {
      failed = 1;
      break;
    }
    G.Xs[i] = (int)PyInt_AS_LONG(elem);
    Py_DECREF(elem);

    elem = PySequence_GetItem(py_Os, i);

    if(NULL == elem || !PyInt_Check(elem)) {
      failed = 1;
      break;
    }
    G.Os[i] = (int)PyInt_AS_LONG(elem);
    Py_DECREF(elem);

    elem = PySequence_GetItem(py_state, i);

    if(NULL == elem || !PyInt_Check(elem)) {
      failed = 1;
      break;
    }
    state[i] = (char)PyInt_AS_LONG(elem);
    Py_DECREF(elem);    
  }

  if(failed || !is_grid(&G) || !is_state(state,&G)) {
    PyErr_SetString(error, "state, Xs, and Os must be lists containing [1,...,N] exactly once with no matching indices between Xs and Os");
    return NULL;
  }
  else if (verbosity < 0 || 2 < verbosity) {
    PyErr_SetString(error, "verbosity must be passed an integer 0, 1, or 2.");
    return NULL;
  }
  else {
    set_verbosity(py_verbosity);
  }

  if(NULL == py_out_stream) {
    PyErr_SetString(error, "An out stream must be specified.");
    return NULL;
  }
  else if (!PyObject_HasAttrString(py_out_stream, "write")){
    PyErr_SetString(error, "The out stream must implement the write method.");
    return NULL;
  }
  else {
    out_stream = py_out_stream;
  }
  
  if(null_homologous_D1Q(state,&G)) {
    Py_RETURN_TRUE;
  }
  else {
    Py_RETURN_FALSE;
  }
}

static char null_homologous_D0Q_doc[] =
  "Returns true if the supplied state is null-homologous for the \
corresponding grid.\n\
\n\
Parameters\n\
----------\n\
state: [int]\n\
    a grid state that represents a homology class\n\
Xs: [int]\n\
    int list specifying the Xs of the grid\n\
Os: [int]\n\
    int list specifying the Os of the grid\n\
out_stream : stream\n\
    An object with a .write method that is used for inner\n\
    printing by the methods. Does nothing if verbosity is 0.\n\
verbosity : int\n\
    An integer specifying the verbosity of the methods. Must\n\
    be 0, 1, or 2. 0 will print no information and 2 will print\n\
    the most. Defaults to 0.\n\n\
\
Note: Xs, Os, and state must be permutations {1,..,N}\n\
where Xs and Os have no overlapping values.";
static char null_homologous_D1Q_doc[] =
  "Returns true if the supplied state is null-homologous after\n\
the d_1 map is applied for the corresponding grid.\n\
\n\
Parameters\n\
----------\n\
state: [int]\n\
    a grid state that represents a homology class\n\
Xs: [int]\n\
    int list specifying the Xs of the grid\n\
Os: [int]\n\
    int list specifying the Os of the grid\n\
out_stream : stream\n\
    An object with a .write method that is used for inner\n\
    printing by the methods. Does nothing if verbosity is 0.\n\
verbosity : int\n\
    An integer specifying the verbosity of the methods. Must\n\
    be 0, 1, or 2. 0 will print no information and 2 will print\n\
    the most. Defaults to 0.\n\n\
\
Note: Xs, Os, and state must be permutations {1,..,N}\n\
where Xs and Os have no overlapping values.";


static PyMethodDef _tHFK_methods[] = {
                                               {"null_homologous_D0Q", (PyCFunction)null_homologous_D0Q_py, METH_VARARGS|METH_KEYWORDS, null_homologous_D0Q_doc},
                                               {"null_homologous_D1Q", (PyCFunction)null_homologous_D1Q_py, METH_VARARGS|METH_KEYWORDS, null_homologous_D1Q_doc},
                                               {NULL, NULL}
};

PyMODINIT_FUNC init_tHFK(void) {
  PyObject *m, *d;
  const char *tHFK_error_name = "tHFK_Error";
  const char *tHFK_dot_error = "tHFK.error";

  m = Py_InitModule("_tHFK", _tHFK_methods);

  d = PyModule_GetDict(m);
  error = PyErr_NewException((char *) tHFK_dot_error, NULL, NULL);
  PyDict_SetItemString(d, tHFK_error_name, error);

  print_ptr = print_py;
}
