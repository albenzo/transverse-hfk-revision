#include <Python.h>
#include "TransverseHFK.c"

static PyObject* error;

static PyObject* null_homologous_D0Q_py(PyObject* self, PyObject* args, PyObject* keywds) {
  PyObject* py_Xs;
  PyObject* py_Os;
  PyObject* py_state;
  
  Grid_t G;
  char state[MAX_INDEX];

  G.arc_index = 0;
  for(int i = 0; i < MAX_INDEX; ++i) {
    state[i] = 0;
    G.Xs[i] = 0;
    G.Os[i] = 0;
  }

  const char* keyword_list[] = {"state","Xs", "Os", 0};

  if(!(PyArg_ParseTupleAndKeywords(args, keywds, "OOO|i:null_homologous_D0Q", (char**)keyword_list, &py_state, &py_Xs, &py_Os))) {
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
  /*
  if(NULL == py_verbose) {
    set_verbosity(false);
  }
  else if (!PyBool_Check(py_verbose)) {
    PyErr_SetString(error, "verbosity must be passed a boolean value");
    return NULL;
  }
  else {
    set_verbosity((int)PyInt_AS_LONG(py_verbose));
  }
  */
  if(null_homologous_D0Q(state,&G)) {
    Py_RETURN_TRUE;
  }
  else {
    Py_RETURN_FALSE;
  }
}

static PyObject* null_homologous_D1Q_py(PyObject* self, PyObject* args, PyObject* keywds) {
  PyObject* py_Xs;
  PyObject* py_Os;
  PyObject* py_state;
  
  Grid_t G;
  char state[MAX_INDEX];

  G.arc_index = 0;
  for(int i = 0; i < MAX_INDEX; ++i) {
    state[i] = 0;
    G.Xs[i] = 0;
    G.Os[i] = 0;
  }

  const char* keyword_list[] = {"state","Xs", "Os", 0};

  if(!(PyArg_ParseTupleAndKeywords(args, keywds, "OOO|i:null_homologous_D0Q", (char**)keyword_list, &py_state, &py_Xs, &py_Os))) {
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
  /*
  if(NULL == py_verbose) {
    set_verbosity(false);
  }
  else if (!PyBool_Check(py_verbose)) {
    PyErr_SetString(error, "verbosity must be passed a boolean value");
    return NULL;
  }
  else {
    set_verbosity((int)PyInt_AS_LONG(py_verbose));
  }
  */
  if(null_homologous_D1Q(state,&G)) {
    Py_RETURN_TRUE;
  }
  else {
    Py_RETURN_FALSE;
  }
}

static char null_homologous_D0Q_doc[] =
  "Returns true if the supplied state is null-homologous for the\
corresponding grid.\
\
Parameters\
----------\
state: [int]\
    a grid state that represents a homology class\
Xs: [int]\
    int list specifying the Xs of the grid\
Os: [int]\
    int list specifying the Os of the grid\
Note: Xs, Os, and state must be permutations {1,..,N}\
where Xs and Os have no overlapping values.";
static char null_homologous_D1Q_doc[] =
  "Returns true if the supplied state is null-homologous after\
the d_1 map is applied for the corresponding grid.\
\
Parameters\
----------\
state: [int]\
    a grid state that represents a homology class\
Xs: [int]\
    int list specifying the Xs of the grid\
Os: [int]\
    int list specifying the Os of the grid\
Note: Xs, Os, and state must be permutations {1,..,N}\
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
}
