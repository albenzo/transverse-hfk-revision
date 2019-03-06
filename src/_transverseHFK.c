#include <Python.h>
#include "TransverseHFK.c"

static PyObject* error;

static PyObject* null_homologous_D0Q_py(PyObject* self, PyObject* args, PyObject* keywds) {
  PyObject* py_Xs;
  PyObject* py_Os;
  PyObject* py_state;
  PyObject* result;

  Grid_t G;
  char[MAX_INDEX] state;

  const char** keyword_list = {"state","Xs", "Os", 0};

  if(!(PySequence_Check(py_Xs) && PySequenceCheck(py_Os))) {
    PyErr_SetString(error, "The Xs and Os must be sequences.");
    return NULL;
  }

  G.arc_index = PySequence_Length(py_Xs);

  if(PySequence_Length(py_Os)) {
    PyErr_SetString(error, "The Xs and Os must be the same length");
      return NULL;
  }

  if(G.arc_index > 30 || G.arc_index < 2) {
    PyErr_SetString(error, "The grid size must be between 2 and 30");
    return NULL;
  }

  
    
  return NULL;
}

static PyObject* null_homologous_D1Q_py(PyObject* self, PyObject* args, PyObject* keywds) {
  return NULL;
}
