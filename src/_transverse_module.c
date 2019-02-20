#include <Python.h>
#include "TransverseHFK.c"

static PyObject* null_homologous_D0Q_py(PyObject* self, PyObject* args, PyObject* keywds) {
  PyObject* py_Xs;
  PyObject* py_Os;
  PyObject* py_state;
  PyObject* result;

  Grid_t G;
  char[MAX_INDEX] state;

  const char** keyword_list = {"state","Xs", "Os", 0};

  if(!(PySequence_Check(py_Xs) && PySequenceCheck(py_Os))) {
    PyErr_SetString()
  }
    
  return NULL;
}

static PyObject* null_homologous_D1Q_py(PyObject* self, PyObject* args, PyObject* keywds) {
  return NULL;
}
