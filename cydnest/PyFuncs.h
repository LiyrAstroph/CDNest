/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */

#ifndef _PyFuncs
#define _PyFuncs

//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <stdio.h>
#include <stdlib.h>
#include <Python.h>
#include <string.h>
#include <numpy/arrayobject.h>

PyObject* py_self_;
int size_;

void set_py_self(PyObject* py_self) 
{ 
  py_self_ = py_self; 
  return;
}

void set_size_(int size)
{
  size_ = size;
  return;
}

PyObject* get_npy_coords (void *params) 
{
  double *pm = (double *)params;
  
  npy_intp shape[] = {size_};
  PyObject* c = PyArray_SimpleNew(1, shape, NPY_DOUBLE);
  if (c == NULL) 
  {
    PyErr_Print();
    exit(0);
  }
  double* data = (double*)PyArray_DATA(c);
  for (int i = 0; i < size_; ++i) data[i] = pm[i];
  return c;
}

void py_from_prior(void *params)
{
  double *pm = (double *)params;
  
  // Call the Python method and get the Python return value.
  PyObject* result = PyObject_CallMethod(py_self_, "from_prior", "");
  if (result == NULL) 
  {
    PyErr_Print();
    Py_XDECREF(result);
    exit(0);
    return;
  }
  
  // Parse that return value as a numpy array.
  PyObject* rarray = PyArray_FROM_OTF(result, NPY_DOUBLE, NPY_IN_ARRAY);
  if (result == NULL || (int)PyArray_NDIM(rarray) != 1) 
  {
    Py_DECREF(result);
    Py_XDECREF(rarray);
    exit(0);
    return;
  }
  
  double* data = (double*)PyArray_DATA(rarray);
  for (int i = 0; i < size_; ++i)
  { 
    pm[i] = data[i];
  }

  Py_DECREF(result);
  Py_DECREF(rarray);
  return;
}

double py_perturb (void *params) 
{
  double *pm = (double *)params;
  
  PyObject* c = get_npy_coords(params);
  
  // Call the Python method and get the Python return value.
  PyObject* result = PyObject_CallMethod(py_self_, "perturb", "O", c);
  
  if (result == NULL || PyErr_Occurred() != NULL) 
  {
    PyErr_Print();
    Py_DECREF(c);
    Py_XDECREF(result);
    exit(0);
    return 0.0;
  }

  double log_H = PyFloat_AsDouble(result);
  Py_DECREF(result);
  if (PyErr_Occurred() != NULL) 
  {
    PyErr_Print();
    Py_DECREF(c);
    exit(0);
    return 0.0;
  }

  double* data = (double*)PyArray_DATA(c);
  for (int i = 0; i < size_; ++i) pm[i] = data[i];
  
  Py_DECREF(c);
  return log_H;
}

 // Likelihood function
double py_log_likelihood(void *params) 
{
  if (size_ == 0) return 0.0;
    
  PyObject* c = get_npy_coords(params);
  
  // Call the Python method and get the Python return value.
  PyObject* result = PyObject_CallMethod(py_self_, "log_likelihood", "O", c);
  Py_DECREF(c);
  if (result == NULL) 
  {
    PyErr_Print();
    Py_XDECREF(result);
    exit(0);
    return -INFINITY;
  }

  // Parse as double.
  double log_like = PyFloat_AsDouble(result);
  if (PyErr_Occurred() != NULL) 
  {
    PyErr_Print();
    Py_DECREF(result);
    exit(0);
    return -INFINITY;
  }
  
  return log_like;
}

void py_print_particle(FILE *fp, void *params)
{
  int i;
  double *pm = (double *)params;
  for(i=0; i<size_; i++)
  {
    fprintf(fp, "%f ", pm[i] );
  }
  fprintf(fp, "\n"); 
  return;
}

void py_restart_action(int iflag)
{
  return;
}

void py_get_param_range(PyObject* py_param_range, double *param_range)
{
  PyListObject *plist = py_param_range;
  PyListObject *pitem = NULL;

  int i; 
  for(i=0; i<PyList_Size(plist); i++)
  {
    pitem = PyList_GetItem(plist, i);
    param_range[i*2+0] = PyFloat_AsDouble(PyList_GetItem(pitem, 0));
    param_range[i*2+1] = PyFloat_AsDouble(PyList_GetItem(pitem, 1));
  }
  return;
}

void py_get_prior_info(PyObject* py_prior_info, double *prior_info)
{
  PyListObject *plist = py_prior_info;
  PyListObject *pitem = NULL;

  int i; 
  for(i=0; i<PyList_Size(plist); i++)
  {
    pitem = PyList_GetItem(plist, i);
    prior_info[i*2+0] = PyFloat_AsDouble(PyList_GetItem(pitem, 0));
    prior_info[i*2+1] = PyFloat_AsDouble(PyList_GetItem(pitem, 1));
  }
  return;
}

void py_get_prior_type(PyObject* py_prior_type, int *prior_type)
{
  PyListObject *plist = py_prior_type;
  PyListObject *pitem = NULL;
  char pystr[100];
  int i;
  for(i=0; i<PyList_Size(plist); i++)
  {
    pitem = PyList_GetItem(plist, i);
    strcpy(pystr, PyUnicode_DATA(pitem));
    if(strcmp(pystr, "Gaussian") == 0 || strcmp(pystr, "gaussian") == 0)
    {
      prior_type[i]=0;
    }
    else if(strcmp(pystr, "Log") == 0 || strcmp(pystr, "log") == 0)
    {
      prior_type[i]=1;
    }
    else if(strcmp(pystr, "Uniform") == 0 || strcmp(pystr, "uniform") == 0)
    {
      prior_type[i]=2;
    }
    else 
    {
      printf("Unable to identify prior type: %s\n", pystr);
      exit(0);
    }
  }
  return;
}
#endif
