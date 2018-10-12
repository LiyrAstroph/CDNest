
from __future__ import division
                                                        
from libc.string cimport *  
from cpython.mem cimport PyMem_Malloc, PyMem_Free, PyMem_Realloc
from cython.operator cimport dereference

cimport mpi4py.MPI as MPI
from mpi4py.libmpi cimport *

import numpy as np
cimport numpy as np
np.import_array()

cdef class sampler:
  cdef int argc               # argc and argv for dnest
  cdef char **argv            #
  cdef DNestFptrSet *fptrset  # function set
  cdef char options_file[200] # options file
  cdef int num_params         # number of paramters
  
  def __cinit__(self, model):
    # Check the model.
    if not hasattr(model, "from_prior") or not callable(model.from_prior):
      raise ValueError("models must have a callable 'from_prior' method")
    if not hasattr(model, "perturb") or not callable(model.perturb):
      raise ValueError("models must have a callable 'perturb' method")
    if not hasattr(model, "log_likelihood") or not callable(model.log_likelihood):
      raise ValueError("models must have a callable 'log_likelihood' method")
    
    set_py_self(model)
    set_size_(model.num_params)

    self.num_params = model.num_params

    strcpy(self.options_file, model.options_file)
    
    cdef int i
    self.argv = <char **>PyMem_Malloc(4*sizeof(char *))
    for i in range(4):
      self.argv[i] = <char *>PyMem_Malloc(200*sizeof(char))
    
    self.argc = 0
    self.argv[self.argc] = 'dnest'
    self.argc += 1
    self.argv[self.argc] = '-s'
    self.argc += 1
    self.argv[self.argc] = 'resetart_dnest.txt'
    self.argc += 1
    
    # setup function set
    self.fptrset = <DNestFptrSet *>PyMem_Malloc(sizeof(DNestFptrSet))  
    self.fptrset.from_prior = py_from_prior
    self.fptrset.log_likelihoods_cal = py_log_likelihood
    self.fptrset.log_likelihoods_cal_initial = py_log_likelihood
    self.fptrset.log_likelihoods_cal_restart = py_log_likelihood
    self.fptrset.perturb = py_perturb
    self.fptrset.print_particle = py_print_particle
    self.fptrset.restart_action = py_restart_action

    return
  
  def __cdealloc__(self):
    cdef int i
    for i in range(4):
      PyMem_Free(self.argv[i])
    
    PyMem_Free(self.argv)
    PyMem_Free(self.fptrset)

    return
  
  def run(self):
    logz = dnest(self.argc, self.argv, self.fptrset, self.num_params, self.options_file)
    return logz