#
# C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
#
# Wrapped to Python
#
# Yan-Rong Li, liyanrong@mail.ihep.ac.cn
# Jun 30, 2016
#
#

from __future__ import division
                                                        
from libc.string cimport *  
from cpython.mem cimport PyMem_Malloc, PyMem_Free, PyMem_Realloc
from cython.operator cimport dereference

#cimport python_unicode
from mpi4py import MPI
cimport mpi4py.MPI as MPI
from mpi4py.libmpi cimport *

import numpy as np
cimport numpy as np
np.import_array()

cdef class sampler:
  """
  a dnest mcmc sampler.

  :param model
    a model that specifies options and function set for dnest.
  """
  cdef int argc               # argc and argv for dnest
  cdef char **argv            #
  cdef DNestFptrSet *fptrset  # function set
  cdef char options_file[200] # options file
  cdef int num_params         # number of paramters
  cdef char sample_tag[200]
  cdef char sample_postfix[200]
  cdef char sample_dir[200]
  cdef int rank, size
  cdef int num_particles, max_num_saves, max_num_levels
  cdef int new_level_interval, save_interval, thread_steps
  cdef double beta, Lambda
  
  def __cinit__(self, model, sample_dir="./", sample_tag="", sample_postfix="", 
                new_level_interval = 10, save_interval = 10, thread_steps = 10, 
                num_particles=1, max_num_saves = 10000, max_num_levels = 10,
                Lambda = 10, beta = 100):
    
    # check model
    if not hasattr(model, "from_prior") or not callable(model.from_prior):
      raise ValueError("models must have a callable 'from_prior' method")
    if not hasattr(model, "perturb") or not callable(model.perturb):
      raise ValueError("models must have a callable 'perturb' method")
    if not hasattr(model, "log_likelihood") or not callable(model.log_likelihood):
      raise ValueError("models must have a callable 'log_likelihood' method")
    
    self.size  = MPI.COMM_WORLD.Get_size()
    self.rank  = MPI.COMM_WORLD.Get_rank()

    # setup options
    self.num_params = model.num_params
    # sample tag
    py_byte_string = sample_tag.encode('UTF-8')
    strcpy(self.sample_tag, py_byte_string)
    # sample postfix
    py_byte_string = sample_postfix.encode('UTF-8')
    strcpy(self.sample_postfix, py_byte_string)
    # default directiory 
    py_byte_string =sample_dir.encode('UTF-8')
    strcpy(self.sample_dir, py_byte_string)
    # options file
    strcpy(self.options_file, self.sample_dir)
    strcat(self.options_file, "/OPTIONS")
    strcat(self.options_file, self.sample_tag)
    
    # options for CDNest
    cdef int num
    num = np.max((self.num_params * 2 * self.num_particles, 10))
    self.new_level_interval = np.max((new_level_interval, num * self.size))
    self.save_interval = np.max((save_interval, num*self.size))
    self.thread_steps = np.max((thread_steps, num))
    self.num_particles = num_particles
    self.max_num_saves = max_num_saves
    self.max_num_levels = max_num_levels
    self.Lambda = Lambda 
    self.beta = beta
    
    # setup functions from model
    set_py_self(model)
    set_size_(model.num_params)
    
    # setup argc and argv
    cdef int i
    cdef int narg = 8
    self.argv = <char **>PyMem_Malloc(narg*sizeof(char *))
    for i in range(narg):
      self.argv[i] = <char *>PyMem_Malloc(200*sizeof(char))
    
    self.argc = 0
    self.argv[self.argc] = 'dnest'
    self.argc += 1
    self.argv[self.argc] = '-s'
    self.argc += 1
    strcpy(self.argv[self.argc], self.sample_dir)
    strcat(self.argv[self.argc], '/restart_dnest.txt')
    self.argc += 1
    self.argv[self.argc] = '-g'
    self.argc += 1
    strcpy(self.argv[self.argc], self.sample_tag)
    self.argc += 1
    self.argv[self.argc] = '-x'
    self.argc += 1
    strcpy(self.argv[self.argc], self.sample_postfix)
    self.argc += 1
    
    # setup function set
    self.fptrset = dnest_malloc_fptrset(); 
    self.fptrset.from_prior = py_from_prior
    self.fptrset.log_likelihoods_cal = py_log_likelihood
    self.fptrset.log_likelihoods_cal_initial = py_log_likelihood
    self.fptrset.log_likelihoods_cal_restart = py_log_likelihood
    self.fptrset.perturb = py_perturb
    self.fptrset.print_particle = py_print_particle
    self.fptrset.restart_action = py_restart_action

    self.save_options_file()
    return
  
  def __cdealloc__(self):
    cdef int i
    for i in range(4):
      PyMem_Free(self.argv[i])
    PyMem_Free(self.argv)

    dnest_free_fptrset(self.fptrset)
    return
  
  def get_sample_tag(self):
    return self.sample_tag.decode('UTF-8')
  
  def get_sample_postfix(self):
    return self.sample_postfix.decode('UTF-8')

  def get_sample_dir(self):
    return self.sample_dir.decode('UTF-8')

  def save_options_file(self):
    """
    write out an OPTIONS file that CDNest needs to read in
    """
    fp = open(self.options_file.decode('UTF-8'), "w")
    fp.write("# File containing parameters for DNest\n")
    fp.write("# Put comments at the top, or at the end of the line\n")
    fp.write("# Do not change the order of lines\n")
    fp.write("# Lines beginning with '#' are regarded as comments\n\n\n")
    fp.write("%d  # Number of particles\n"%(self.num_particles))
    fp.write("%d  # New level interval\n"%(self.new_level_interval))
    fp.write("%d  # Save interval\n"%(self.save_interval))
    fp.write("%d  # ThreadSteps\n"%(self.thread_steps))
    fp.write("%d  # Maximum number of levels\n"%(self.max_num_levels))
    fp.write("%d  # Backtracking scale length\n"%(self.Lambda))
    fp.write("%d  # Strength of effect to force histogram to equal push\n"%(self.beta))
    fp.write("%d  # Maximum number of saves\n"%(self.max_num_saves))
    fp.close()
  
  def run(self):
    """
    run dnest sampling.

    the value of evidence is returned.
    """
    logz = dnest(self.argc, self.argv, self.fptrset, self.num_params, self.sample_dir, self.options_file)
    return logz