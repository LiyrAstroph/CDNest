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
  cdef double *param_range, *prior_info
  cdef int *prior_type
  cdef void *args
  cdef char sample_tag[200]
  cdef char sample_postfix[200]
  cdef char sample_dir[200]
  cdef int rank, size
  cdef int num_particles, max_num_saves, max_num_levels
  cdef int new_level_interval, save_interval, thread_steps
  cdef int thread_steps_factor, new_level_interval_factor, save_interval_factor
  cdef double beta, Lambda, max_ptol
  
  def __cinit__(self, model, sample_dir="./", sample_tag="", sample_postfix="", 
                num_particles=1, thread_steps_factor = 10, 
                max_num_saves = 10000, max_num_levels = 0,
                new_level_interval_factor = 2, save_interval_factor = 2,
                Lambda = 10, beta = 100, ptol = 0.1):
    
    cdef int i

    # check model
    if not hasattr(model, "num_params"):
      raise ValueError("model must has a member variable 'num_params' to specify number of parameters")
    if not hasattr(model, "log_likelihood") or not callable(model.log_likelihood):
      raise ValueError("models must have a callable 'log_likelihood' method")
    
    self.size  = MPI.COMM_WORLD.Get_size()
    self.rank  = MPI.COMM_WORLD.Get_rank()

    # setup options
    self.num_params = model.num_params

    # prior 
    self.param_range = <double *>PyMem_Malloc(self.num_params*2*sizeof(double))
    self.prior_info = <double *>PyMem_Malloc(self.num_params*2*sizeof(double))
    self.prior_type = <int *>PyMem_Malloc(self.num_params*sizeof(int))
    for i in range(self.num_params):
      self.param_range[i*2+0] = model.param_range[i][0]
      self.param_range[i*2+1] = model.param_range[i][1]
      self.prior_info[i*2+0] = model.prior_info[i][0]
      self.prior_info[i*2+1] = model.prior_info[i][1]
      if model.prior_type[i] == 'Uniform' or model.prior_type[i] == 'uniform':
        self.prior_type[i] = UNIFORM
      elif model.prior_type[i] == 'Gaussian' or model.prior_type[i] == 'gaussian':
        self.prior_type[i] = GAUSSIAN
      elif model.prior_type[i] == 'Log' or model.prior_type[i] == 'log':
        self.prior_type[i] = LOG
      else:
        raise ValueError("Unable to identify prior type.")

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
    if hasattr(model, "options_file"):
      py_byte_string = model.options_file.encode('UTF-8')
      strcpy(self.options_file, py_byte_string)
    else:
      strcpy(self.options_file, self.sample_dir)
      strcat(self.options_file, "/OPTIONS")
      strcat(self.options_file, self.sample_tag)
      # options for CDNest
      self.num_particles = num_particles
      self.thread_steps_factor = thread_steps_factor
      self.new_level_interval_factor = new_level_interval_factor
      self.save_interval_factor = save_interval_factor
      self.thread_steps = thread_steps_factor * num_particles * self.num_params
      self.new_level_interval = self.thread_steps * self.size * new_level_interval_factor
      self.save_interval = self.thread_steps * self.size * save_interval_factor
      self.max_num_saves = max_num_saves
      self.max_num_levels = max_num_levels
      self.Lambda = Lambda 
      self.beta = beta
      self.max_ptol = ptol
      # save option file
      self.save_options_file()
    
    # setup functions from model
    set_py_self(model)
    set_size_(model.num_params)
    
    # setup argc and argv
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
    self.fptrset.log_likelihoods_cal = py_log_likelihood
    self.fptrset.log_likelihoods_cal_initial = py_log_likelihood
    self.fptrset.log_likelihoods_cal_restart = py_log_likelihood
    self.fptrset.print_particle = py_print_particle
    self.fptrset.restart_action = py_restart_action
    if hasattr(model, "from_prior") and callable(model.from_prior):
      self.fptrset.from_prior = py_from_prior
    if hasattr(model, "perturb") and callable(model.perturb):
      self.fptrset.perturb = py_perturb
    
    return
  
  def __cdealloc__(self):
    cdef int i
    for i in range(4):
      PyMem_Free(self.argv[i])
    PyMem_Free(self.argv)

    dnest_free_fptrset(self.fptrset)

    PyMem_Free(self.param_range)
    PyMem_Free(self.prior_type)
    PyMem_Free(self.prior_info)
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
    fp.write("# Lines beginning with '#' are regarded as comments\n\n\n")
    fp.write("NumberParticles          %d  # Number of particles\n"%(self.num_particles))
    fp.write("NewLevelIntervalFactor   %d  # New level interval factor\n"%(self.new_level_interval_factor))
    fp.write("SaveIntervalFactor       %d  # Save interval factor\n"%(self.save_interval_factor))
    fp.write("ThreadStepsFactor        %d  # ThreadSteps factor\n"%(self.thread_steps_factor))
    fp.write("MaxNumberLevels          %d  # Maximum number of levels\n"%(self.max_num_levels))
    fp.write("BacktrackingLength       %.1f  # Backtracking scale length\n"%(self.Lambda))
    fp.write("StrengthEqualPush        %.1f  # Strength of effect to force histogram to equal push\n"%(self.beta))
    fp.write("MaxNumberSaves           %d    # Maximum number of saves\n"%(self.max_num_saves))
    fp.write("PTol                     %.1e  # Likelihood tolerance in loge"%(self.max_ptol))
    fp.close()
  
  def run(self):
    """
    run dnest sampling.

    the value of evidence is returned.
    """
    logz = dnest(self.argc, self.argv, self.fptrset, self.num_params, 
                 self.param_range, self.prior_type, self.prior_info, 
                 self.sample_dir, self.options_file, self.args)
    return logz