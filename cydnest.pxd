#
# C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
#
# Wrapped to Python
#
# Yan-Rong Li, liyanrong@mail.ihep.ac.cn
# Jun 30, 2016
#
#

#!/usr/bin/python
#cython: initializedcheck=False, boundscheck=False, wraparound=False, cdivision=True, profile=False

from libc.stdio cimport *                                                                

# include declarations in PyFuncs.h
cdef extern from "PyFuncs.h":
  void set_py_self (object py_self) 
  void set_size_(int size)
  void py_from_prior(void *pm)
  double py_perturb(void *pm)
  double py_log_likelihood(void *pm) 
  void py_print_particle(FILE *fp, void *pm)
  void py_restart_action(int iflag)

# typedef functions for dnest
ctypedef void (*from_prior_type)(void *)
ctypedef double (*log_likelihoods_cal_type)(void *)
ctypedef double (*perturb_type)(void *)
ctypedef void (*print_particle_type)(FILE *, void*)
ctypedef void (*restart_action_type)(int)

# include declarations from dnestvars.h
cdef extern from "dnestvars.h":
  double dnest_randh()
  double dnest_rand()
  double dnest_randn()
  double wrap(double *x, double xmin, double xmax)
  int dnest_rand_int(int size);

  # function set
  ctypedef struct DNestFptrSet:
    from_prior_type from_prior
    log_likelihoods_cal_type log_likelihoods_cal
    log_likelihoods_cal_type log_likelihoods_cal_initial
    log_likelihoods_cal_type log_likelihoods_cal_restart
    perturb_type perturb
    print_particle_type print_particle
    restart_action_type restart_action
  
  DNestFptrSet * dnest_malloc_fptrset()
  void dnest_free_fptrset(DNestFptrSet *fptrset)
  double dnest(int argc, char **argv, DNestFptrSet *fptrset,  int num_params, char *optfile)