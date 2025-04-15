*****
Usage
*****
**CDNest** can be called in both C/C++ and Python languages.

C/C++
=====
The statements in C/C++ look like 

.. code-block:: C 
  
  /* include dnest header */
  #include "dnest.h"
  
  /* functions used for dnest */
  DNestFptrSet *fptrset_thismodel;
  /* Bayesian evidence */
  double logz;


  /* allocate function memory used by cdnest */
  fptrset_thismodel = dnest_malloc_fptrset();

  /* setup functions used for cdnest*/
  fptrset_thismodel->from_prior = from_prior_thismodel;
  fptrset_thismodel->log_likelihoods_cal = log_likelihoods_cal_thismodel;
  fptrset_thismodel->perturb = perturb_thismodel;
  fptrset_thismodel->print_particle = print_particle_thismodel;
  fptrset_thismodel->restart_action = restart_action_model1;
  
  /* run dnest */
  logz=dnest(argc, argv, fptrset_thismodel, num_params, NULL, NULL, NULL, "./", "OPTIONS1", NULL, NULL);
    
  /* free function memory */
  dnest_free_fptrset(fptrset_thismodel);

Here, ``OPTIONS1`` specifies the name of the OPTIONS file, which usually looks like,::

  # File containing parameters for DNest
  # Lines beginning with '#' are regarded as comments
  
  
  NumberParticles          1      # Number of particles
  NewLevelIntervalFactor   2.00   # New level interval factor
  ThreadStepsFactor        10.00  # ThreadSteps factor
  MaxNumberSaves           20000    # Maximum number of saves
  PTol                     1.0e-01  # Likelihood tolerance in loge

  # Full options and their default values (if not specified) are:
  # MaxNumberSaves           10000 #maximum number of saving
  # PTol                     0.1   #likelihood tolerance in loge
  # NumberParticles          1  #number of particles
  # NewLevelIntervalFactor   2  #new level interval
  # SaveIntervalFactor       2  #particular saving interval
  # ThreadStepsFactor        10 #thread steps before communications between cores
  # MaxNumberLevels          0  #maximum number of levels; unlimited for 0
  # BacktrackingLength       10.0  #backforward tracking length (lambda)
  # StrengthEqualPush        100.0 #strength to force equal push (beta)

All lines are optional and if not specified, the default value will be used. 
 
  * ``NumberParticles`` constrols the number of particles for each core.
  * ``MaxNumberSaves`` constrols the number of saves, namely, the number of output parameter samples
  * ``PTol`` constrols likelihood tolerance in loge 
  * ``ThreadStepsFactor`` constrols the steps that each core run before all cores communicate, i.e., 
    ThreadStep=ThreadStepsFactor * NumberParticles * NumberParameters
  * ``NewLevelIntervalFactor`` constrol the steps for creating a new likelihood level, i.e., 
    NewLevelInterval=NewLevelIntervalFactor * NumberCores * ThreadSteps
  * ``SaveIntervalFactor`` controls the steps that all cores run for saving a sample, i.e.,
    SaveInterval = SaveIntervalFactor * NumberCores * ThreadSteps
  * ``BacktrackingLength`` controls the backforward tracking length (lambda in Brewer's paper)
  * ``StrengthEqualPush`` controls strength to force equal push (beta in Brewer's paper)

The full arguments to the function ``dnest`` are 

.. code-block:: C

  double dnest(int argc, char** argv, DNestFptrSet *fptrset, int num_params, 
             double *param_range, int *prior_type, double *prior_info,
             char *sample_dir, char *optfile, DNestOptions *opts, void *args)
  
  /* ========================================================
  * arguments:
  * argc:        number of command-line options (mandatory)
  * argv:        command-line options  (mandatory)
  * fptrset:     function set pointers required  (mandatory)
  * num_params:  number of parameters (mandatory)
  * param_range: parameter ranges  (mandatory)
  * prior_type:  prior types  (optional)
  * prior_info:  prior informations  (optional)
  * sample_dir:  output directory for sampling  (mandatory)
  * optfile:     option files  (optional)
  * opts:        options struct. if optfile is empty, use this struct (optional)
  * args:        any other arguments transferred to cdnest.  (optional)
  *              useful when referring to external variables or calling external functions (optional)
  * ========================================================
  * optional arguments can be set to "NULL"
  *
  * return:
  * Bayesian evidence, double 
  */

Python
======
The statements in Python look like 

.. code-block:: Python

  from mpi4py import MPI
  import cydnest
  
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()

  class Model(object):
    """
    model passed to cdnest.
    """
    def __init__(self, num_params=1):
      """
      intialize the model
      """
      # number of parameters
      self.num_params = num_params 
  
      # parameter ranges, a list
      self.param_range = [[-5.0, 5.0]]*num_params
  
      # parameter prior type.
      # three types: Uniform, Gaussian, Log 
      self.prior_type = ["Uniform"]*num_params
  
      # parameter prior information. used when the prior is Gaussian
      # indicate the mean and standard deviation of the Gaussian prior
      # note:
      # 1) prior_info is a list consisting of num_params elements
      #    echo element is a list specifying the mean and standard deviation 
      # 2) the element is used only when the corresponding prior_type is gaussian,
      #    otherwise, it will be ignored.
      self.prior_info = [[0.0, 1.0]]*num_params
  
    def log_likelihood(self, coords):
      """
      calculate likelihood
      coords is a set of model parameters 
      """
      return -0.5*np.sum(coords**2) + self.num_params * (-0.5*np.log(2*np.pi))

  # create a model
  model = Model()
  
  # create a dnest sampler
  # max_num_save is the number of samples to generate
  # ptol is the likelihood tolerance in loge()
  sampler = cydnest.sampler(model, sample_dir="./", max_num_saves = 10000, ptol=0.1)
  #
  # The full argument lists look like:
  # sampler = cydnest.sampler(model, sample_dir="./", max_num_saves = 10000, ptol=0.1, 
  #               num_particles=1, thread_steps_factor = 10, 
  #               max_num_levels = 0, Lambda = 10, beta = 100
  #               new_level_interval_factor = 2, save_interval_factor = 2)
  #
  
  # run sampler
  logz = sampler.run()
  comm.Barrier()

  if rank == 0:
    
    print("Evidence:", logz)

    # load posterior sample 
    sample = np.loadtxt(sampler.get_sample_dir() +"posterior_sample" + sampler.get_sample_tag() + ".txt")

    # do postprocess, plot, show the properties of sampling 
    cydnest.postprocess(sampler.get_sample_dir(), sampler.get_sample_tag(), doplot=True)

In the folder **tests/** packaged along with the source code, a number of python examples are provided.