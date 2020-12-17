#
# an example borrowed from the DNest4 package by Brendon J. Brewer, with minor modifications
# 
#

from mpi4py import MPI
import numpy as np
import scipy.special as sp
import cydnest

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def analytic_log_Z(num_params):
  """
  return evidence for multi-dimensional Gaussian
  """
  width = 10.0
  return (
    num_params * np.log(sp.erf(0.5*width/np.sqrt(2))) - num_params * np.log(width)
  )

class Model(object):
  """
  model input to cdnest
  """
  def __init__(self, num_params=5):
    """
    intialize the model.
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
    self.prior_info = [[0.0, 1.0]]*num_params

  def log_likelihood(self, coords):
    """
    calculate likelihood
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

# ouput evidence
if rank == 0:
  # print evidence
  print("Caclulated Evidence:", logz, ", Real Evidence:", analytic_log_Z(model.num_params))
  
  # do postprocess, plot, show the properties of sampling 
  cydnest.postprocess(sampler.get_sample_dir(), sampler.get_sample_tag(), doplot=True)
