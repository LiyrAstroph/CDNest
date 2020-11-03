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
    self.num_params = num_params  # number of parameters
    self.param_range = [[-5.0, 5.0]]*num_params
    self.prior_type = ["Uniform"]*num_params
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
# max_num_levels is the number of levels 
sample = cydnest.sampler(model, sample_dir="./", max_num_saves = 10000, ptol=0.1, thread_steps_factor=50, 
                                new_level_interval_factor = 5, save_interval_factor = 5)

# run sampler
logz = sample.run()
comm.Barrier()

# ouput evidence
if rank == 0:
  print("Caclulated Evidence:", logz, ", Real Evidence:", analytic_log_Z(model.num_params))
