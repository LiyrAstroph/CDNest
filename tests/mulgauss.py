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

def randh(N=1):
  """
  generate from the heavy-tailed distribution.
  """
  if N==1:
      return 10.0**(1.5 - 3*np.abs(np.random.randn()/np.sqrt(-np.log(np.random.rand()))))*np.random.randn()
  return 10.0**(1.5 - 3*np.abs(np.random.randn(N)/np.sqrt(-np.log(np.random.rand(N)))))*np.random.randn(N)


def wrap(x, a, b):
  """
  wrap x into a range [a, b]
  """
  assert b > a
  return (x - a)%(b - a) + a


def analytic_log_Z(num_params):
  """
  return evidence for multi-dimensional Gaussian
  """
  width = 10.0
  return (
    num_params * np.log(sp.erf(0.5*width/np.sqrt(2))) - num_params * np.log(width)
  )

class Model(object):
  def __init__(self, num_params=5):
    """
    intialize the model.
    """
    self.num_params = num_params  # number of parameters

  def from_prior(self):
    """
    generate initial values of model parameters from priors
    """
    return np.random.uniform(-5.0, 5.0, size=(self.num_params,))

  def perturb(self, coords):
    """
    perturb the parameters
    """
    width = 10.0
    i = np.random.randint(self.num_params)
    coords[i] += width*randh()
    coords[i] = wrap(coords[i], -0.5*width, 0.5*width)
    return 0.0

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
sample = cydnest.sampler(model, sample_dir="./", max_num_saves = 10000, max_num_levels=20)

# run sampler
logz = sample.run()
comm.Barrier()

# ouput evidence
if rank == 0:
  print("Caclulated Evidence:", logz, ", Real Evidence:", analytic_log_Z(model.num_params))
