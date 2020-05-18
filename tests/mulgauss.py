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
  assert b > a
  return (x - a)%(b - a) + a


class Model(object):
  def __init__(self, num_params=5, width=10.0):
    """
    intialize the model.
    """
    self.num_params = num_params
    self.width = width
    self.options_file = "OPTIONS5"  # options file for dnest

  def analytic_log_Z(self):
    return (
      self.num_params * np.log(sp.erf(0.5*self.width/np.sqrt(2))) -
      self.num_params * np.log(self.width)
    )

  def from_prior(self):
    """
    generate initial values of model parameters from priors
    """
    return np.random.uniform(-0.5*self.width, 0.5*self.width,
                                    size=(self.num_params,))

  def perturb(self, coords):
    """
    perturb the parameters
    """
    i = np.random.randint(self.num_params)
    coords[i] += self.width*randh()
    # Note: use the return value of wrap, unlike in C++
    coords[i] = wrap(coords[i], -0.5*self.width, 0.5*self.width)
    
    return 0.0

  def log_likelihood(self, coords, const=-0.5*np.log(2*np.pi)):
    """
    calculate likelihood
    """
    return -0.5*np.sum(coords**2) + self.num_params * const

# create a model
model = Model()

# create a dnest sampler
sample = cydnest.sampler(model, sample_dir="./", sample_tag="5")

# run sampler
logz = sample.run()
comm.Barrier()

# ouput evidence
if rank == 0:
  print("Caclulated Evidence:", logz, ", Real Evidence:", model.analytic_log_Z())
