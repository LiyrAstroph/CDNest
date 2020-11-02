#
# an example borrowed from the DNest4 package by Brendon J. Brewer, with minor modifications 
#
# sample from a Gaussian distribution
#

from mpi4py import MPI
import numpy as np
import cydnest
import matplotlib.pyplot as plt

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

  def __init__(self, num_params=1):
    """
    intialize the model
    """
    self.num_params = num_params # number of parameters
    #self.options_file = "OPTIONS" # optional, if not set, use the default options

  def from_prior(self):
    """
    generate initial values of model parameters from priors
    """
    return np.random.uniform(-5.0, 5.0,size=(self.num_params,))
                                    
  def perturb(self, coords):
    """
    perturb the parameters
    """
    i = np.random.randint(self.num_params)
    coords[i] += 10.0*randh()
    coords[i] = wrap(coords[i], -5.0, 5.0)
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
sample = cydnest.sampler(model, sample_dir="./", max_num_saves = 10000, ptol=0.1, thread_steps_factor=100, 
                                new_level_interval_factor = 5, save_interval_factor = 5)

# run sampler
logz = sample.run()
comm.Barrier()

# ouput evidence
if rank == 0:
  print("Evidence:", logz)

  sample = np.loadtxt(sample.get_sample_dir() +"posterior_sample" + sample.get_sample_tag() + ".txt")
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.hist(sample, bins=20, density=True, label='CDNest sample')
  x = np.linspace(-5, 5, 100)
  y = 1.0/np.sqrt(2.0*np.pi) * np.exp(-0.5*x**2)
  ax.plot(x, y, color='red', label='Gaussian')
  ax.set_xlabel('x')
  ax.set_ylabel('probability density')
  ax.legend(frameon=False)
  fig.savefig("fig_gau.jpg")
  plt.show()


