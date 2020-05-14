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

  def __init__(self, num_params=1, width=10.0):
    self.num_params = num_params
    self.width = width
    self.options_file = "OPTIONS4" # option file for cydnest  

  def from_prior(self):
    return np.random.uniform(-0.5*self.width, 0.5*self.width,
                                  size=(self.num_params,))  
  def perturb(self, coords):
    i = np.random.randint(self.num_params)
    coords[i] += self.width*randh()
    coords[i] = wrap(coords[i], -0.5*self.width, 0.5*self.width)
    return 0.0  

  def log_likelihood(self, coords, const=-0.5*np.log(2*np.pi)):
    return -0.5*np.sum(coords**2) + self.num_params * const

# create a model
model = Model()

# create a dnest sampler
sample = cydnest.sampler(model)

# run sampler
logz = sample.run()
comm.Barrier()

# ouput evidence
if rank == 0:
  print(logz)

  sample = np.loadtxt("posterior_sample4.txt")
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


