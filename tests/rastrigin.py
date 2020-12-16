#
# sample from a Rastrigin test function
#

from mpi4py import MPI
import numpy as np
import cydnest
import matplotlib.pyplot as plt
from matplotlib import cm

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
    # number of parameters
    self.num_params = num_params 

    # parameter ranges, a list
    self.param_range = [[-5.12, 5.12]]*num_params

    # parameter prior type.
    # three types: Uniform, Gaussian, Log 
    self.prior_type = ["Uniform"]*num_params

    # parameter prior information. used when the prior is Gaussian
    # indicate the mean and standard deviation of the Gaussian prior
    self.prior_info = [[0.0, 1.0]]*num_params


  def from_prior(self):
    """
    generate initial values of model parameters from priors
    """
    return np.random.uniform(-5.12, 5.12,size=(self.num_params,))
                                    
  def perturb(self, coords):
    """
    perturb the parameters
    """
    i = np.random.randint(self.num_params)
    coords[i] += 10.0*randh()
    coords[i] = wrap(coords[i], -5.12, 5.12)
    return 0.0  

  def log_likelihood(self, coords):
    """
    calculate likelihood
    """
    return -(10.0*2 + (coords[0]**2 - 10*np.cos(2.0*np.pi*coords[0])) + (coords[1]**2 - 10*np.cos(2.0*np.pi*coords[1])) )

# create a model
model = Model(num_params=2)

# create a dnest sampler
# max_num_save is the number of samples to generate
# max_num_levels is the number of levels 
sampler = cydnest.sampler(model, sample_dir="./", max_num_saves = 10000, ptol=0.1)

# run sampler
logz = sampler.run()
comm.Barrier()

# ouput evidence
if rank == 0:
  print("Evidence:", logz)

  psample = np.loadtxt(sampler.get_sample_dir() +"posterior_sample" + sampler.get_sample_tag() + ".txt")
  psample_info = np.loadtxt(sampler.get_sample_dir() +"posterior_sample_info" + sampler.get_sample_tag() + ".txt")

  fig = plt.figure(figsize=(15, 12))
  ax = fig.add_subplot(111, projection='3d')
  

  X = np.arange(-1.5, 1.5, 0.01)
  Y = np.arange(-1.5, 1.5, 0.01)
  X, Y = np.meshgrid(X, Y)
  Z = -(10.0*2 + (X**2 - 10*np.cos(2.0*np.pi*X)) + (Y**2 - 10*np.cos(2.0*np.pi*Y)) )
  ax.plot_surface(X, Y, Z, cmap=cm.ocean, rstride=2, cstride=2, linewidth=0, antialiased=False, zorder=0)

  idx = np.where((np.abs(psample[:, 0]) <1.4) & (np.abs(psample[:, 1]) <1.4))
  ax.plot(psample[idx[0], 0], psample[idx[0], 1], psample_info[idx[0]], ls='none', marker='+', zorder=10)
  ax.set_xlim(-1.5, 1.5)
  ax.set_ylim(-1.5, 1.5)
  ax.set_xlabel(r'$\theta_1$')
  ax.set_ylabel(r'$\theta_2$')
  ax.set_zlabel(r'$\log L$')
  fig.savefig("fig_rastrigin.jpg", bbox_inches='tight')
  plt.show()

  # do postprocess, plot 
  cydnest.postprocess(sampler.get_sample_dir(), sampler.get_sample_tag(), doplot=True)
  

  

