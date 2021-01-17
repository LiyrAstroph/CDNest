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

  
  # users can define their own functions to generate 
  # the initial parameter values 
  # this is optinal. if not defined, cydnest will use the internal 
  # function.  
  def from_prior(self):
    """
    generate initial values of model parameters from priors
    """
    coords = np.zeros(self.num_params)
    for i in range(self.num_params):
      if self.prior_type[i] == "Uniform":
        coords[i] = np.random.uniform(self.param_range[i][0], self.param_range[i][1])
      elif self.prior_type[i] == "Gaussian":
        coords[i] = np.random.randn() * self.prior_info[i][1] + self.prior_info[0]
        wrap(coords[i], self.param_range[i][0], self.param_range[i][1])
      else: # LOG prior
        coords[i] = np.random.uniform(np.log(self.param_range[i][0]), np.log(self.param_range[i][1]))
        coords[i] = np.exp(coords[i])

    return coords

  # users can define their own functions to perturb 
  # parameter values for sampling 
  # this is optinal. if not defined, cydnest will use the internal 
  # function.                                  
  def perturb(self, coords):
    """
    perturb the parameters
    """
    i = np.random.randint(self.num_params)
    
    LogH = 0.0   # prior ratio: ln(prior(new)/prior(old)) = ln(prior(new)) - ln(prior(old))
    width = (self.param_range[i][1]-self.param_range[i][0])
    if self.prior_type[i] == "Uniform":
      coords[i] += width*randh()
      coords[i] = wrap(coords[i], self.param_range[i][0], self.param_range[i][1])
    elif self.prior_type[i] == "Gaussian":  
      LogH -= ( -0.5* (coords[i] - self.prior_info[i][0])**2/self.prior_info[i][1]**2  ) # ln(Gaussian)
      coords[i] += width*randh()
      coords[i] = wrap(coords[i], self.param_range[i][0], self.param_range[i][1])
      LogH += ( -0.5* (coords[i] - self.prior_info[i][0])**2/self.prior_info[i][1]**2  )
    elif self.prior_type[i] == "Log":
      LogH -= ( -np.log(coords[i]) )   # ln(1/x) = -ln(x)
      coords[i] += width*randh()
      coords[i] = wrap(coords[i], self.param_range[i][0], self.param_range[i][1])
      LogH += ( -np.log(coords[i]) )
    return LogH 
  
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
# ptol is the likelihood tolerance in loge()
sampler = cydnest.sampler(model, sample_dir="./", max_num_saves = 10000, ptol=0.1)
#
# The full argument lists look like:
# sampler = cydnest.sampler(model, sample_dir="./", max_num_saves = 10000, ptol=0.1, 
#               num_particles=1, thread_steps_factor = 10, 
#               max_num_levels = 0, lam = 10, beta = 100
#               new_level_interval_factor = 2, save_interval_factor = 2)
#


# run sampler
logz = sampler.run()
comm.Barrier()

# ouput evidence
if rank == 0:
  print("Evidence:", logz)

  psample = np.loadtxt(sampler.get_sample_dir() +"/posterior_sample" + sampler.get_sample_tag() + ".txt")
  psample_info = np.loadtxt(sampler.get_sample_dir() +"/posterior_sample_info" + sampler.get_sample_tag() + ".txt")

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

  # do postprocess, plot, show the properties of sampling 
  cydnest.postprocess(sampler.get_sample_dir(), sampler.get_sample_tag(), temperature=1.0, doplot=True)