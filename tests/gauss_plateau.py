#
# 2d clipped Gaussian distribution
# 
# L(x, y) = 1/(2pi) exp( - (x^2+y^2)/2) )  for x^2+y^2 >= 4
# L(x, y) = 1/(2pi) exp( - 2 )             for x^2+y^2 < 4
# 
# for uniform prior of x, y in (-5, 5)
# the evidence is  Z = 3*exp(-2)/100
#                  log(Z) = -2 + log(3/100) = -5.5
#

from mpi4py import MPI
import numpy as np
import cydnest
import matplotlib.pyplot as plt
from matplotlib import cm

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

class Model(object):
  """
  model passed to cdnest.
  """
  def __init__(self, num_params=2):
    """
    intialize the model
    """
    # number of parameters
    self.num_params = num_params 

    # parameter ranges, a list
    self.param_range = [[-5.0, 5.0], [-5.0, 5.0]]

    # parameter prior type.
    # three types: Uniform, Gaussian, Log 
    self.prior_type = ["Uniform", "Uniform"]
    
    # if using Gaussian prior
    # parameter prior information. used when the prior is Gaussian
    # indicate the mean and standard deviation of the Gaussian prior
    #self.prior_info = [[0.0, 1.0]]*num_params

  def log_likelihood(self, coords):
    """
    calculate likelihood
    """
    sig = 1.0
    Const = -0.5*( 2.0**2/sig**2 ) - np.log(2.0*np.pi*sig**2)
    logL = -0.5*( (coords[0]**2+coords[1]**2)/sig**2) - np.log(2.0*np.pi*sig**2)
    return np.min((logL, Const))
     

# create a model
model = Model()

# create a dnest sampler
# max_num_save is the number of samples to generate
# ptol is the likelihood tolerance in loge()
sampler = cydnest.sampler(model, sample_dir="./", max_num_saves = 20000, ptol=0.1)
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
  sample = np.loadtxt(sampler.get_sample_dir() +"/posterior_sample" + sampler.get_sample_tag() + ".txt")
  sample_info = np.loadtxt(sampler.get_sample_dir() +"/posterior_sample_info" + sampler.get_sample_tag() + ".txt")

  fig = plt.figure(figsize=(15, 12))
  ax = fig.add_subplot(111, projection='3d')
  
  X = np.arange(-5, 5, 0.05)
  Y = np.arange(-5, 5, 0.05)
  X, Y = np.meshgrid(X, Y)
  Const = np.zeros(X.shape)
  sig = 1.0
  Const[:, :] = -0.5*( 2.0**2/sig**2) - np.log(2.0*np.pi*sig**2)
  Z = np.where(X**2+Y**2<2.0**2, Const,-0.5*( (X**2+Y**2)/sig**2) - np.log(2.0*np.pi*sig**2))
  Z = np.exp(Z)
  ax.plot_surface(X, Y, Z, cmap=cm.ocean, rstride=2, cstride=2, linewidth=0, antialiased=False, zorder=0)
  ax.plot(sample[:, 0], sample[:, 1], np.exp(sample_info[:]), ls='none', marker='+', zorder=10)
  ax.set_xlabel(r'$x_1$')
  ax.set_ylabel(r'$x_2$')
  ax.set_zlabel(r'$\log L$')
  ax.set_xlim((-5, 5))
  ax.set_ylim((-5, 5))
  fig.savefig("fig_gauss_plateau.jpg", bbox_inches='tight', dpi=200)
  plt.show()

  # do postprocess, plot, show the properties of sampling 
  cydnest.postprocess(sampler.get_sample_dir(), sampler.get_sample_tag(), temperature=1.0, doplot=True)