#
# sample from a Rastrigin test function
# this is to illustrate how to use accept_action in CDNest to avoid repeat calculations.
#
# A 2D Rastrigin function looks
# 
# logL=-(10.0*2 + (coords[0]**2 - 10*np.cos(2.0*np.pi*coords[0])) + (coords[1]**2 - 10*np.cos(2.0*np.pi*coords[1])) ) 
#
# Every perturb, only one parameter is updated, so that the terms related to the rest parameters 
# do not need to recalculate, just use the values in the previous step.
#
# In this example, we use an array to record values of the term "(coords[0]**2 - 10*np.cos(2.0*np.pi*coords[0]))"
# in every accepted perturb.
#

from mpi4py import MPI
import numpy as np
import cydnest
import matplotlib.pyplot as plt
from matplotlib import cm

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def randh():
  """
  generate from the heavy-tailed distribution.
  """
  return 10.0**(1.5 - 3*np.abs(np.random.randn()/np.sqrt(-np.log(np.random.rand()))))*np.random.randn()

def wrap(x, a, b):
  assert b > a
  return (x - a)%(b - a) + a

class Model(object):

  def __init__(self, num_params=1, num_particles=1):
    """
    intialize the model
    """
    # number of particles each core holds
    self.num_particles = num_particles

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
    
    # which parameter being perturbed 
    # which particle being perturbed
    self.which_param_update = 0
    self.which_particle_update = 0

    # perturbed values and accepted values for all particles
    self.value_perturb = [0.0]*self.num_particles
    self.value_accept = [0.0]*self.num_particles

  def accept_action(self):
    """
    action taken when a perturb is accepted
    record the accepted values from the perturbed values
    """

    # note "which_particle_update" is updated and "which_param_update" is updated
    if self.which_param_update < 1:
      self.value_accept[self.which_particle_update] = self.value_perturb[self.which_particle_update]
  
  def kill_action(self, i, i_copy):
    """
    cdnest kill a particle when it is not updated for a long time.
    action taken when a particle is killed: i particle is killed,
    copy i_copy particle's values to i particle's values
    this function is needed, since we record some accepted values 
    """
    self.value_accept[i] = self.value_accept[i_copy]
    return
    
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
      elif self.prior_type[i] == "Log": # LOG prior
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
    
    # record which parameter is updated
    self.which_param_update = i

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
  
  def log_likelihood_initial(self, coords):
    """
    calculate likelihood at initial start
    """    
    self.which_particle_update = cydnest.get_which_particle_update()
    self.value_accept[self.which_particle_update] = coords[0]**2 - 10*np.cos(2.0*np.pi*coords[0])
    value =  self.value_accept[self.which_particle_update]
    return -(10.0*2 + (value) + (coords[1]**2 - 10*np.cos(2.0*np.pi*coords[1])) )

  def log_likelihood(self, coords):
    """
    calculate likelihood
    """
    # get which particle is being updated, and save it to self model

    self.which_particle_update = cydnest.get_which_particle_update()
    
    value = 0.0
    if self.which_param_update < 1: # when 0-th parameter update,  recalculate
      self.value_perturb[self.which_particle_update] = coords[0]**2 - 10*np.cos(2.0*np.pi*coords[0])
      value =  self.value_perturb[self.which_particle_update]
    else:  # otherwise, use the accepted value
      value =  self.value_accept[self.which_particle_update]

    return -(10.0*2 + (value) + (coords[1]**2 - 10*np.cos(2.0*np.pi*coords[1])) )

# create a model
model = Model(num_params=2, num_particles=2)

# create a dnest sampler
# max_num_save is the number of samples to generate
# max_num_levels is the number of levels 
# ptol is the likelihood tolerance in loge()
sampler = cydnest.sampler(model, sample_dir="./", max_num_saves = 10000, ptol=0.1, num_particles=model.num_particles)
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

  # do postprocess, plot, show the properties of sampling 
  cydnest.postprocess(sampler.get_sample_dir(), sampler.get_sample_tag(), temperature=1.0, doplot=True)