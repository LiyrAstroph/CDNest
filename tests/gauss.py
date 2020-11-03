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

class Model(object):

  def __init__(self, num_params=1):
    """
    intialize the model
    """
    self.num_params = num_params # number of parameters
    self.param_range = [[-5.0, 5.0]]
    self.prior_type = ["Uniform"]
    self.prior_info = [[0.0, 1.0]]

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


