#
# the diamond ring problem
# see Buchner, J. 2023, arXiv:2101.09675
#
from mpi4py import MPI
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
import corner
import cydnest

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
 
def analytic_log_Z():
  """
  return analytic evidence
  """
  r1 = 1.0
  r2 = r1/40.0
  w1 = 0.4
  w2 = w1/40.0

  Z = np.sqrt(2*np.pi)/4*w1*np.exp(-r1**2/2/w1**2)+np.pi/4*r1*sp.erfc(-r1/np.sqrt(2)/w1) \
     +100*(np.sqrt(2*np.pi)/4*w2*np.exp(-r2**2/2/w2**2)+np.pi/4*r2*sp.erfc(-r2/np.sqrt(2)/w2))
  
  return np.log(Z) + np.log(1.0e-11)

def logL(x, y):
  r1 = 1.0e-11
  r2 = r1/40.0
  w1 = 0.4*r1 
  w2 = w1/40.0
  d1 = np.sqrt(x**2+y**2)
  d2 = np.sqrt((x+r1)**2+y**2)
  
  n1_log = -0.5*(d1-r1)**2/w1**2 - 0.5*np.log(2.0*np.pi*w1**2)
  n2_log = -0.5*(d2-r2)**2/w1**2 - 0.5*np.log(2.0*np.pi*w2**2)
  return np.logaddexp(n1_log, n2_log)
  
class Model(object):
  """
  model input to cdnest
  """
  def __init__(self, num_params=2):
    """
    intialize the model.
    """
    # number of parameters
    self.num_params = num_params 

    # parameter ranges, a list
    self.param_range = [[-1.0, 1.0]]*num_params

    # parameter prior type.
    # three types: Uniform, Gaussian, Log
    self.prior_type = ["Uniform"]*num_params
    
    # if Guassian,
    # parameter prior information. used when the prior is Gaussian
    # indicate the mean and standard deviation of the Gaussian prior
    # self.prior_info = [[0.0, 1.0]]*num_params

  def log_likelihood(self, coords):
    """
    calculate likelihood
    """
    r1 = 1.0e-11
    r2 = r1/40.0
    w1 = 0.4*r1 
    w2 = w1/40.0
    d1 = np.sqrt(coords[0]**2+coords[1]**2)
    d2 = np.sqrt((coords[0]+r1)**2+coords[1]**2)
      
    n1_log = -0.5*(d1-r1)**2/w1**2 - 0.5*np.log(2.0*np.pi*w1**2)
    ex = -0.5*(d2-r2)**2/w2**2 + 0.5*(d1-r1)**2/w1**2
    
    if ex >= 100.0:
      return n1_log + np.log(100.0) + ex + np.log(w1/w2)
    if ex <= -100.0:
      n2rn1 = w1/w2 * np.exp(ex)
      return n1_log + (100*n2rn1)
    else:
      n2rn1 = w1/w2 * np.exp(ex)
      return n1_log + np.log(1.0+100.0*n2rn1) 
      

# create a model
model = Model()

# create a dnest sampler
# max_num_save is the number of samples to generate
# ptol is the likelihood tolerance in loge()
sampler = cydnest.sampler(model, sample_dir="./", max_num_saves = 20000, ptol=1.0e-3, num_particles=10, limits_on=True,
                          compression = np.exp(0.1))
#
# The full argument lists look like:
# sampler = cydnest.sampler(model, sample_dir="./", max_num_saves = 10000, ptol=0.1, 
#               num_particles=1, thread_steps_factor = 10, 
#               max_num_levels = 0, lam = 10, beta = 100
#               new_level_interval_factor = 2, save_interval_factor = 2, 
#               limits_on = False, compression = None)
#

# run sampler
logz = sampler.run()
comm.Barrier()

# ouput evidence
if rank == 0:
  # print evidence
  print("Caclulated Evidence:", logz, ", Real Evidence:", analytic_log_Z())
  
  # do postprocess, plot, show the properties of sampling 
  cydnest.postprocess(sampler.get_sample_dir(), sampler.get_sample_tag(), temperature=1.0, doplot=True)

  # plot posterior sample 
  sample = np.loadtxt("posterior_sample.txt")

  x = np.linspace(-3e-11, 3e-11, 500)
  mx, my = np.meshgrid(x, x)
  mz = np.exp(logL(mx, my))
  
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.contourf(mx, my, mz, levels=np.exp(-0.5*np.arange(30.0, 0.0,-0.1)**2)*np.max(mz), cmap=cm.ocean, \
              norm=colors.LogNorm(vmin=mz.min(), vmax=mz.max()))
  ax.plot(sample[:, 0], sample[:, 1], ls='none', marker='o', markersize=0.5)
  ax.set_aspect('equal', adjustable='box')

  axins = ax.inset_axes((0.05, 0.65, 0.3, 0.3))
  axins.contourf(mx, my, mz, levels=np.exp(-0.5*np.arange(30.0, 0.0,-0.1)**2)*np.max(mz), cmap=cm.ocean, \
              norm=colors.LogNorm(vmin=mz.min(), vmax=mz.max()))
  axins.plot(sample[:, 0], sample[:, 1], ls='none', marker='o', markersize=0.5)
  
  axins.minorticks_on()
  
  x1, x2, y1, y2 = -1.1e-11, -0.9e-11, -0.1e-11, 0.1e-11
  axins.set_xlim(x1, x2)
  axins.set_ylim(y1, y2)
  axins.set_xticklabels([])
  axins.set_yticklabels([])
  ax.indicate_inset_zoom(axins, edgecolor="black")
  ax.set_xlabel('x')
  ax.set_ylabel('y')
  fig.savefig("diamond_ring.jpg", bbox_inches='tight', dpi=200)
  plt.show()
  

