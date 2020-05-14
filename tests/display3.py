import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D


sample = np.loadtxt("sample3.txt", skiprows=1)
sample_info = np.loadtxt("sample_info3.txt", skiprows=1)

fig = plt.figure()

ax = fig.gca(projection = '3d')

ax.plot(sample[:, 0], sample[:, 1], np.exp(sample_info[:, 1]), ls='none', marker='.', markersize=1.0, color='r')

ax.set_xlim(-6, 6)
ax.set_ylim(-6, 6)
ax.set_zlim(0, 5)

ax.set_zlabel(r'Likelihood')
ax.set_xlabel(r'x')
ax.set_ylabel(r'y')

fig.savefig("fig3.jpg", bbox_inches='tight')

plt.show()
