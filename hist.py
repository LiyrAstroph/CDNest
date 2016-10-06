import numpy as np
import matplotlib.pyplot as plt
import corner

data = np.loadtxt("posterior_sample.txt")

figure = corner.corner(data)

plt.show()

