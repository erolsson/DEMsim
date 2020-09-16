import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

particle_data = np.genfromtxt(os.path.expanduser('~/DEMsim/results/deformable_surfaces/particle_1.dou'), delimiter=',')
time = particle_data[:, 0]
x = particle_data[:, 2]
ry = particle_data[:, 6]
plt.figure(0)
plt.plot(time, x)

plt.figure(1)
plt.plot(time, ry)
plt.show()
