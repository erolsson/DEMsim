import numpy as np

import matplotlib.pyplot as plt
import matplotlib

data = np.genfromtxt("contact_forces.dat", delimiter=',')
plt.plot(data[:, 1] + 1e-3, data[:, 0])
r = 0.00999276 + 1e-3
E = 1.7e9
v = 0.49
E0 = E/(1-v**2)
F = 4*E0*np.sqrt(r)/3*(data[:, 1] + 1e-3)**1.5
plt.plot(data[:, 1] + 1e-3, F, '--')
plt.show()
