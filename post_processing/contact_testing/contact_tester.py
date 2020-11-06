import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

force_data = np.genfromtxt('../../results/porous_electrode/contact_testing.dou',
                           delimiter=',')

plt.figure(1)
Ep = 139e9
vp = 0.3
E0 = Ep/(1-vp**2)/2
R0 = 0.5
h = force_data[:, 0]
plt.plot(h, force_data[:, 1], '-')
# F = 4./3*E0*R0**0.5*h**1.5
# plt.plot(h, F, '--')

plt.figure(2)
plt.plot(force_data[:, 1])

plt.figure(3)
plt.plot(force_data[:, 4], force_data[:, 3])

plt.figure(4)
plt.plot(force_data[:, 3])

plt.show()
