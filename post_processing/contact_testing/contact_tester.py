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

force_data = np.genfromtxt('../../results/viscoelastic/contact_force_control',
                           delimiter=',')

plt.figure(1)
Ep = 139e9
vp = 0.25
E0 = Ep/(1-vp**2)/2
R0 = 0.03/2
h = force_data[:, 0]
plt.plot(h, force_data[:, 1], '-')
F = 4./3*E0*R0**0.5*h**1.5
plt.plot(h, F, '--')

plt.figure(2)
plt.plot(force_data[:, 1])

plt.show()
