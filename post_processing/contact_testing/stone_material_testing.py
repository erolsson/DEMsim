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

force_data = np.genfromtxt('../../results/contact_testing/stone_material_contact/contact_forces.dat',
                           delimiter=',')
h = (force_data[:, 0] - force_data[0, 0])/2
plt.plot(h, force_data[:, 1])

plt.show()
