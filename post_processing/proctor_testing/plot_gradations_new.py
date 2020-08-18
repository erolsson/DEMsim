import os
import pickle

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

base_directory = os.path.expanduser('~/DEMsim/post_processing/proctor_testing/gradation/')
simulations = ['8-16', 'Fuller']
colors = ['b', 'r']
labels = ['8-16 mm', 'Fuller curve']
bg_index = {'8-16 mm': {}, r'Fuller curve': {}}
for i, sim in enumerate(simulations):
    for material, line, in zip(['_start', '_strong', '_weak'], ['', '--', ':']):
        simulation = sim + material
        data = np.genfromtxt(base_directory + simulation + '.csv', delimiter=',')
        data = data[np.argsort(data[:, 0]), :]
        if material == '_start':
            plt.semilogx(data[:, 0], data[:, 1], line + colors[i], lw=2, label=labels[i])
        else:
            plt.semilogx(data[:, 0], data[:, 1], line + colors[i], lw=2)


for line, label in zip(['w', '-k', '--k', ':k'], [r'$\quad$', 'Starting',
                                                  r'$\sigma_w=387.5$ MPa', r'$\sigma_w=200$ MPa']):
    plt.plot([1, 2], [-1, -1], line, lw=2, label=label)
plt.xlim(1, 16)
plt.ylim(0, 100)
sieves = [1., 2., 4., 8., 11.2, 16.]

plt.xticks(sieves, [str(s) for s in sieves])
plt.xlabel('Size [mm]')
plt.ylabel(r'Weight passing [\%]')
ax = plt.gca()
box = ax.get_position()
# plt.text(0.05, 0.9, '(a)', horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
plt.legend()
plt.tight_layout()
plt.savefig('gradations_after_test_leg.png')

plt.show()
