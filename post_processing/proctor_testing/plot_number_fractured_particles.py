from collections import namedtuple
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

base_directory = os.path.expanduser('~/DEMsim/results/proctor_test')
Simulation = namedtuple('Simulation', ['directory', 'line', 'name'])
simulations = [Simulation(directory=base_directory + '/8-16mm/', line='b',
                          name=r'Fraction 8-16 mm $\sigma_w=387.5$ MPa'),
               Simulation(directory=base_directory + '/8-16mm_weak/', line='--b',
                          name=r'Fraction 8-16 mm $\sigma_w=200$ MPa'),
               Simulation(directory=base_directory + '/fuller1/', line='r',
                          name=r'F{\"u}ller curve  $0 \% < 2$ mm $\sigma_w=387.5$ MPa'),
               Simulation(directory=base_directory + '/fuller_weak/', line='--r',
                          name=r'F{\"u}ller curve  $0 \% < 2$ mm $\sigma_w=200$ MPa')]
fig = plt.figure(0)
legend_handles = []
for simulation in simulations:
    with open(simulation.directory + '/fractured_particles.pkl') as pickle_file:
        fracture_data = pickle.load(pickle_file)
    leg_h = plt.plot(fracture_data[:, 0], fracture_data[:, 3], simulation.line, lw=2, label=simulation.name)
    legend_handles.append(leg_h[0])

for layer_count in [25, 50, 75, 100]:
    plt.plot([layer_count, layer_count], [0, 400], '--k')
plt.xlabel('Blows', fontsize=24)
plt.ylabel('Particle cracks', fontsize=24)
fig.set_size_inches(8., 10., forward=True)
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([0.125, 0.5, box.width, box.height/2])
box = ax.get_position()
legend = ax.legend(handles=legend_handles, numpoints=1, loc='upper left', bbox_to_anchor=(-0.1, -0.2))
plt.gca().add_artist(legend)
plt.xlim(0, 125)
# plt.tight_layout()
plt.savefig('particle_cracks.png')
plt.show()
