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
Simulation = namedtuple('Simulation', ['directory', 'line', 'fig', 'name'])
Figure = namedtuple('Figure', ['xlim', 'ylim', 'leg_handles'])
simulations = [Simulation(directory=base_directory + '/8-16mm/', line='b', fig=0,
                          name=r'Fraction 8-16 mm $\sigma_w=387.5$ MPa'),
               Simulation(directory=base_directory + '/8-16mm_weak/', line='--b', fig=0,
                          name=r'Fraction 8-16 mm $\sigma_w=200$ MPa'),
               Simulation(directory=base_directory + '/fuller1/', line='r', fig=0,
                          name=r'F{\"u}ller curve  $0 \% < 2$ mm $\sigma_w=387.5$ MPa'),
               Simulation(directory=base_directory + '/fuller_weak/', line='--r', fig=0,
                          name=r'F{\"u}ller curve  $0 \% < 2$ mm $\sigma_w=200$ MPa'),
               Simulation(directory=base_directory + '/8-16mm_continued/', line='b', fig=1,
                          name=r'Fraction 8-16 mm $\sigma_w=387.5$ MPa'),
               Simulation(directory=base_directory + '/8-16mm_continued_weak/', line='--b', fig=1,
                          name=r'Fraction 8-16 mm $\sigma_w=200$ MPa'),
               Simulation(directory=base_directory + '/8-16mm_one_layer/', line='b', fig=2,
                          name=r'Fraction 8-16 mm $\sigma_w=387.5$ MPa'),
               Simulation(directory=base_directory + '/8-16mm_one_layer_weak/', line='--b', fig=2,
                          name=r'Fraction 8-16 mm $\sigma_w=200$ MPa')]
figures = [Figure(xlim=125, ylim=400, leg_handles=[]),
           Figure(xlim=200, ylim=600, leg_handles=[]),
           Figure(xlim=125, ylim=300, leg_handles=[])]

for simulation in simulations:
    plt.figure(simulation.fig)
    with open(simulation.directory + '/fractured_particles.pkl') as pickle_file:
        fracture_data = pickle.load(pickle_file)
    leg_h = plt.plot(fracture_data[:, 0], fracture_data[:, 3], simulation.line, lw=2, label=simulation.name)
    figures[simulation.fig].leg_handles.append(leg_h[0])


for fig_number, figure in enumerate(figures):
    plt.figure(fig_number)
    if fig_number != 2:
        for layer_count in [25, 50, 75, 100]:
            plt.plot([layer_count, layer_count], [0, figure.ylim], '--k')
    plt.xlabel('Blows', fontsize=24)
    plt.ylabel('Particle cracks', fontsize=24)
    plt.xlim(0, figure.xlim)
    plt.ylim(0, figure.ylim)

plt.figure(0)
fig = plt.figure(0)
fig.set_size_inches(8., 10., forward=True)
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([0.125, 0.5, box.width, box.height/2])
box = ax.get_position()
legend = ax.legend(handles=figures[0].leg_handles, numpoints=1, loc='upper left', bbox_to_anchor=(-0.1, -0.2))
plt.gca().add_artist(legend)
plt.savefig('particle_cracks.png')

plt.figure(1)
plt.legend(handles=figures[1].leg_handles, numpoints=1, loc='best')
plt.savefig('continued.png')

plt.figure(2)
plt.legend(handles=figures[2].leg_handles, numpoints=1, loc='best')
plt.savefig('one_layer.png')
plt.show()
