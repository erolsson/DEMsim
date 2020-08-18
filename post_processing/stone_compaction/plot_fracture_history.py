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

base_directory = os.path.expanduser('~/DEMsim/results/stone_compaction/')
Simulation = namedtuple('Simulation', ['directory', 'line', 'fig', 'name'])
Figure = namedtuple('Figure', ['xlim', 'ylim', 'leg_handles'])
simulations = [Simulation(directory=base_directory + '/8-16mm/', line='b', fig=0,
                          name=r'Fraction 8-16 mm $\sigma_w=387.5$ MPa'),
               Simulation(directory=base_directory + '/8-16mm_weak/', line='--b', fig=0,
                          name=r'Fraction 8-16 mm $\sigma_w=200$ MPa'),
               Simulation(directory=base_directory + '/fuller/', line='r', fig=0,
                          name=r'F{\"u}ller curve  $0 \% < 2$ mm $\sigma_w=387.5$ MPa'),
               Simulation(directory=base_directory + '/fuller_weak/', line='--r', fig=0,
                          name=r'F{\"u}ller curve  $0 \% < 2$ mm $\sigma_w=200$ MPa')]

figures = [Figure(xlim=100, ylim=6000, leg_handles=[]),
           Figure(xlim=100, ylim=20, leg_handles=[]),
           Figure(xlim=20, ylim=10e3, leg_handles=[])]

for simulation in simulations:
    with open(simulation.directory + "/crack_history.pkl", 'r') as crack_pickle:
        data = pickle.load(crack_pickle)
    plt.figure(0)
    leg_h = plt.plot((data[0, 1] - data[:, 1])*1000, data[:, 3], simulation.line, lw=2, label=simulation.name)
    figures[0].leg_handles.append(leg_h[0])

    plt.figure(1)
    leg_h = plt.semilogx(data[:, 2]/np.pi/(2*25.4)**2, data[:, 3], simulation.line, lw=2, label=simulation.name)
    figures[1].leg_handles.append(leg_h[0])

    plt.figure(2)
    leg_h = plt.plot((data[0, 1] - data[:, 1])*1000, data[:, 2], simulation.line, lw=2, label=simulation.name)
    figures[2].leg_handles.append(leg_h[0])

for fig_nr in [0, 1, 2]:
    fig = plt.figure(fig_nr)
    fig.set_size_inches(8., 10., forward=True)
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([0.125, 0.5, box.width, box.height/2])
    legend = ax.legend(handles=figures[0].leg_handles, numpoints=1, loc='upper left', bbox_to_anchor=(-0.1, -0.2))
    plt.gca().add_artist(legend)
    plt.xlim(0.1, figures[fig_nr].xlim)
    plt.ylim(0, figures[fig_nr].ylim)

plt.figure(0)
plt.xlabel('Plate travel [mm]')
plt.ylabel('Particle cracks', fontsize=24)
plt.savefig('particle_cracks_displacement.png')

plt.figure(1)
plt.xlabel('Compaction pressure [MPa]')
plt.ylabel('Particle cracks', fontsize=24)
plt.savefig('particle_cracks_force.png')

plt.figure(2)
plt.xlabel('Plate travel [mm]')
plt.ylabel('Compaction force [kN]')
plt.savefig('force_displacement.png')

plt.show()
