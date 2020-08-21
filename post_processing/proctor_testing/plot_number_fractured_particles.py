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

base_directory = os.path.expanduser('~/DEMsim/post_processing/proctor_testing/particle_cracks/')
Simulation = namedtuple('Simulation', ['file', 'line', 'fig', 'name'])
Figure = namedtuple('Figure', ['xlim', 'ylim', 'leg_handles', 'label', 'name'])
simulations = [Simulation(file=base_directory + '/8-16_strong.csv', line='--b', fig=0,
                          name=r'8-16 mm $\sigma_W=387.5$ MPa'),
               Simulation(file=base_directory + '/8-16_weak.csv', line=':b', fig=0,
                          name=r'8-16 mm $\sigma_W=200$ MPa'),
               Simulation(file=base_directory + '/fuller_strong.csv', line='--r', fig=0,
                          name=r'Fuller curve  $\sigma_W=387.5$ MPa'),
               Simulation(file=base_directory + '/fuller_weak.csv', line=':r', fig=0,
                          name=r'Fuller curve $\sigma_W=200$ MPa')]
figures = [Figure(xlim=125, ylim=200, leg_handles=[], label='(a)', name='proctor_leg.png')]
           # Figure(xlim=200, ylim=300, leg_handles=[], label='(a)', name='proctor_continued.png'),
           # Figure(xlim=125, ylim=150, leg_handles=[], label='(b)', name='proctor_one_layer.png')]

for simulation in simulations:
    plt.figure(simulation.fig)
    fracture_data = np.genfromtxt(simulation.file, delimiter=",")
    leg_h = plt.plot(fracture_data[:, 0]/120*125, fracture_data[:, 1], simulation.line, lw=2, label=simulation.name)
    figures[simulation.fig].leg_handles.append(leg_h[0])


for fig_number, figure in enumerate(figures):
    plt.figure(fig_number)
    if fig_number in range(0, 2):
        for layer_count in [25, 50, 75, 100]:
            plt.plot([layer_count, layer_count], [0, figure.ylim], '--k')
    plt.xlabel('Blows', fontsize=24)
    plt.ylabel('Particle cracks', fontsize=24)
    plt.xlim(0, figure.xlim)
    plt.ylim(0, figure.ylim)
    ax = plt.subplot(111)
    # plt.text(0.05, 0.9, figure.label, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    plt.legend(loc='upper left', bbox_to_anchor=(0., 0.89), framealpha=1.)
    plt.tight_layout()
    plt.savefig(figure.name)

plt.show()
