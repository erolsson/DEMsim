from collections import namedtuple
import os
import pickle

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
Figure = namedtuple('Figure', ['xlim', 'ylim', 'leg_handles', 'label', 'name'])
simulations = [Simulation(directory=base_directory + '/8-16mm_continued/', line='b', fig=0,
                          name=r'8-16 mm $\sigma_w=387.5$ MPa'),
               Simulation(directory=base_directory + '/8-16mm_continued_weak/', line='--b', fig=0,
                          name=r'8-16 mm $\sigma_w=200$ MPa'),
               Simulation(directory=base_directory + '/fuller/', line='r', fig=1,
                          name=r'F{\"u}ller curve  $\sigma_w=387.5$ MPa'),
               Simulation(directory=base_directory + '/fuller_weak/', line='--r', fig=1,
                          name=r'F{\"u}ller curve $\sigma_w=200$ MPa'),
               Simulation(directory=base_directory + '/8-16mm_continued/', line='b', fig=2,
                          name=r'8-16 mm $\sigma_w=387.5$ MPa'),
               Simulation(directory=base_directory + '/8-16mm_continued_weak/', line='--b', fig=2,
                          name=r'8-16 mm $\sigma_w=200$ MPa'),
               Simulation(directory=base_directory + '/8-16mm_one_layer/', line='b', fig=3,
                          name=r'8-16 mm $\sigma_w=387.5$ MPa'),
               Simulation(directory=base_directory + '/8-16mm_one_layer_weak/', line='--b', fig=3,
                          name=r'8-16 mm $\sigma_w=200$ MPa')]
figures = [Figure(xlim=125, ylim=400, leg_handles=[], label='(a)', name='proctor_8-16mm.png'),
           Figure(xlim=125, ylim=200, leg_handles=[], label='(a)', name='proctor_fuller.png'),
           Figure(xlim=200, ylim=600, leg_handles=[], label='(a)', name='proctor_continued.png'),
           Figure(xlim=125, ylim=300, leg_handles=[], label='(b)', name='proctor_one_layer.png')]

for simulation in simulations:
    plt.figure(simulation.fig)
    with open(simulation.directory + '/fractured_particles.pkl') as pickle_file:
        fracture_data = pickle.load(pickle_file)
    leg_h = plt.plot(fracture_data[:, 0], fracture_data[:, 3], simulation.line, lw=2, label=simulation.name)
    figures[simulation.fig].leg_handles.append(leg_h[0])


for fig_number, figure in enumerate(figures):
    plt.figure(fig_number)
    if fig_number in range(0, 3):
        for layer_count in [25, 50, 75, 100]:
            plt.plot([layer_count, layer_count], [0, figure.ylim], '--k')
    plt.xlabel('Blows', fontsize=24)
    plt.ylabel('Particle cracks', fontsize=24)
    plt.xlim(0, figure.xlim)
    plt.ylim(0, figure.ylim)
    ax = plt.subplot(111)
    plt.text(0.05, 0.9, figure.label, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    plt.legend(loc=(0.05, 0.6), bbox_transform=ax.transAxes)
    plt.tight_layout()
    plt.savefig(figure.name)

plt.show()
