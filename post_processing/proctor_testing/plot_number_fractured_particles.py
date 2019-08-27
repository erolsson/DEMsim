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
Simulation = namedtuple('Simulation', ['directory', 'color', 'name'])
simulations = [Simulation(directory=base_directory + '/8-16mm/', color='b', name='Fraction 8-16 mm'),
               Simulation(directory=base_directory + '/fuller1/', color='r', name=r'F{\"u}ller curve  0\% < 2 mm')]

for simulation in simulations:
    with open(simulation.directory + '/fractured_particles.pkl') as pickle_file:
        fracture_data = pickle.load(pickle_file)

    plt.plot(fracture_data[:, 0], fracture_data[:, 3], simulation.color, lw=2, label=simulation.name)
plt.show()
