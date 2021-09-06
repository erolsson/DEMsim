import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 24})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def main():
    simulation_directory = os.path.expanduser('~/DEMsim/results/stone_compaction/13mm/')
    particles = np.genfromtxt(simulation_directory + '/animation/particles/particles_0.01.dou', delimiter=',')
    compaction_force = np.genfromtxt(simulation_directory + '/compaction/surface_forces.dou', delimiter=',')
    surface_positions = np.genfromtxt(simulation_directory + '/compaction/surface_positions.dou', delimiter=',')
    F = compaction_force[:, -2]
    r = particles[:, 7]
    Vp = np.sum(4*np.pi*r**3/3)
    z = surface_positions[:, -2]
    R = surface_positions[:, 2]
    V = np.pi*R*R*z
    delta = surface_positions[0, -2] - z
    plt.figure(0)
    plt.plot((delta-0.01)*1000, F/1000, lw=2)
    plt.ylim(0, 40)
    plt.xlim(0, 10)
    plt.xlabel('Displacement [mm]')
    plt.ylabel('Force [kN]')
    plt.tight_layout()

    plt.figure(1)
    plt.plot(1 - Vp/V, F/1000, lw=2)
    plt.ylim(0, 40)
    plt.xlabel('Porosity [-]')
    plt.ylabel('Force [kN]')
    plt.tight_layout()

    plt.show()


if __name__ == '__main__':
    main()
