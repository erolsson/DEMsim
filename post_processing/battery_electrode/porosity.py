import os

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


def main():
    directory = os.path.expanduser('~/DEMsim/results/battery_rve/Cathod/')
    box_edges = np.genfromtxt(directory + 'periodic_bc.dou', delimiter=',')
    time = box_edges[:, 0]
    box_side = 2*box_edges[:, 1]
    surface_positions = np.genfromtxt(directory + 'surface_positions.dou', delimiter=',')
    box_height = surface_positions[:, -2]
    volume = box_side**2*box_height
    particles = np.genfromtxt(directory + 'particles/particles_20.096.dou', delimiter=',')
    r = particles[:, 7]
    v_part = np.sum(4*np.pi*r**3/3)
    relative_density = v_part/volume

    surface_forces = np.genfromtxt(directory + 'surface_forces.dou', delimiter=',')
    f = surface_forces[:, 6]
    plt.plot(relative_density, f)
    plt.show()


if __name__ == '__main__':
    main()
