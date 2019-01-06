import glob

from math import pi

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


def particle_volume(particle_directory):
    particle_files = glob.glob(particle_directory + '/particles_*.dat')
    particle_data = np.genfromtxt(particle_files[0], delimiter=', ')
    r = particle_data[:, 7]
    return np.sum(4*pi*r**3/3)


def cylinder_volume(data_directory):
    with open(data_directory + '/surface_positions.dat', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx+1][5:] for idx in id_idx]
    if surface_types.count('Cylinder') == 1 and surface_types.count('PointSurface') == 2:
        id_idx.sort(key=lambda x: first_line[x+1])
        wall_data = np.genfromtxt(data_directory + '/surface_positions.dat', delimiter=', ')
        r = wall_data[:, id_idx[0]+2]
        z0 = np.mean(wall_data[:, id_idx[1]+5:id_idx[1]+15:3], 1)
        z1 = np.mean(wall_data[:, id_idx[2]+5:id_idx[2]+15:3], 1)
        return pi*r**2*np.abs(z0-z1)

    else:
        raise ValueError("A cylinder could not be defined from the data in " + data_directory +
                         '/surface_positions.dat')


def cylinder_pressures(data_directory):
    with open(data_directory + '/surface_positions.dat', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx+1][5:] for idx in id_idx]
    if surface_types.count('Cylinder') == 1 and surface_types.count('PointSurface') == 2:
        id_idx.sort(key=lambda index: first_line[index+1])
        wall_data = np.genfromtxt(data_directory + '/surface_positions.dat', delimiter=', ')
        force_data = np.genfromtxt(data_directory + '/surface_forces.dat', delimiter=', ')
        r = wall_data[:, id_idx[0]+2]
        z0 = np.mean(wall_data[:, id_idx[1]+5:id_idx[1]+15:3], 1)
        z1 = np.mean(wall_data[:, id_idx[2]+5:id_idx[2]+15:3], 1)

        cyl_idx = [i for i, surface_type in enumerate(surface_types) if surface_type == 'Cylinder'][0]
        surface_indices = [i for i, surface_type in enumerate(surface_types) if surface_type == 'PointSurface']
        cylinder_force = force_data[:, cyl_idx*5+1]
        surface1_force = force_data[:, surface_indices[0]*5+1]
        surface2_force = force_data[:, surface_indices[1]*5+1]

        p = np.zeros((r.shape[0], 3))
        p[:, 0] = cylinder_force/(2*pi*r*np.abs(z0-z1))
        p[:, 1] = surface1_force/(2*pi*r**2)
        p[:, 2] = surface2_force/(2*pi*r**2)
        return p

    else:
        raise ValueError("A cylinder could not be defined from the data in " + data_directory +
                         '/surface_positions.dat')


if __name__ == '__main__':
    simulation_directory = '../results/closed_die_compaction/3'
    v_particles = particle_volume(simulation_directory)

    v_cylinder = cylinder_volume(simulation_directory)
    pressures = cylinder_pressures(simulation_directory)
    relative_density = v_particles/v_cylinder
    plt.plot(relative_density, pressures)
    plt.show()
