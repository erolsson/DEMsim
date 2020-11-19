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
    particle_files = glob.glob(particle_directory + '/particles_*.dou')
    particle_data = np.genfromtxt(particle_files[0], delimiter=', ')
    r = particle_data[:, 7]
    return 0.0267963


def dimensions_box(data_directory):
    with open(data_directory + '/surface_positions.dou', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx+1][5:] for idx in id_idx]
    if surface_types.count('Cylinder') == 1 and surface_types.count('PointSurface') == 2:
        id_idx.sort(key=lambda x: first_line[x+1])
        wall_data = np.genfromtxt(data_directory + '/surface_positions.dou', delimiter=', ')
        data = np.zeros((wall_data.shape[0], 4))
        data[:, 0] = wall_data[:, -1]
        data[:, 1] = wall_data[:, id_idx[0]+2]
        data[:, 2] = np.mean(wall_data[:, id_idx[1]+5:id_idx[1]+15:3], 1)
        data[:, 3] = np.mean(wall_data[:, id_idx[2]+5:id_idx[2]+15:3], 1)
        return data
    else:
        raise ValueError("A cylinder could not be defined from the data in " + data_directory +
                         '/surface_positions.dou')


def volume_cylinder(data_directory):
    dimensions = dimensions_box(data_directory)
    volume = pi*dimensions[:, 1]**2*np.abs(dimensions[:, 3] - dimensions[:, 2])
    return volume


def relative_density_cylinder(data_directory):
    data = volume_cylinder(data_directory)
    data[:, 1] = particle_volume(data_directory) / volume_cylinder(data_directory)[:, 1]
    return data


def pressures_cylinder(data_directory):
    dimensions = dimensions_box(data_directory)
    force_data = np.genfromtxt(data_directory + '/surface_forces.dou', delimiter=', ')
    r = dimensions[:, 1]
    z0 = dimensions[:, 2]
    z1 = dimensions[:, 3]

    with open(data_directory + '/surface_positions.dou', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
        id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx+1][5:] for idx in id_idx]

    cyl_idx = [i for i, surface_type in enumerate(surface_types) if surface_type == 'Cylinder'][0]
    surface_indices = [i for i, surface_type in enumerate(surface_types) if surface_type == 'PointSurface']
    cylinder_force = force_data[:, cyl_idx*5+1]
    surface1_force = force_data[:, surface_indices[0]*5+1]
    surface2_force = force_data[:, surface_indices[1]*5+1]

    p = 1*dimensions
    p[:, 1] = cylinder_force/(2*pi*r*np.abs(z0-z1))
    p[:, 2] = surface1_force/(pi*r**2)
    p[:, 3] = surface2_force/(pi*r**2)
    return p


if __name__ == '__main__':
    simulation_directory = '../results/viscoelastic/11-18-kvall/'
    v_particles = particle_volume(simulation_directory)
    for directory in ["", "unload_D=0.65", "unload_D=0.7", "unload_D=0.75", "unload_D=0.8"]:
        v_cylinder = volume_cylinder(simulation_directory + directory)
        pressures = pressures_cylinder(simulation_directory + directory)
        relative_density = v_particles/v_cylinder
        plt.plot(relative_density, pressures[:, 1:])
    plt.show()
