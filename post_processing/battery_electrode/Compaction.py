import glob

from math import pi

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')



def particle_volume():
    p_volume = 0.0267963
    return p_volume


def dimensions_box(data_directory):
    with open(data_directory + '/surface_positions.dou', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx+1][5:] for idx in id_idx]
    if surface_types.count('ID=') == 0 and surface_types.count('PointSurface') == 2:
        id_idx.sort(key=lambda x: first_line[x+1])
        wall_data = np.genfromtxt(data_directory + '/surface_positions.dou', delimiter=', ')
        data = np.zeros((wall_data.shape[0], 1))
        data = wall_data[1:79700,  id_idx[0]+32]
        print(data)
        return data
    else:
        raise ValueError("A box could not be defined from the data in " + data_directory +
                         '/surface_positions.dou')


def pressures_box(data_directory):
    force_data = np.genfromtxt(data_directory + '/surface_forces.dou', delimiter=', ')
    with open(data_directory + '/surface_positions.dou', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
        id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx+1][5:] for idx in id_idx]

    surface_indices = [i for i, surface_type in enumerate(surface_types) if surface_type == 'PointSurface']

    surface1_force = force_data[1:79700, surface_indices[0]+1]

    p = surface1_force/(0.34*0.34)
    print(p)
    return p


if __name__ == '__main__':
    simulation_directory = '../../results/viscoelastic/procent70bt4E34height45-1000/'

    volume_box = (2*0.172726)**2 * dimensions_box(simulation_directory)
    porosity = (1-(particle_volume()*(1+0.07/(0.33+0.07)))/volume_box)
    pressures = pressures_box(simulation_directory)
    plt.plot(pressures, porosity*100)
    plt.xlabel("Pressure [Pa]")
    plt.ylabel("Porosity")
    plt.show()

