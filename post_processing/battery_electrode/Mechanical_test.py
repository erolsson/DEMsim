
import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')


def dimensions_box(data_directory):
    file = open(data_directory + '/periodic_bc.dou',  'r')
    first_line = file.readlines()[0]
    first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line))]

    data = file[:,  2]
    print(data)
    return data


def box_height(data_directory):
    file = np.genfromtxt(data_directory + '/surface_positions.dou',  delimiter=', ')
    first_line = file.readlines()[0]
    first_line = first_line.split(', ')
    height = first_line[:, 34]
    return height


def pressures_box(data_directory):
    force_data = np.genfromtxt(data_directory + '/force_fabric_tensor.dou', delimiter=', ')
    with open(data_directory + '/force_fabric_tensor.dou', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')

    surface_tensor = first_line[:, 2]
    print(surface_tensor)
    return surface_tensor


if __name__ == '__main__':
    simulation_directory = '../../results/viscoelastic/procent65bt4E34PelleUtanFT/'
    box_width = 0.172726
    strain = (box_width - dimensions_box(simulation_directory))/box_width
    stress = pressures_box(simulation_directory)/(box_width * box_height() * dimensions_box()*2)
    plt.plot(stress, strain)
    plt.xlabel("Pressure [Pa]")
    plt.ylabel("Porosity")
    plt.show()
