
import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')


def dimensions_box(data_directory):
    with open(data_directory + '/periodic_bc.dou', 'r') as periodic_bc:
        first_line = periodic_bc.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line))]
    wall_data = np.genfromtxt(data_directory + '/periodic_bc.dou', delimiter=', ')
    data = np.zeros((wall_data.shape[0], 1))
    data = wall_data[:,  id_idx[0]+2]
    print(data)
    return data


def pressures_box(data_directory):
    with open(data_directory + '/force_fabric_tensor.dou', 'r') as force_fabric_tensor:
        first_line = force_fabric_tensor.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line))]
    force_data = np.genfromtxt(data_directory + '/force_fabric_tensor.dou', delimiter=', ')
    force = np.zeros((force_data.shape[0], 1))
    force = force_data[:,  id_idx[0]+1]
    return force


def time_box(data_directory):
    with open(data_directory + '/force_fabric_tensor.dou', 'r') as force_fabric_tensor:
        first_line = force_fabric_tensor.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line))]
    force_data = np.genfromtxt(data_directory + '/force_fabric_tensor.dou', delimiter=', ')
    time = np.zeros((force_data.shape[0], 1))
    time = force_data[:,  id_idx[0]]
    return time


if __name__ == '__main__':
    simulation_directory = '../../results/viscoelastic/relaxation_test/'
    box_width = 0.172726
    box_height = 0.733144
    strain = (box_width - dimensions_box(simulation_directory))[:]/box_width
    Stress = pressures_box(simulation_directory)[:]/(box_width * box_height * box_width *2)
    stress = pressures_box(simulation_directory)[:]/(box_width * box_height *
                                                               dimensions_box(simulation_directory)[:] *2)
    #inkompresibelt på binder, isotropiskt material
    # inelastic strain
    time = time_box(simulation_directory)[:]
    plt.plot(time, Stress)
    plt.xlabel("time[s]")
    plt.ylabel("Stress [Pa]")
    plt.show()
    plt.plot(strain, stress)
    plt.xlabel("Strain")
    plt.ylabel("Stress [Pa]")
    plt.show()

