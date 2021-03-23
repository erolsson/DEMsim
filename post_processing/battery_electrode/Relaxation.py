import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')





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
    simulation_directory = 'C:/DEMsim/results/relaxaiton'
    box_width = 0.448928 *2
    surface_height = 0.882022
    E = 1.6e9

    Stress = pressures_box(simulation_directory)[:]/(box_width * surface_height * box_width *2)



    time = time_box(simulation_directory)[:]
    plt.plot(time, Stress/1e6, label='DEM')
    plt.xlabel("time[s]")
    plt.ylabel("Stress [MPa]")
    epsilon = 0.0175
    t = np.arange(14,300)
    relaxation = 0.117+0.065 * np.exp(-1*t/211)+ 0.057* np.exp(-1*t/4807)
    #Sigma_DEM = 0.239+0.272*np.exp(-1*t/211)+0.2385*np.exp(-1*t/4807)
    Sigma = relaxation *epsilon*3e9
    plt.plot(t,Sigma/1e6, label='Theory')
    plt.legend()
    plt.show()


