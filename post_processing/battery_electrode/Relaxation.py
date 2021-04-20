import numpy as np
import os

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')


def epsilon_box_compression(simulation_directory_compression):

    force_data = np.genfromtxt(simulation_directory_compression + '/force_fabric_tensor.dou', delimiter=', ')
    time= force_data[1007:,0]-14.0543
    pressures_box= force_data[1007:, 1]

    Stress_compression= (pressures_box[:]/pressures_box[0])+0.01
    print(Stress_compression)

    plt.plot(time, Stress_compression, 'b', label='DEM data-Compression')
    plt.xlabel("time[s]")
    plt.ylabel("Stress")

    epsilon = 0.0175
    t = np.arange(0,215)
    relaxation_compression = 0.344+0.097 * np.exp(-1*t/283)+ 0.074* np.exp(-1*t/6770)
    Sigma_compression = relaxation_compression/relaxation_compression[0]

    plt.plot(t,Sigma_compression, 'bx',label='Fit for Prony series-Compression')


def epsilon_box_tension(simulation_directory_tension):

    force_data = np.genfromtxt(simulation_directory_tension + '/force_fabric_tensor.dou', delimiter=', ')
    time= force_data[1007:,0]-14.0543
    pressures_box= force_data[1007:, 1]
    Stress_tension= pressures_box[:]/pressures_box[0]

    plt.plot(time, Stress_tension, 'r',label='DEM data-Tension')
    plt.xlabel("time[s]")
    plt.ylabel("Normalized Stress ")
    epsilon = 0.0175
    t = np.arange(0,215)
    relaxation_tension = 0.117+0.065 * np.exp(-1*t/211)+ 0.057* np.exp(-1*t/4807)
    Sigma_tension = relaxation_tension /relaxation_tension[0]

    plt.plot(t,Sigma_tension, 'rx',label='Fit for Prony series-Tension')





if __name__ == '__main__':


    simulation_directory = os.path.expanduser(r'C:/DEMsim/results/relaxation-compression/')
    epsilon_box_compression(simulation_directory)

    simulation_directory_tension= os.path.expanduser(r'C:/DEMsim/results/relaxation-tension/')
    epsilon_box_tension(simulation_directory_tension)

    plt.legend(loc='upper right', numpoints =1)
    plt.show()



