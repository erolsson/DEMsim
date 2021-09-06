import numpy as np
import os

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')


def epsilon_box_compression(simulation_directory_compression):

    force_data = np.genfromtxt(simulation_directory_compression + '/force_fabric_tensor.dou', delimiter=', ')
    time= force_data[1609:,0]-9
    pressures_box= force_data[1609:, 1]
    print(pressures_box)
    print(pressures_box[0])

    Stress_compression= (pressures_box[:]/pressures_box[0])
    print(Stress_compression)

    plt.plot(time, Stress_compression, '--k', label='DEM data-Compression $Fraction of binder contacts=0.8$')
    plt.xlabel("time[s]")
    plt.ylabel("Stress")
    plt.legend()




def epsilon_box_tension(simulation_directory_tension):

    force_data = np.genfromtxt(simulation_directory_tension + '/force_fabric_tensor.dou', delimiter=', ')
    time= force_data[1422:,0]-9
    pressures_box= force_data[1422:, 1]
    Stress_tension= pressures_box[:]/pressures_box[0]-0.009

    plt.plot(time, Stress_tension+0.01, '--y',label='DEM data-tension $Fraction of binder contacts=0.8$')
    plt.xlabel("time[s]")
    plt.ylabel("Normalized Stress ")

    plt.legend()



def epsilon_box_compression2(simulation_directory_compression2):

    force_data = np.genfromtxt(simulation_directory_compression2 + '/force_fabric_tensor.dou', delimiter=', ')
    time= force_data[1725:,0]-9
    pressures_box= force_data[1725:, 1]

    Stress_compression= (pressures_box[:]/pressures_box[0])-0.01895
    print(Stress_compression)

    plt.plot(time, Stress_compression, '--b', label='DEM data-compression $Fraction of binder contacts=0.4$')
    plt.xlabel("time[s]")
    plt.ylabel("Stress")




def epsilon_box_tension2(simulation_directory_tension2):

    force_data = np.genfromtxt(simulation_directory_tension2 + '/force_fabric_tensor.dou', delimiter=', ')
    time= force_data[1500:,0]-9
    pressures_box= force_data[1500:, 1]
    Stress_tension= pressures_box[:]/pressures_box[0]-0.018

    plt.plot(time, Stress_tension+0.017, '--r',label='DEM data-tension $Fraction of binder contacts=0.4$')
    plt.xlabel("time[s]")
    plt.ylabel("Normalized Stress ")


if __name__ == '__main__':


    simulation_directory_compression = os.path.expanduser(r'C:/DEMsim/results/relaxation-compression/E34bt01fraction08/')
    epsilon_box_compression(simulation_directory_compression)

    simulation_directory_tension= os.path.expanduser(r'C:/DEMsim/results/relaxation-tension/E34bt01fraction08/')
    epsilon_box_tension(simulation_directory_tension)


    simulation_directory2 = os.path.expanduser(r'C:/DEMsim/results/relaxation-compression/E34bt01fraction04/')
    epsilon_box_compression2(simulation_directory2)

    simulation_directory_tension2= os.path.expanduser(r'C:/DEMsim/results/relaxation-tension/E34bt01fraction04/')
    epsilon_box_tension2(simulation_directory_tension2)


    epsilon = 0.003544
    t = np.arange(0,215)
    relaxation_compression = 0.344+0.097 * np.exp(-1*t/283)+ 0.074* np.exp(-1*t/6770)
    Sigma_compression = relaxation_compression/relaxation_compression[0]
    plt.plot(t,Sigma_compression, 'bx',label='Fit for Prony series-Compression')

    epsilon = 0.003544
    t = np.arange(0,215)
    relaxation_tension = 0.117+0.065 * np.exp(-1*t/211)+ 0.057* np.exp(-1*t/4807)
    Sigma_tension = relaxation_tension/relaxation_tension[0]
    plt.plot(t,Sigma_tension, 'rx',label='Fit for Prony series-Tension')

    plt.legend(loc='upper right', numpoints=1)
    plt.xlim(-0.4,238)
    plt.show()




