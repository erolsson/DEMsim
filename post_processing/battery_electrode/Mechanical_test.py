
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
    #print(data)
    return data


def position_zz(data_directory):
    with open(data_directory + '/surface_positions.dou', 'r') as periodic_bc:
        first_line = periodic_bc.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx+1][5:] for idx in id_idx]
    if surface_types.count('ID=') == 0 and surface_types.count('PointSurface') == 2:
        id_idx.sort(key=lambda x: first_line[x+1])
        wall_data = np.genfromtxt(data_directory + '/surface_positions.dou', delimiter=', ')
        time = np.zeros((wall_data.shape[0], 1))
        time = wall_data[:, id_idx[0]+32]

    return time



def pressures_box(data_directory):
    with open(data_directory + '/force_fabric_tensor.dou', 'r') as force_fabric_tensor:
        first_line = force_fabric_tensor.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line))]
    force_data = np.genfromtxt(data_directory + '/force_fabric_tensor.dou', delimiter=', ')
    force = np.zeros((force_data.shape[0], 1))
    force = force_data[:,  id_idx[0]+1]
    return force

def pressures_box_yy(data_directory):
    with open(data_directory + '/force_fabric_tensor.dou', 'r') as force_fabric_tensor:
        first_line = force_fabric_tensor.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line))]
    force_data = np.genfromtxt(data_directory + '/force_fabric_tensor.dou', delimiter=', ')
    force = np.zeros((force_data.shape[0], 1))
    force = force_data[:,  id_idx[0]+5]
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
    simulation_directory = '../../results/viscoelastic/cubic_box/'
    box_width = 0.141403*2
    surface_height =  0.48568 # when the mechanical testing begins
    E = 2e9
    print(dimensions_box(simulation_directory))
    strain = -( 2* dimensions_box(simulation_directory)[:]-box_width)/box_width
    Stress = pressures_box(simulation_directory)[:]/(box_width * surface_height * box_width *2)
    Stress_y = pressures_box_yy(simulation_directory)[:]/(box_width * position_zz(simulation_directory)[:] *
                                                                   dimensions_box(simulation_directory)[:] *2)
    start_presure=1315
    stress = pressures_box(simulation_directory)[:]/(box_width * position_zz(simulation_directory)[:] *
                                                              dimensions_box(simulation_directory)[:] *2)
    stress1 = pressures_box(simulation_directory)[start_presure:start_presure+400]-pressures_box(simulation_directory)[start_presure] /(box_width * position_zz(simulation_directory)[start_presure:start_presure+400] *
                                                     dimensions_box(simulation_directory)[start_presure:start_presure+400] *2)
    strain1 = -( 2* dimensions_box(simulation_directory)[start_presure:start_presure+400]-box_width)/box_width

    stress2 = pressures_box(simulation_directory)[start_presure+400:start_presure+400+440]-pressures_box(simulation_directory)[start_presure+400]/(box_width * position_zz(simulation_directory)[start_presure+400:start_presure+400+440] *
                                                     dimensions_box(simulation_directory)[start_presure+400:start_presure+400+440] *2)
    strain2 = -( 2* dimensions_box(simulation_directory)[start_presure+400:start_presure+400+440]-box_width)/box_width

    stress3 = pressures_box(simulation_directory)[start_presure+400+440:start_presure+400+440+480]-pressures_box(simulation_directory)[start_presure+400+440]/(box_width * position_zz(simulation_directory)[start_presure+400+440:start_presure+400+440+480] *
                                                     dimensions_box(simulation_directory)[start_presure+400+440:start_presure+400+440+480] *2)
    strain3 = -( 2* dimensions_box(simulation_directory)[start_presure+400+440:start_presure+400+440+480]-box_width)/box_width

    stress4 = pressures_box(simulation_directory)[start_presure+400+440+480:start_presure+400+440+480+560]-pressures_box(simulation_directory)[start_presure+400+440+480]/(box_width * position_zz(simulation_directory)[start_presure+400+440+480+560] *
                                                     dimensions_box(simulation_directory)[start_presure+400+440+480+560] *2)
    strain4 = -( 2* dimensions_box(simulation_directory)[start_presure+400+440+480+560]-box_width)/box_width

    np.savetxt('stress.dat', stress)


    #inkompresibelt på binder, isotropiskt material+
    # inelastic strain
    epsilon_zz = -(position_zz(simulation_directory)[:]-surface_height)/surface_height
    print(position_zz(simulation_directory)[:])
    print(epsilon_zz)
    total_stress = Stress+Stress_y
    nu = -(E*epsilon_zz)/ total_stress
    print(nu)

    time = time_box(simulation_directory)[:]
    plt.plot(time, stress)
    plt.xlabel("time[s]")
    plt.ylabel("Stress [Pa]")
    plt.show()

    plt.plot(strain1, stress1)
    plt.plot(strain2, stress2)
    plt.plot(strain3, stress3)
    plt.plot(strain4, stress4)

    plt.xlabel("Strain")
    plt.ylabel("Stress [Pa]")
    plt.show()
    plt.plot(time, Stress_y)
    plt.xlabel("Time")
    plt.ylabel("Stress_y[Pa]")

    plt.show()
    plt.plot(time, nu)
    plt.xlabel("Time")
    plt.ylabel("nu")

    plt.show()


