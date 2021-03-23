
import numpy as np
from sympy import symbols, solve
import matplotlib.pyplot as plt
import matplotlib
from sympy import symbols

matplotlib.style.use('classic')


def dimensions_xx(data_directory):
    with open(data_directory + '/periodic_bc.dou', 'r') as periodic_bc:
        first_line = periodic_bc.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line))]
    wall_data = np.genfromtxt(data_directory + '/periodic_bc.dou', delimiter=', ')
    data = np.zeros((wall_data.shape[0], 1))
    data = wall_data[:,  id_idx[0]+2]
    #print(data)
    return data


def dimensions_zz(data_directory):
    with open(data_directory + '/surface_positions.dou', 'r') as periodic_bc:
        first_line = periodic_bc.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx+1][5:] for idx in id_idx]
    if surface_types.count('ID=') == 0 and surface_types.count('PointSurface') == 2:
        id_idx.sort(key=lambda x: first_line[x+1])
        wall_data = np.genfromtxt(data_directory + '/surface_positions.dou', delimiter=', ')
        Height = np.zeros((wall_data.shape[0], 1))
        Height = wall_data[:, id_idx[0]+32]

    return Height



def pressures_xx(data_directory):
    with open(data_directory + '/force_fabric_tensor.dou', 'r') as force_fabric_tensor:
        first_line = force_fabric_tensor.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line))]
    force_data = np.genfromtxt(data_directory + '/force_fabric_tensor.dou', delimiter=', ')
    force = np.zeros((force_data.shape[0], 1))
    force = force_data[:,  1]
    return force

def pressures_yy(data_directory):
    with open(data_directory + '/force_fabric_tensor.dou', 'r') as force_fabric_tensor:
        first_line = force_fabric_tensor.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line))]
    force_data = np.genfromtxt(data_directory + '/force_fabric_tensor.dou', delimiter=', ')
    force = np.zeros((force_data.shape[0], 1))
    force = force_data[:,  5]
    return force

def pressures_zz(data_directory):
    with open(data_directory + '/force_fabric_tensor.dou', 'r') as force_fabric_tensor:
        first_line = force_fabric_tensor.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line))]
    force_data = np.genfromtxt(data_directory + '/force_fabric_tensor.dou', delimiter=', ')
    force = np.zeros((force_data.shape[0], 1))
    force = force_data[:, 9]
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
    simulation_directory = 'C:/DEMsim/results/compression-tension'
    #165 495  775 1003 1249 1460 1669 1819 1919
    lln=3061
    cln=1919+25

    start=lln-cln

    volym= dimensions_xx(simulation_directory)[start:start+20]*dimensions_xx(simulation_directory)[1]*dimensions_zz(simulation_directory)[start:start+20]*4
    time = time_box(simulation_directory)[start:start+20]
    sigma_x = pressures_xx(simulation_directory)[start:start+20]/volym
    delta_sigma_x=(sigma_x[:]-sigma_x[0])
    sigma_y = pressures_yy(simulation_directory)[start:start+20]/volym
    delta_sigma_y=sigma_y[:]-sigma_y[0]
    sigma_z= (pressures_zz(simulation_directory)[start:start+20])/volym
    delta_sigma_z=sigma_z[:]-sigma_z[0]
    print(delta_sigma_z)
    print(delta_sigma_y)
    print(delta_sigma_x)

    nu = delta_sigma_y/(delta_sigma_x+delta_sigma_z)


    epsilon_x = (dimensions_xx(simulation_directory)[start:start+20]-dimensions_xx(simulation_directory)[start-1] )/dimensions_xx(simulation_directory)[start-1]
    print(epsilon_x)
    delta_epsilon_x= (epsilon_x[:]-epsilon_x[0])

    sigma=sigma_x[1:19]- nu[1:19]*(sigma_y[1:19]+sigma_z[1:19])
    delta_sigma= -(sigma-sigma[0])
    print(delta_sigma)


    E = (delta_sigma/delta_epsilon_x[1:19])
    print(nu)
    print (E)
    plt.plot(delta_epsilon_x[1:19],nu[1:19])

    plt.xlabel("delta_strain")
    plt.ylabel("nu")
    plt.show()
    end= 1000

    volym_tot = dimensions_xx(simulation_directory)[end:start]*dimensions_xx(simulation_directory)[1]*dimensions_zz(simulation_directory)[end:start]*4
    stress = pressures_xx(simulation_directory)[end:start]/volym_tot
    tid = time_box(simulation_directory)[end:start]
    strain = (dimensions_xx(simulation_directory)[end:start]-dimensions_xx(simulation_directory)[start-1] )/dimensions_xx(simulation_directory)[start-1]
    plt.plot(strain,stress/10**6,label='DEM')
    plt.xlabel("strain")
    plt.ylabel("stress x-direction[MPa]")
    plt.xlabel("delta_strain")
    plt.ylabel("E")
    plt.show()





    volym_tot = dimensions_xx(simulation_directory)[end:start]*dimensions_xx(simulation_directory)[1]*dimensions_zz(simulation_directory)[end:start]*4
    stress = pressures_xx(simulation_directory)[end:start]/volym_tot
    tid = time_box(simulation_directory)[end:start]
    strain = -(dimensions_xx(simulation_directory)[end:start]-dimensions_xx(simulation_directory)[3200] )/dimensions_xx(simulation_directory)[3200]
    plt.plot(strain,stress/10**6,label='DEM')
    plt.xlabel("strain")
    plt.ylabel("stress x-direction[MPa]")



    plt.show()







    sigma_x_tot = pressures_xx(simulation_directory)[1277:1297]/(dimensions_xx(simulation_directory)[1277] *2* dimensions_zz(simulation_directory)[1277:1297] *
                                                              dimensions_xx(simulation_directory)[1277:1297] *2)
    strain = -( dimensions_xx(simulation_directory)[1277:1297]-dimensions_xx(simulation_directory)[1277])/dimensions_xx(simulation_directory)[1277]

    plt.plot(strain,sigma_x_tot)
    plt.xlabel("Strain")
    plt.ylabel("Stress [Pa]")
    plt.show()

    time_tot = time_box(simulation_directory)[1277:1297]


    plt.plot(time_tot, sigma_x_tot)
    plt.xlabel("time[s]")
    plt.ylabel("Stress [Pa]")
    plt.show()

    quit()



























    start_presure=5321

    stress1 = pressures_xx(simulation_directory)[start_presure:start_presure+200]/(dimensions_xx(simulation_directory)[1277]  * dimensions_zz(simulation_directory)[start_presure:start_presure+200] *
                                                                                    dimensions_xx(simulation_directory)[start_presure:start_presure+200] *2)
    strain1 = -( 2* dimensions_xx(simulation_directory)[start_presure:start_presure+200]-dimensions_xx(simulation_directory)[1277] )/dimensions_xx(simulation_directory)[1277]

    stress2 = pressures_xx(simulation_directory)[start_presure+400:start_presure+400+220]/(dimensions_xx(simulation_directory)[1277]  * dimensions_zz(simulation_directory)[start_presure+400:start_presure+400+220] *
                                                                                            dimensions_xx(simulation_directory)[start_presure+400:start_presure+400+220] *2)
    strain2 = -( 2* dimensions_xx(simulation_directory)[start_presure+400:start_presure+400+220]-dimensions_xx(simulation_directory)[1277] )/dimensions_xx(simulation_directory)[1277]

    stress3 = pressures_xx(simulation_directory)[start_presure+400+440:start_presure+400+440+240]/(dimensions_xx(simulation_directory)[1277]  * dimensions_zz(simulation_directory)[start_presure+400+440:start_presure+400+440+240] *
                                                                                                    dimensions_xx(simulation_directory)[start_presure+400+440:start_presure+400+440+240] *2)
    strain3 = -( 2* dimensions_xx(simulation_directory)[start_presure+400+440:start_presure+400+440+240]-dimensions_xx(simulation_directory)[1277] )/dimensions_xx(simulation_directory)[1277]

    stress4 = pressures_xx(simulation_directory)[start_presure+400+440+480:start_presure+400+440+480+280]/(dimensions_xx(simulation_directory)[1277]  * dimensions_zz(simulation_directory)[start_presure+400+440+480:start_presure+400+440+480+280] *
                                                                                                            dimensions_xx(simulation_directory)[start_presure+400+440+480:start_presure+400+440+480+280] *2)
    strain4 = -( 2* dimensions_xx(simulation_directory)[start_presure+400+440+480:start_presure+400+440+480+280]-dimensions_xx(simulation_directory)[1277] )/dimensions_xx(simulation_directory)[1277]

    plt.plot(strain1, stress1[1277:1297]-stress1[0])
    plt.plot(strain2, stress2)
    plt.plot(strain3, stress3)
    plt.plot(strain4, stress4)

    plt.xlabel("Strain")
    plt.ylabel("Stress [Pa]")
    plt.show()


    plt.show()


