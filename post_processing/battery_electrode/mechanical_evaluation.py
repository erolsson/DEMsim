import glob
import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib
from multiprocesser.multiprocesser import multi_processer

matplotlib.style.use('classic')


def plot_mechanical_data_for_simulation(directory):
    periodic_bc = np.genfromtxt(directory + 'periodic_bc.dou', delimiter=',')
    surface_positions = np.genfromtxt(directory + 'surface_positions.dou', delimiter=',')
    force_fabric_tensor = np.genfromtxt(directory + 'force_fabric_tensor.dou', delimiter=',')
    time = surface_positions[:, -1]
    d = 2*periodic_bc[:, 2]
    d0 = d[0]
    w = 2*periodic_bc[:, 4]

    # Finding the point where the compaction force has decreased to zero after its maximum value
    # The thickness at this point will be the thickness of the electrode t0
    t0 = 0.650098

    t_start = time[d == d0][-1]
    volume = (d*w*t0)[time > t_start]
    sxx = force_fabric_tensor[time > t_start, 1]/volume
    syy = force_fabric_tensor[time > t_start, 5]/volume
    szz = force_fabric_tensor[time > t_start, 9]/volume

    sxx *= -1
    syy *= -1
    szz *= -1

    linear_strain = (d[time > t_start] - d0)/d0
    e0_vec=[]

    E_vec= []
    v_vec=[]

    plt.figure(0)
    if sxx[-1] > 0:
        plt.plot(linear_strain, sxx/1e6, 'r', lw=2, label=r'$\sigma_{xx}$')
        plt.plot(linear_strain, syy/1e6, 'b', lw=2, label=r'$\sigma_{yy}$')
        plt.plot(linear_strain, szz/1e6, 'g', lw=2, label=r'$\sigma_{zz}$')
    else:
        plt.plot(linear_strain, sxx/1e6, 'r', lw=2)
        plt.plot(linear_strain, syy/1e6, 'b', lw=2)
        plt.plot(linear_strain, szz/1e6, 'g', lw=2)

    idx = [1]
    delta_strain = 3e-3
    while len(idx):
        # This finds the turning points
        bool_arr = np.logical_and(np.diff(np.abs(linear_strain)) < 0,
                                  np.abs(linear_strain[1:]) > np.abs(linear_strain[idx[0]]))

        idx = np.where(bool_arr)[0]
        if len(idx):
            e0 = linear_strain[idx[0]]
            e1 = e0 + np.abs(e0)/e0*delta_strain           # increasing the magnitude of the strain with delta_strain
            idx1 = np.argmin(np.abs(linear_strain - e1))   # Finding the index of the strain point closest to e1
            e1 = linear_strain[idx1]                       # Grabbing the correct e1
            dsxx = sxx[idx1] - sxx[idx[0]]
            dsyy = syy[idx1] - syy[idx[0]]
            dszz = szz[idx1] - szz[idx[0]]

            dexx = e1 - e0
            if dexx != 0:

                v = dsyy/(dsxx + dszz)
                plt.figure(1)
                plt.plot(e0, v, '-kx', ms=12, mew=2)

                E = (dsxx - v*(dsyy + dszz))/dexx
                plt.figure(2)
                e0_vec.append(e0)

                E_vec.append(E/1e9)
                v_vec.append(v)

    plt.plot(e0_vec, E_vec, '--kx',ms= 12, mew=2)
    return e0_vec[0],E_vec[0]

#                e0_exp_comp = np.array([-0.01,-0.011,-0.0123,-0.0141,-0.0165])
 #               E_exp_comp = np.array([1.20, 1.43, 1.51, 1.55, 1.99])
  #              plt.plot(e0_exp_comp, E_exp_comp, 'bo', ms=12)
   #             e0_exp_ten = np.array([0.01,0.011,0.0123,0.0141,0.0165])
    #            E_exp_ten = np.array([0.95, 0.78, 0.84, 0.76, 1.06])
     #           plt.plot(e0_exp_ten, E_exp_ten, 'bo', ms=12 )
    #plt.plot(e0_exp_ten, E_exp_ten, 'bo', ms=12, label=r'Experiment data')


def plot_mechanical_data_for_simulation1(directory1):
    periodic_bc = np.genfromtxt(directory1 + 'periodic_bc.dou', delimiter=',')
    surface_positions = np.genfromtxt(directory1 + 'surface_positions.dou', delimiter=',')
    force_fabric_tensor = np.genfromtxt(directory1 + 'force_fabric_tensor.dou', delimiter=',')
    time = surface_positions[:, -1]
    d = 2*periodic_bc[:, 2]
    d0 = d[0]
    w = 2*periodic_bc[:, 4]

    # Finding the point where the compaction force has decreased to zero after its maximum value
    # The thickness at this point will be the thickness of the electrode t0
    t0 = 0.650098

    t_start = time[d == d0][-1]
    volume = (d*w*t0)[time > t_start]
    sxx = force_fabric_tensor[time > t_start, 1]/volume
    syy = force_fabric_tensor[time > t_start, 5]/volume
    szz = force_fabric_tensor[time > t_start, 9]/volume

    sxx *= -1
    syy *= -1
    szz *= -1

    linear_strain = (d[time > t_start] - d0)/d0
    e0_vec=[]

    E_vec= []
    v_vec=[]

    plt.figure(0)
    if sxx[-1] > 0:
        plt.plot(linear_strain, sxx/1e6, 'r', lw=2, label=r'$\sigma_{xx}$')
        plt.plot(linear_strain, syy/1e6, 'b', lw=2, label=r'$\sigma_{yy}$')
        plt.plot(linear_strain, szz/1e6, 'g', lw=2, label=r'$\sigma_{zz}$')
    else:
        plt.plot(linear_strain, sxx/1e6, 'r', lw=2)
        plt.plot(linear_strain, syy/1e6, 'b', lw=2)
        plt.plot(linear_strain, szz/1e6, 'g', lw=2)

    idx = [1]
    delta_strain = 3e-3
    while len(idx):
        # This finds the turning points
        bool_arr = np.logical_and(np.diff(np.abs(linear_strain)) < 0,
                                  np.abs(linear_strain[1:]) > np.abs(linear_strain[idx[0]]))

        idx = np.where(bool_arr)[0]
        if len(idx):
            e0 = linear_strain[idx[0]]
            e1 = e0 + np.abs(e0)/e0*delta_strain           # increasing the magnitude of the strain with delta_strain
            idx1 = np.argmin(np.abs(linear_strain - e1))   # Finding the index of the strain point closest to e1
            e1 = linear_strain[idx1]                       # Grabbing the correct e1
            dsxx = sxx[idx1] - sxx[idx[0]]
            dsyy = syy[idx1] - syy[idx[0]]
            dszz = szz[idx1] - szz[idx[0]]

            dexx = e1 - e0
            if dexx != 0:

                v = dsyy/(dsxx + dszz)
                plt.figure(1)
                plt.plot(e0, v, '--rx', ms= 12, mew=2)

                E = (dsxx - v*(dsyy + dszz))/dexx
                plt.figure(2)
                e0_vec.append(e0)

                E_vec.append(E/1e9)
                v_vec.append(v)
    plt.plot(e0_vec, E_vec, '--rx',ms= 12, mew=2)
    return e0_vec[0],E_vec[0]


def plot_mechanical_data_for_simulation2(directory2):
    periodic_bc = np.genfromtxt(directory2 + 'periodic_bc.dou', delimiter=',')
    surface_positions = np.genfromtxt(directory2 + 'surface_positions.dou', delimiter=',')
    force_fabric_tensor = np.genfromtxt(directory2 + 'force_fabric_tensor.dou', delimiter=',')
    time = surface_positions[:, -1]
    d = 2*periodic_bc[:, 2]
    d0 = d[0]
    w = 2*periodic_bc[:, 4]

    # Finding the point where the compaction force has decreased to zero after its maximum value
    # The thickness at this point will be the thickness of the electrode t0
    t0 = 0.650098

    t_start = time[d == d0][-1]
    volume = (d*w*t0)[time > t_start]
    sxx = force_fabric_tensor[time > t_start, 1]/volume
    syy = force_fabric_tensor[time > t_start, 5]/volume
    szz = force_fabric_tensor[time > t_start, 9]/volume

    sxx *= -1
    syy *= -1
    szz *= -1

    linear_strain = (d[time > t_start] - d0)/d0
    e0_vec=[]

    E_vec= []
    v_vec=[]

    plt.figure(0)
    if sxx[-1] > 0:
        plt.plot(linear_strain, sxx/1e6, 'r', lw=2, label=r'$\sigma_{xx}$')
        plt.plot(linear_strain, syy/1e6, 'b', lw=2, label=r'$\sigma_{yy}$')
        plt.plot(linear_strain, szz/1e6, 'g', lw=2, label=r'$\sigma_{zz}$')
    else:
        plt.plot(linear_strain, sxx/1e6, 'r', lw=2)
        plt.plot(linear_strain, syy/1e6, 'b', lw=2)
        plt.plot(linear_strain, szz/1e6, 'g', lw=2)

    idx = [1]
    delta_strain = 3e-3
    while len(idx):
        # This finds the turning points
        bool_arr = np.logical_and(np.diff(np.abs(linear_strain)) < 0,
                                  np.abs(linear_strain[1:]) > np.abs(linear_strain[idx[0]]))

        idx = np.where(bool_arr)[0]

        if len(idx):
            e0= linear_strain[idx[0]]
            e1 = e0 + np.abs(e0)/e0*delta_strain           # increasing the magnitude of the strain with delta_strain
            idx1 = np.argmin(np.abs(linear_strain - e1))   # Finding the index of the strain point closest to e1
            e1 = linear_strain[idx1]                       # Grabbing the correct e1
            dsxx = sxx[idx1] - sxx[idx[0]]
            dsyy = syy[idx1] - syy[idx[0]]
            dszz = szz[idx1] - szz[idx[0]]
            dexx = e1 - e0

            if dexx != 0:

                v= dsyy/(dsxx + dszz)


                plt.figure(1)
                plt.plot(e0, v, 'bx',ms= 12, mew=2)
                E = (dsxx - v*(dsyy + dszz))/dexx
                plt.figure(2)
                e0_vec.append(e0)

                E_vec.append(E/1e9)
                v_vec.append(v)

    plt.plot(e0_vec, E_vec, '--bx',ms= 12, mew=2)
    return e0_vec[0],E_vec[0]


def main():

    directory = os.path.expanduser(r'C:/DEMsim/results/E34bt005rb05/Compression/')
    plot_mechanical_data_for_simulation(directory)
    directory = os.path.expanduser(r'C:/DEMsim/results/E34bt005rb05/Tension/')
    plot_mechanical_data_for_simulation(directory)

    directory1 = os.path.expanduser(r'C:/DEMsim/results/E34bt02rb05/Compression/')
    plot_mechanical_data_for_simulation1(directory1)
    directory1 = os.path.expanduser(r'C:/DEMsim/results/E34bt02rb05/Tension/')
    plot_mechanical_data_for_simulation1(directory1)

    directory2 = os.path.expanduser(r'C:/DEMsim/results/E34bt01rb05/tension/')
    plot_mechanical_data_for_simulation2(directory2)
    directory2 = os.path.expanduser(r'C:/DEMsim/results/E34bt01rb05/compression/')
    plot_mechanical_data_for_simulation2(directory2)


    plt.figure(0)
    plt.xlabel('Strain [-]')
    plt.ylabel('Stress [MPa]')
    plt.legend(loc='best')

    plt.figure(1)
    plt.xlabel('Strain [-]')
    plt.ylabel(r'$\nu$ [-]')

    plt.figure(2)

    plt.xlabel('Strain [-]')
    plt.ylabel('$E$ [GPa]')
    plt.plot(label='Exp')


    e0= plot_mechanical_data_for_simulation(directory)[0]
    E= plot_mechanical_data_for_simulation(directory)[1]
    plt.plot(e0, E, '--kx', ms= 12, mew=2, label='DEM data $b_{t}/R=0.05$')


    e0_2= plot_mechanical_data_for_simulation2(directory2)[0]
    E_2= plot_mechanical_data_for_simulation2(directory2)[1]
    plt.plot(e0_2, E_2, '--bx', ms= 12, mew=2, label='DEM data $b_{t}/R=0.1$')

    e01= plot_mechanical_data_for_simulation1(directory1)[0]
    E1= plot_mechanical_data_for_simulation1(directory1)[1]
    plt.plot(e01, E1, '-rx',ms= 12, mew=2 ,label='DEM data $b_{t}/R=0.2$')
    plt.legend(loc='best', numpoints=1)


    plt.show()


if __name__ == '__main__':
    main()