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
                plt.plot(e0, v, 'kx', ms=12, mew=2)

                E = (dsxx - v*(dsyy + dszz))/dexx
                plt.figure(2)
                plt.plot(e0, E/1e9, 'kx', ms=12, mew=2)

                e0_exp_comp = np.array([-0.01,-0.011,-0.0123,-0.0141,-0.0165])
                E_exp_comp = np.array([1.20, 1.43, 1.51, 1.55, 1.99])
                plt.plot(e0_exp_comp, E_exp_comp, 'bo', ms=12)
                e0_exp_ten = np.array([0.01,0.011,0.0123,0.0141,0.0165])
                E_exp_ten = np.array([0.95, 0.78, 0.84, 0.76, 1.06])
                plt.plot(e0_exp_ten, E_exp_ten, 'bo', ms=12)




def main():

    directory = os.path.expanduser(r'C:/DEMsim/results/Compression-E34br05bt01/')
    plot_mechanical_data_for_simulation(directory)


    directory = os.path.expanduser(r'C:/DEMsim/results/tension-E34bt01Rbr05/')
    plot_mechanical_data_for_simulation(directory)



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
    plt.legend(['DEM','EXP'], loc='upper right')

    plt.show()





if __name__ == '__main__':
    main()