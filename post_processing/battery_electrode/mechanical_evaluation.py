import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def plot_mechanical_data_for_simulation(directory):
    periodic_bc = np.genfromtxt(directory + 'periodic_bc.dou', delimiter=',')
    surface_forces = np.genfromtxt(directory + 'surface_forces.dou', delimiter=',')
    surface_positions = np.genfromtxt(directory + 'surface_positions.dou', delimiter=',')
    force_fabric_tensor = np.genfromtxt(directory + 'force_fabric_tensor.dou', delimiter=',')
    time = surface_positions[:, -1]
    t = surface_positions[:, -2]
    d = 2*periodic_bc[:, 2]
    d0 = d[0]
    w = 2*periodic_bc[:, 4]
    compaction_force = surface_forces[:, -2]

    # Finding the point where the compaction force has decreased to zero after its maximum value
    # The thickness at this point will be the thickness of the electrode t0
    idx = np.logical_and(compaction_force == 0, time > time[np.argmax(compaction_force)])
    t0 = 0.882022

    t_start = time[d == d0][-1]
    volume = (d*w*t0)[time > t_start]
    sxx = force_fabric_tensor[time > t_start, 1]/volume
    sxx -= sxx[0]

    syy = force_fabric_tensor[time > t_start, 5]/volume
    syy -= syy[0]

    szz = force_fabric_tensor[time > t_start, 9]/volume
    szz -= szz[0]
    sxx *= -1
    syy *= -1
    szz *= -1

    linear_strain = (d[time > t_start] - d0)/d0
    time_mechanical = time[time > t_start]
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
            if dsxx != 0:

                v = dsyy/(dsxx + dszz)
                plt.figure(1)
                plt.plot(e0, v, 'x', ms=12, mew=2)

                E = (dsxx - v*(dsyy + dszz))/dexx
                plt.figure(2)
                plt.plot(e0, E/1e9, 'x', ms=12, mew=2)


def main():
    directory = os.path.expanduser(r'C:/DEMsim/results/relaxaiton/')
    plot_mechanical_data_for_simulation(directory)

    directory = os.path.expanduser(r'C:/DEMsim/results/relaxaiton/')
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

    plt.show()


if __name__ == '__main__':
    main()
