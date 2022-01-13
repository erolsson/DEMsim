import pathlib
import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 24})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{xcolor}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def calculate_mechanical_data(directory):
    periodic_bc = np.genfromtxt(directory / 'periodic_bc.dou', delimiter=',')
    surface_positions = np.genfromtxt(directory / 'surface_positions.dou', delimiter=',')
    force_fabric_tensor = np.genfromtxt(directory / 'force_fabric_tensor.dou', delimiter=',')
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

    e0_vec = []
    E_vec = []
    v_vec = []

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

                E = (dsxx - v*(dsyy + dszz))/dexx
                e0_vec.append(e0)

                E_vec.append(E/1e9)
                v_vec.append(v)

    return np.array([e0_vec, E_vec, v_vec]).T


"""
Experimental data
e0_exp_comp = np.array([-0.01,-0.011,-0.0123,-0.0141,-0.0165])
E_exp_comp = np.array([1.20, 1.43, 1.51, 1.55, 1.99])
plt.plot(e0_exp_comp, E_exp_comp, 'bo', ms=12)
e0_exp_ten = np.array([0.01,0.011,0.0123,0.0141,0.0165])
E_exp_ten = np.array([0.95, 0.78, 0.84, 0.76, 1.06])
plt.plot(e0_exp_ten, E_exp_ten, 'bo', ms=12 )
plt.plot(e0_exp_ten, E_exp_ten, 'bo', ms=12, label=r'Experiment data')
"""


def main():
    for bt in [0.05, 0.1, 0.2]:
        data = []
        for direction in ["tension", "compression"]:
            directory = pathlib.Path("~/batteries/DEM_elaheh") / (direction + "-E34bt" + str(bt).replace('.', '')
                                                                  + "Rbr05")
            directory = directory.expanduser()
            data.extend(calculate_mechanical_data(directory))
        data = np.array(data)
        data = data[np.argsort(data[:, 0]), :]
        plt.plot(data[:, 0]*100, data[:, 1], '--x', lw=3, ms=12, mew=3, label="$b_t=" + str(bt) + "R$")
    e0_exp_comp = np.array([-0.01, -0.011, -0.0123, -0.0141, -0.0165])
    E_exp_comp = np.array([1.20, 1.43, 1.51, 1.55, 1.99])
    plt.plot(e0_exp_comp*100, E_exp_comp, 'kx', ms=12, mew=3)
    e0_exp_ten = np.array([0.01, 0.011, 0.0123, 0.0141, 0.0165])
    E_exp_ten = np.array([0.95, 0.78, 0.84, 0.76, 1.06])
    plt.plot(e0_exp_ten*100, E_exp_ten, 'kx', ms=12, mew=3, label=r'Exp.')

    plt.xlabel(r'Strain [\%]')
    plt.ylabel('$E$ [GPa]')
    plt.legend(numpoints=1)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
