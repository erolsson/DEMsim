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


def calculate_stress_tensor(directory):
    force_tensor = np.genfromtxt(directory + '/force_fabric_tensor.dou', delimiter=',')
    periodic_box = np.genfromtxt(directory + '/periodic_bc.dou', delimiter=',')
    volume = 1 + 0*periodic_box[:, 0]
    time = force_tensor[:, 0]
    for i in range(3):
        volume *= (periodic_box[:, 2*i+2] - periodic_box[:, 2*i+1])

    return time, force_tensor[:, 1:]/volume[:, None]


def main():
    time, s = calculate_stress_tensor(os.path.expanduser('~/DEMsim/results/periodic_bc_test/sim_1'))
    plt.plot(time, s)
    time, s = calculate_stress_tensor(os.path.expanduser('~/DEMsim/results/periodic_bc_test/sim_1/original'))
    plt.plot(time, s)

    time, s = calculate_stress_tensor(os.path.expanduser('~/DEMsim/results/periodic_bc_test/sim_1/restart'))
    plt.plot(time, s, '--')
    plt.show()


if __name__ == '__main__':
    main()
