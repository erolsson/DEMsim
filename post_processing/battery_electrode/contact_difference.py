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


def main():
    directory = os.path.expanduser('~/DEMsim/results/battery_rve/Cathod/contacts/')
    file_1 = np.genfromtxt(directory + 'contacts_24.254.dou', delimiter=',')
    file_2 = np.genfromtxt(directory + 'contacts_24.257.dou', delimiter=', ')
    data_1 = {}
    data_2 = {}
    for data_file, data_dict in zip([file_1, file_2], [data_1, data_2]):
        for c in data_file:
            data_dict[(c[0], c[1])] = c[6]

    difference = {}
    for c_id, f2 in data_2.items():
        if f2 != 0:
            print(f2)
            difference[c_id] = f2 - data_1[c_id]
    print(difference)


if __name__ == '__main__':
    main()
