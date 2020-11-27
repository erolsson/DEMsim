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
    directory = os.path.expanduser('~/DEMsim/results/battery_rve/electrode_elaheh/contacts/')
    loaded = np.genfromtxt(directory + 'contacts_62.931.dou', delimiter=',')
    unloaded = np.genfromtxt(directory + 'contacts_65.743.dou', delimiter=', ')
    loaded_contacts = {}
    for c in loaded:
        if c[6] != 0.:
            loaded_contacts[(int(c[0]), int(c[1]))] = [c[5], c[6]]

    unloaded_contacts = {}
    for c in unloaded:
        if c[6] != 0.:
            unloaded_contacts[(int(c[0]), int(c[1]))] = [c[5], c[6]]

    force_difference = {}
    overlap_difference = {}
    for contact_pair, data in loaded_contacts.items():
        force_difference[contact_pair] = data[1]
        overlap_difference[contact_pair] = data[0]
        if contact_pair in unloaded_contacts:
            overlap_difference[contact_pair] -= unloaded_contacts[contact_pair][0]
            force_difference[contact_pair] -= unloaded_contacts[contact_pair][1]
    force_difference = np.sort(np.array(list(force_difference.values())))
    overlap_difference = np.sort(np.array(list(overlap_difference.values())))
    plt.figure(1)
    plt.plot(overlap_difference)
    plt.figure(2)
    plt.plot(force_difference)
    print(np.sum(force_difference))
    plt.show()


if __name__ == '__main__':
    main()
