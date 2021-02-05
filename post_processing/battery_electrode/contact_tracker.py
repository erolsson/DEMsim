import glob
import os
import re

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


def get_contact_data(time, directory, contact_pair):
    t = time
    if t.is_integer():
        t = int(t)
    contact_data = np.genfromtxt(directory + '/contacts_' + str(t) + '.dou', delimiter=',')
    n = contact_data.shape[1]
    id1 = contact_data[:, 0]
    id2 = contact_data[:, 1]
    h_max = np.max(contact_data[:, 5])
    print(h_max, t)
    idx = np.logical_or(np.logical_and(id1 == contact_pair[0], id2 == contact_pair[1]),
                        np.logical_and(id1 == contact_pair[1], id2 == contact_pair[0]))
    data = np.zeros(n)
    data[1] = time
    if idx.any():
        data[0] = 1
        data[2:] = contact_data[idx, 2:].flatten()
    return data


def main():
    directory = os.path.expanduser('~/DEMsim/results/elaheh/cubic_box_FT/')
    contact_files = glob.glob(directory + '/contacts/*.dou')
    times = []
    for contact_filename in contact_files:
        time = contact_filename[len(directory + '/contacts/contacts_'):-4]
        times.append(float(time))
    times = np.array(sorted(times))
    contact_data = np.array([get_contact_data(t, directory + '/contacts/', contact_pair=(22, 33)) for t in times])
    idx = contact_data[:, 0] == 1
    plt.figure(0)
    plt.plot(contact_data[idx, 1], contact_data[idx, 5], '-*')
    plt.figure(1)
    plt.plot(contact_data[idx, 5], contact_data[idx, 6], '-*')
    plt.figure(2)
    plt.plot(contact_data[idx, 1], contact_data[idx, 6], '-*')

    plt.show()


if __name__ == '__main__':
    main()
