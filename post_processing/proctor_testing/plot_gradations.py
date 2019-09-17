import os
import pickle

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

base_directory = os.path.expanduser('~/DEMsim/results/proctor_test/')
simulation = '8-16mm'

with open(base_directory + '/' + simulation + '/gradation_pickle.pkl', 'r') as pickle_handle:
    r0 = pickle.load(pickle_handle)
    r1 = pickle.load(pickle_handle)

for r in [r0, r1]:
    number_cdf = np.linspace(0, 1, r.shape[0], endpoint=False)
    plt.figure(0)
    plt.plot(2*r, number_cdf)
    plt.xlabel('Size [mm]')
    plt.ylabel('Number passing')

    c = r[-1]**3 - 3*np.trapz(number_cdf*r**2, r)
    print c
    volume_cdf = number_cdf*r**3
    print volume_cdf
    for i in range(volume_cdf.shape[0]):
        volume_cdf[i] -= 3*np.trapz(number_cdf[r <= r[i]]*r[r <= r[i]]**2, r[r <= r[i]])
    volume_cdf /= c
    plt.figure(1)
    plt.plot(2*r, volume_cdf)
    plt.xlabel('Size [mm]')
    plt.ylabel('Weight passing')

plt.show()
