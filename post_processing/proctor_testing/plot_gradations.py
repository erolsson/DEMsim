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
simulations = ['8-16mm_continued', 'fuller']
colors = ['b', 'r']
for i, sim in enumerate(simulations):
    plt.figure(i)
    for material, line, label in zip(['', '_weak'], ['-', '--'], [r'$\sigma_w=387.5$ MPa', r'$\sigma_w=200$ MPa']):
        simulation = sim + material

        with open(base_directory + '/' + simulation + '/gradation_pickle.pkl', 'r') as pickle_handle:
            r0 = pickle.load(pickle_handle)
            r1 = pickle.load(pickle_handle)

        if material == '':
            particle_volume = 4*np.pi*r0**3/3
            volume_cdf = np.cumsum(particle_volume)/np.sum(particle_volume)
            plt.plot(2*r0, volume_cdf, ':' + colors[i], lw=3, label='Starting')

        particle_volume = 4*np.pi*r1**3/3
        volume_cdf = np.cumsum(particle_volume)/np.sum(particle_volume)
        plt.plot(2*r1, volume_cdf, line + colors[i], lw=2, label=label)

        cdf_start = np.zeros(r0.shape[0] + 1)
        particle_volume_start = 4*np.pi*r0**3/3
        cdf_start[1:] = np.cumsum(particle_volume_start)/np.sum(particle_volume_start)
        start_size = np.zeros(r0.shape[0] + 1)
        start_size[1:] = 2*r0
        start_size[0] = 0.0
        end_size = np.zeros(r1.shape[0] + 1)
        end_size[1:] = 2*r1
        end_size[0] = 0.0
        cdf_end = np.zeros(r1.shape[0] + 1)

        particle_volume_end = 4*np.pi*r1**3/3
        cdf_end[1:] = np.cumsum(particle_volume_end)/np.sum(particle_volume_end)

        sieve_sizes = [0.063, 0.125, 0.25, 0.5, 1., 2, 4, 8, 11.2, 16]
        bg = 0
        for j in range(1, len(sieve_sizes)):
            dv_start = (np.interp(sieve_sizes[j], start_size, cdf_start) -
                        np.interp(sieve_sizes[j-1], start_size, cdf_start))
            dv_end = (np.interp(sieve_sizes[j], end_size, cdf_end) -
                      np.interp(sieve_sizes[j-1], end_size, cdf_end))
            diff = dv_end - dv_start
            print diff
            if diff > 0:
                bg += diff
        print "Bg for", sim, material, "is", bg*100
    plt.xlabel('Size [mm]')
    plt.ylabel('Weight passing')
    plt.legend(loc='best')
plt.figure(0)
plt.savefig('gradations_8-16mm.png')

plt.figure(1)
plt.savefig('gradations_fuller.png')

plt.show()
