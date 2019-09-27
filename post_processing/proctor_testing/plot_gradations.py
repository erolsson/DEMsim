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
labels = ['8-16 mm', r'F{\"u}ller curve']
bg_index = {'8-16 mm': {}, r'F{\"u}ller curve': {}}
for i, sim in enumerate(simulations):
    for material, line, in zip(['', '_weak'], ['--', ':']):
        simulation = sim + material

        with open(base_directory + '/' + simulation + '/gradation_pickle.pkl', 'r') as pickle_handle:
            r0 = pickle.load(pickle_handle)
            r1 = pickle.load(pickle_handle)

        if material == '':
            particle_volume = 4*np.pi*r0**3/3
            volume_cdf = np.cumsum(particle_volume)/np.sum(particle_volume)
            plt.semilogx(2*r0, volume_cdf*100, '' + colors[i], lw=2, label=labels[i])

        particle_volume = 4*np.pi*r1**3/3
        volume_cdf = np.cumsum(particle_volume)/np.sum(particle_volume)
        plt.semilogx(2*r1, volume_cdf*100, line + colors[i], lw=2)

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
            if diff > 0:
                bg += diff
        print "Bg for", sim, material, "is", bg*100
        bg_index[labels[i]][material] = bg*100

print bg_index
for line, label in zip(['w', '-k', '--k', ':k'], [r'$\quad$', 'Starting',
                                                  r'$\sigma_w=387.5$ MPa', r'$\sigma_w=200$ MPa']):
    plt.plot([1, 2], [-1, -1], line, lw=2, label=label)
plt.xlim(1, 16)
plt.ylim(0, 100)
sieves = [1., 2., 4., 8., 11.2, 16.]
plt.xticks(sieves, [str(s) for s in sieves])
plt.xlabel('Size [mm]')
plt.ylabel(r'Weight passing [\%]')
ax = plt.gca()
box = ax.get_position()
plt.text(0.05, 0.9, '(a)', horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
plt.legend(loc='upper left', bbox_to_anchor=(0., 0.89), framealpha=0.9)
plt.tight_layout()
plt.savefig('gradations_after_test.png')

plt.figure(2)
plt.bar(0, bg_index['8-16 mm'][''], 1, color='b', label=labels[0])
plt.bar(1, bg_index[r'F{\"u}ller curve'][''], 1, color='r', label=labels[1])
print bg_index['8-16 mm']['']/bg_index[r'F{\"u}ller curve']['']
plt.bar(3, bg_index['8-16 mm']['_weak'], 1, color='b')
plt.bar(4, bg_index[r'F{\"u}ller curve']['_weak'], 1, color='r')
print bg_index['8-16 mm']['_weak']/bg_index[r'F{\"u}ller curve']['_weak']

plt.legend(loc='upper left', bbox_to_anchor=(0., 0.89), framealpha=0.9)

plt.xticks([1, 4], [r'$\sigma_w=387.5$ MPa', r'$\sigma_w=200$ MPa'])
plt.ylabel('Modified Bg. Index')

ax = plt.gca()
ax.tick_params('x', pad=15)
plt.tight_layout()
plt.text(0.05, 0.9, '(b)', horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
plt.savefig('sim_bg.png')
plt.show()
