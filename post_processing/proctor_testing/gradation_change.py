from collections import defaultdict

import pickle
import sys

import numpy as np

base_directory = sys.argv[-1]
particle_data = np.genfromtxt(base_directory + "/animation/particles_84.59.dou", delimiter=",")
r0 = particle_data[:, -3]*1000

particle_dict = dict(zip(particle_data[:, 0], r0))

r0 = np.sort(r0)

filename = base_directory + '/layer_4/stroke_24/particle_cracks.dou'
full_data = np.genfromtxt(filename, delimiter=',', dtype=None, usecols=(0, 1, 6, 7, 8))

time_data = np.genfromtxt(filename, delimiter=',', dtype=None, usecols=(-1,))
data = full_data[time_data == time_data[-1]]
data = np.unique(data)

cracked_particles = defaultdict(int)
for crack in data:
    particle_id = int(crack[0][3:])
    cracked_particles[particle_id] += 1
for i, (particle_id, cracks) in enumerate(cracked_particles.iteritems()):
    r = particle_dict[particle_id]
    vol = 4*np.pi*r**3/3
    particle_dict.pop(particle_id)
    for j in range(cracks+1):
        particle_dict[(i+1)*100000 + j] = (3*vol/(cracks + 1)/4/np.pi)**(1./3)

r1 = np.array(sorted(particle_dict.values()))

with open(base_directory + "/gradation_pickle.pkl", 'w') as pickle_handle:
    pickle.dump(r0, pickle_handle)
    pickle.dump(r1, pickle_handle)
