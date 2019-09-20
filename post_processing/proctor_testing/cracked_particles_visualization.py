from collections import namedtuple
import os

import numpy as np

from mayavi import mlab

from post_processing.visualization_functions_3d import SpheresPlotter

SimulationData = namedtuple('SimulationData', ['particles', 'fracture_data'])
simulations = []
layer = 4
stroke = 24
mlab.figure(size=(1920, 1200), bgcolor=(1., 1., 1.))
max_cracks = []
simulation_name = '8-16mm_continued'
for d, post_fix in zip([-0.06, 0.06], ['', '_weak']):
    base_directory = os.path.expanduser('~/DEMsim/results/proctor_test/' + simulation_name + post_fix)
    fracture_file = base_directory + '/layer_' + str(layer) + '/stroke_' + str(stroke) + '/particle_cracks.dou'
    time_data = np.genfromtxt(fracture_file, delimiter=',', dtype=None, usecols=(-1,))
    time = time_data[-1]
    full_data = np.genfromtxt(fracture_file, delimiter=',', dtype=None)

    data = full_data[time_data == time]
    particle_time = np.round(time, 2)

    particle_file = base_directory + '/animation/particles_' + str(particle_time) + '.dou'
    particle_data = np.genfromtxt(particle_file, delimiter=',')
    particle_data[:, 1] += d
    fracture_count = 0*particle_data[:, 0]
    for data_line in data:
        fracture_count[int(data_line[0][3:])] += 1
    sim = SimulationData(particles=particle_data,
                         fracture_data=fracture_count)
    simulations.append(sim)

max_cracks = max([np.max(sim.fracture_data) for sim in simulations])

for sim in simulations:
    plotter_cracked = SpheresPlotter()
    plotter_cracked.plot(sim.particles[sim.fracture_data > 0, :],
                         function=0*sim.fracture_data[sim.fracture_data > 0] +1, vmin=0, vmax=1)

    plotter_uncracked = SpheresPlotter(opacity=0.6, color=(128./255, 128./255, 128./255))
    plotter_uncracked.plot(sim.particles[sim.fracture_data == 0, :])

view_data = mlab.view()
mlab.view(azimuth=-90, elevation=view_data[1], distance=0.5, focalpoint=(0, 0, 0.075))

# bar = mlab.scalarbar(nb_labels=int(max_cracks), label_fmt='%.0f')
# bar.label_text_property.color = (0., 0., 0.)
# bar.scalar_bar_representation.position = [0.25, 0.05]
# bar.scalar_bar_representation.position2 = [0.5, 0.05]
mlab.savefig('fractured_particles_8-16mm.png')
mlab.show()
