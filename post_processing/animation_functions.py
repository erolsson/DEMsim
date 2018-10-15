import glob
import os
import re

import numpy as np
from mayavi import mlab

from visualization_functions_3d import SpheresPlotter
from visualization_functions_3d import SurfacesPlotter


def animate_simulation(directory, delay=10, start_time=0., end_time=None, save_frames=False, save_directory=None,
                       image_file_prefix='frame', image_file_extension='png'):
    @mlab.animate(delay=delay)
    def animation(figure_directory, file_prefix, file_ext):
        particle_files = glob.glob(directory + '/particles_*.dat')
        particle_files = [os.path.basename(particle_file) for particle_file in particle_files]
        frame_times = [float(re.findall('\d+.\d+', particle_file)[0]) for particle_file in particle_files]
        frame_times = np.array(sorted(frame_times))
        frame_times = frame_times[frame_times >= start_time]
        if end_time:
            frame_times = frame_times[frame_times<end_time]

        spheres_plotter = SpheresPlotter()
        surfaces_plotter = SurfacesPlotter(directory + '/surface_positions.dat')
        f = mlab.gcf()
        n = len(frame_times)
        for i, t in enumerate(frame_times):
            particle_data = np.genfromtxt(directory + 'particles_' + str(t) + '.dat', delimiter=',')
            spheres_plotter.plot(particle_data)
            surfaces_plotter.plot()
            f.scene.render()
            mlab.show()
            if save_frames:
                if figure_directory is None:
                    figure_directory = directory
                if not os.path.isdir(figure_directory):
                    os.makedirs(figure_directory)
                # each file has name frame_00x, _0xx, xxx etc
                name = '/' + file_prefix + '0'*(len(str(n))-len(str(i))) + str(i) + '.' + file_ext
                mlab.savefig(filename=figure_directory + name)
            yield

    animation(save_directory, image_file_prefix, image_file_extension)
