import glob
import os
import re
import time

import numpy as np
from mayavi import mlab

from visualization_functions_3d import SpheresPlotter
from visualization_functions_3d import SurfacesPlotter


class Animation:
    def __init__(self, directory):
        self.directory = directory
        self.delay = 0.0
        self.start_time = 0,
        self.end_time = None
        self.save_frames = False
        self.save_directory = None
        self.image_file_prefix = 'frame'
        self.image_file_extension = '.png'
        self.figure_directory = None

    def run(self):
        a = self._animation()
        for _ in a:
            pass

    def _animation(self):
        particle_files = glob.glob(self.directory + '/particles_*.dat')
        particle_files = [os.path.basename(particle_file) for particle_file in particle_files]
        frame_times = [float(re.findall('\d+.\d+', particle_file)[0]) for particle_file in particle_files]
        frame_times = np.array(sorted(frame_times))
        frame_times = frame_times[frame_times >= self.start_time]
        if self.end_time:
            frame_times = frame_times[frame_times < self.end_time]

        spheres_plotter = SpheresPlotter()
        surfaces_plotter = SurfacesPlotter(self.directory + '/surface_positions.dat')
        n = len(frame_times)
        for i, t in enumerate(frame_times):
            print t
            particle_data = np.genfromtxt(self.directory + 'particles_' + str(t) + '.dat', delimiter=',')
            spheres_plotter.plot(particle_data)
            surfaces_plotter.plot(t)
            if self.save_frames:
                if self.figure_directory is None:
                    self.figure_directory = self.directory
                if not os.path.isdir(self.figure_directory):
                    os.makedirs(self.figure_directory)
                # each file has name frame_00x, _0xx, xxx etc

                name = '/' + self.image_file_prefix + '0'*(len(str(n))-len(str(i+1))) + str(i+1) + '.' \
                       + self.image_file_extension
                mlab.savefig(filename=self.figure_directory + name)
            time.sleep(self.delay)
            yield