from collections import defaultdict
import glob
import os
import re
import time

import numpy as np
from mayavi import mlab

from visualization_functions_3d.plotting_functions import SpheresPlotter
from visualization_functions_3d.plotting_functions import SurfacesPlotter
from visualization_functions_3d.plotting_functions import BoundingBox
from visualization_functions_3d.periodic_bc import PeriodicBC
from visualization_functions_3d import colors


class Animation:
    def __init__(self, directory):
        self.directory = directory
        self.delay = 0.0
        self.start_time = 0
        self.end_time = None
        self.save_frames = False
        self.frame_times = None
        self.save_directory = ''
        self.image_file_prefix = 'frame'
        self.image_file_extension = 'png'
        self.figure_directory = None
        self.surfaces_colors = defaultdict(lambda: (0., 0., 1.))
        self.surfaces_opacities = defaultdict(lambda: 0.5)
        self.plot_order = None
        self.bounding_boxes = defaultdict(BoundingBox)
        self.visible_functions = defaultdict(lambda: lambda t: True)
        self.dpi = 500, 800
        self.zoom_settings = None
        self.plot_periodic_bc = False
        self.periodic_bc_plotter = None
        self.mirror_particles = False

        self.spheres_plotter = SpheresPlotter()
        self.mirror_particles_plotter = SpheresPlotter(color=colors.silver)

        self.initialized = False
        self.surfaces_plotter = None
        self.view_surfaces = True

    def run(self):
        a = self._animation()
        for _ in a:
            pass

    def initialize(self):
        if self.view_surfaces:
            self.surfaces_plotter = SurfacesPlotter(self.directory + '/surface_positions.dou', self.surfaces_colors,
                                                    self.surfaces_opacities, self.plot_order, self.bounding_boxes,
                                                    self.visible_functions)
            self.surfaces_plotter.surfaces_opacities = self.surfaces_opacities
            self.surfaces_plotter.surfaces_colors = self.surfaces_colors
            self.surfaces_plotter.plot_order = self.plot_order
            self.surfaces_plotter.bounding_boxes = self.bounding_boxes
            self.surfaces_plotter.visible_times = self.visible_functions

        if self.plot_periodic_bc:
            self.periodic_bc_plotter = PeriodicBC(self.directory + '/periodic_bc.dou')

        particle_files = glob.glob(self.directory + '/particles/particles_*.dou')
        particle_files = [os.path.basename(particle_file) for particle_file in particle_files]

        self.frame_times = []
        for p_file in particle_files:
            self.frame_times.append(re.findall(r"[-+]?\d*\.\d+|\d+", p_file)[0])
        self.frame_times = np.array(sorted(self.frame_times, key=lambda x: float(x)), dtype=str)
        frame_times_np = np.array([float(t) for t in self.frame_times])
        self.frame_times = self.frame_times[frame_times_np >= self.start_time]
        frame_times_np = frame_times_np[frame_times_np >= self.start_time]
        if self.end_time:
            self.frame_times = self.frame_times[frame_times_np < self.end_time]
        print(self.frame_times)

        if self.save_frames and not os.path.isdir(self.save_directory):
            os.makedirs(self.save_directory)
        self.initialized = True

    def _animation(self):
        if not self.initialized:
            self.initialize()
        n = len(self.frame_times)
        for i, t in enumerate(self.frame_times):
            print(i)
            particle_data = np.genfromtxt(self.directory + '/particles/particles_' + t + '.dou', delimiter=',')
            self.spheres_plotter.plot(particle_data)
            if self.mirror_particles:
                mirror_particle_data = np.genfromtxt(self.directory + '/mirror_particles/mirror_particles_'
                                                     + t + '.dou', delimiter=',')
                self.mirror_particles_plotter.plot(mirror_particle_data)
            if self.surfaces_plotter:
                self.surfaces_plotter.plot(float(t))
            if self.plot_periodic_bc:
                self.periodic_bc_plotter.plot(float(t))

            f = mlab.gcf()
            f.scene.render()

            if self.save_frames:
                if self.figure_directory is None:
                    self.figure_directory = self.directory
                if not os.path.isdir(self.figure_directory):
                    os.makedirs(self.figure_directory)

                # each file has name frame_00x, _0xx, xxx etc
                name = '/' + self.image_file_prefix + '0'*(len(str(n))-len(str(i))) + str(i) + '.' \
                       + self.image_file_extension
                mlab.savefig(filename=self.save_directory + name)

            time.sleep(self.delay)
            yield
