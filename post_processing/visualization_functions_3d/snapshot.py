from collections import defaultdict
import os

import numpy as np
from mayavi import mlab

from visualization_functions_3d.plotting_functions import SpheresPlotter, SurfacesPlotter, BoundingBox
from visualization_functions_3d.periodic_bc import PeriodicBC
from visualization_functions_3d.battery_contact_plotter import BatteryContactPlotter
from visualization_functions_3d import colors


class Snapshot:
    def __init__(self, directory, contact_plotter_class=None):
        self.plot_periodic_bc = True
        self.periodic_bc_plotter = None
        self.mirror_particles = False
        self.bounding_boxes = defaultdict(BoundingBox)

        self.directory = directory
        self.surfaces_colors = defaultdict(lambda: colors.blue)
        self.surfaces_opacities = defaultdict(lambda: 0.5)
        self.plot_order = None
        self.visible_functions = defaultdict(lambda: lambda t: True)

        self.spheres_plotter = SpheresPlotter()
        self.mirror_particles_plotter = SpheresPlotter(color=colors.silver)
        self.surfaces_plotter = None
        if os.path.isfile(self.directory + '/surface_positions.dou'):
            self.surfaces_plotter = SurfacesPlotter(self.directory + '/surface_positions.dou', self.surfaces_colors,
                                                    self.surfaces_opacities, self.plot_order, self.bounding_boxes,
                                                    self.visible_functions)
        if contact_plotter_class is not None:
            self.contact_plotter = contact_plotter_class(directory)
        else:
            self.contact_plotter = None

    def plot(self, time):
        if self.plot_periodic_bc:
            self.periodic_bc_plotter = PeriodicBC(self.directory + '/periodic_bc.dou')

        particle_data = np.genfromtxt(self.directory + '/particles/particles_' + str(time) + '.dou', delimiter=',')
        self.spheres_plotter.plot(particle_data)
        if self.mirror_particles:
            mirror_particle_data = np.genfromtxt(self.directory + '/mirror_particles/mirror_particles_'
                                                 + str(time) + '.dou', delimiter=',')
            self.mirror_particles_plotter.plot(mirror_particle_data)
        if self.surfaces_plotter:
            self.surfaces_plotter.plot(time)
        if self.plot_periodic_bc:
            self.periodic_bc_plotter.plot(time)
        if self.contact_plotter:
            self.contact_plotter.plot(time)

        f = mlab.gcf()
        f.scene.render()


def main():
    snapshot = Snapshot('C:/DEMsim/results/viscoelastic',
                        BatteryContactPlotter)
    snapshot.mirror_particles = True
    snapshot.contact_plotter.color = colors.red
    snapshot.contact_plotter.binder_radius = np.sqrt(0.3*0.01**2/np.pi)
    snapshot.plot(134.69)
    mlab.show()


if __name__ == '__main__':
    main()
