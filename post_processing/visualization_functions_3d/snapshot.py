import pathlib
from collections import defaultdict
import os

from mayavi import mlab
import numpy as np

from visualization_functions_3d.plotting_functions import SpheresPlotter, SurfacesPlotter, BoundingBox
from visualization_functions_3d.periodic_bc import PeriodicBC
from visualization_functions_3d.battery_contact_plotter import BatteryContactPlotter
from visualization_functions_3d import colors


class Snapshot:
    def __init__(self, directory, contact_plotter_class=None):
        directory = pathlib.Path(directory)
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
        if os.path.isfile(self.directory / 'surface_positions.dou'):
            self.surfaces_plotter = SurfacesPlotter(self.directory / 'surface_positions.dou', self.surfaces_colors,
                                                    self.surfaces_opacities, self.plot_order, self.bounding_boxes,
                                                    self.visible_functions)
        if contact_plotter_class is not None:
            self.contact_plotter = contact_plotter_class(directory)
        else:
            self.contact_plotter = None

    def create_periodic_bc_plotter(self):
        self.periodic_bc_plotter = PeriodicBC(self.directory / 'periodic_bc.dou')

    def plot(self, time):
        if time.is_integer():
            time = int(time)

        particle_data = np.genfromtxt(self.directory / ('particles/particles_' + str(time) + '.dou'), delimiter=',')
        self.spheres_plotter.plot(particle_data)
        if self.mirror_particles:
            mirror_particle_data = np.genfromtxt(self.directory / ('mirror_particles/mirror_particles_'
                                                 + str(time) + '.dou'), delimiter=',')
            self.mirror_particles_plotter.plot(mirror_particle_data)
        if self.surfaces_plotter:
            self.surfaces_plotter.bounding_boxes = self.bounding_boxes
            self.surfaces_plotter.surfaces_colors = self.surfaces_colors
            self.surfaces_plotter.plot(time)
        if self.plot_periodic_bc:
            self.periodic_bc_plotter.plot(time)
        if self.contact_plotter:
            self.contact_plotter.plot(time)

        f = mlab.gcf()
        f.scene.render()


def main():
    mlab.figure(size=(1024, 768), bgcolor=(1., 1., 1.), fgcolor=(0, 0., 0.))
    snapshot = Snapshot(os.path.expanduser('~/DEMsim/results/asphalt_shear_box/bonded/mu=0.8_mu_wall=0.0/big_big_400kPa'))
    snapshot.plot_periodic_bc = False
    bbox = BoundingBox()
    bbox.z_max = lambda t: 0.04
    snapshot.bounding_boxes[3] = bbox
    snapshot.surfaces_colors[1] = colors.red
    snapshot.surfaces_colors[2] = colors.red
    snapshot.plot(2.31)
    mlab.view(10, -60, distance=0.25)
    mlab.savefig("shear.png")
    mlab.show()


if __name__ == '__main__':
    main()
