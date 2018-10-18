import sys

from collections import namedtuple
from collections import OrderedDict
from math import pi

import numpy as np

from mayavi import mlab


class BoundingBox:
    def __init__(self):
        self.x_min = -1e99,
        self.x_max = 1e99
        self.y_min = -1e99,
        self.y_max = 1e99,
        self.z_min = -1e99
        self.z_max = 1e99


# Todo, fix unused parameter opacity
class SpheresPlotter:
    def __init__(self, opacity=1., color=(184./255, 115./255., 51./255.)):
        self.ms = None
        self.color = color
        self.color = color
        del opacity

    def plot(self, data, ):
        x = data[:, 1]
        y = data[:, 2]
        z = data[:, 3]
        r = data[:, 7]
        if self.ms is None:
            self.ms = mlab.points3d(x, y, z, 2 * r,
                                    color=self.color,
                                    resolution=32,
                                    scale_factor=1.,
                                    scale_mode='scalar').mlab_source
        else:
            self.ms.set(x=x, y=y, z=z)


def fulfill_bounding_box(bounding_box, x, y, z):
    x[x < bounding_box.x_min] = bounding_box.x_min
    x[x > bounding_box.x_max] = bounding_box.x_max
    y[y < bounding_box.y_min] = bounding_box.y_min
    y[y > bounding_box.y_max] = bounding_box.y_max
    z[z < bounding_box.z_min] = bounding_box.z_min
    z[z > bounding_box.z_max] = bounding_box.z_max


class PointSurfacePlotter:
    def __init__(self, opacity=0.5, color=(0., 0., 1.)):
        self.ms = None
        self.opacity = opacity
        self.color = color
        self.bounding_box = BoundingBox()

    def plot(self, data):
        n = data.shape[0]/3

        x = data[0:3*n-2:3]
        y = data[1:3*n-1:3]
        z = data[2:3*n:3]
        fulfill_bounding_box(self.bounding_box, x, y, z)

        pts = mlab.points3d(x, y, z, z)
        mesh = mlab.pipeline.delaunay2d(pts)

        pts.remove()

        if self.ms is None:
            self.ms = mlab.pipeline.surface(mesh, color=self.color, opacity=self.opacity, transparent=True).mlab_source
        else:
            self.ms.set(x=x, y=y, z=z)


class CylinderPlotter:
    def __init__(self, opacity=0.5, color=(0., 0., 1)):
        self.ms = None
        self.opacity = opacity
        self.color = color
        self.bounding_box = BoundingBox()

    def plot(self, data):
        r = data[0]
        # axis = data[1:4]   # Todo use axis parameter
        point = data[4:7]
        length = data[7]

        q, z = np.meshgrid(np.linspace(0, 2*pi, 100), np.linspace(point[2], point[2] + length, 100))

        x = r*np.cos(q) + point[0]
        y = r*np.sin(q) + point[1]
        fulfill_bounding_box(self.bounding_box, x, y, z)
        if self.ms is None:
            self.ms = mlab.mesh(x, y, z, color=self.color, opacity=self.opacity, transparent=True).mlab_source
        else:
            self.ms.set(x=x, y=y, z=z)


PlotObject_ = namedtuple('PlotObject', ['start_idx', 'end_idx'])


class SurfacesPlotter:
    def __init__(self, surface_file_name):
        self.plotters = {}
        self.plotter_data = {}
        self.data = OrderedDict()
        self.counter = 0
        self.set_data_file(surface_file_name)

    def set_data_file(self, surface_file_name):
        with open(surface_file_name) as data_file:
            data_lines = data_file.readlines()

        # Inspect the first line
        line = data_lines[0]
        words = line.split(", ")

        id_idx = [i for i in range(len(words)) if words[i].startswith('ID')]
        for idx in id_idx:
            surface_id = int(words[idx][3:])
            surface_type = words[idx+1][5:]
            if surface_type == 'Cylinder':
                self.plotter_data[surface_id] = PlotObject_(idx+2, idx+10)
                self.plotters[surface_id] = CylinderPlotter()
            elif surface_type == 'PointSurface':
                num_points = int(words[idx+2])
                self.plotter_data[surface_id] = PlotObject_(idx+3, idx+3+num_points*3)
                self.plotters[surface_id] = PointSurfacePlotter()

        data = np.genfromtxt(surface_file_name, delimiter=',')
        for i in range(data.shape[0]):
            self.data[data[i, -1]] = data[i, :-1]

    def plot(self, t=None):
        if t is None:
            data_line = self.data[self.data.keys()[self.counter]]
            self.counter += 1
        else:
            try:
                data_line = self.data[t]
                self.counter = 0
            except LookupError:
                print "No surface data at time ", t
                sys.exit(1)

        for surface_id, plotter in self.plotters.iteritems():
            plotter_data = self.plotter_data[surface_id]
            data = data_line[plotter_data.start_idx:plotter_data.end_idx]
            plotter.plot(data)

    def set_bounding_box(self, bounding_box):
        for plotter in self.plotters.itervalues():
            plotter.bounding_box = bounding_box


if __name__ == '__main__':
    simulation_directory = '../results/gyratory_compaction/1/'
    surfaces_plotter = SurfacesPlotter(simulation_directory + 'surface_positions.dat')
    bbox = BoundingBox()
    bbox.z_min = -0.01
    bbox.z_max = 0.05
    surfaces_plotter.set_bounding_box(bbox)
    surfaces_plotter.plot()

    # mlab.show()
