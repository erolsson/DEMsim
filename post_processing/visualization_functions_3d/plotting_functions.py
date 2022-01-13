import sys

from collections import namedtuple
from collections import OrderedDict
from math import pi

import numpy as np

from mayavi import mlab

from visualization_functions_3d import colors

class BoundingBox:
    def __init__(self):
        self.x_min = lambda t: -1e99
        self.x_max = lambda t: 1e99
        self.y_min = lambda t: -1e99
        self.y_max = lambda t: 1e99
        self.z_min = lambda t: -1e99
        self.z_max = lambda t: 1e99

    def values(self, time):
        return [self.x_min(time), self.x_max(time),
                self.y_min(time), self.y_max(time),
                self.z_min(time), self.z_max(time)]


# Todo, fix unused parameter opacity
class SpheresPlotter:
    def __init__(self, opacity=1., color=colors.copper):
        self.ms = None
        self.color = color
        self.opacity = opacity

    def plot(self, data, function=None, vmin=None, vmax=None):
        if len(data) > 0:
            if len(data.shape) == 1:
                data = np.expand_dims(data, 0)
            x = data[:, 1]
            y = data[:, 2]
            z = data[:, 3]
            r = data[:, 7]
            if vmin is None:
                vmin = np.min(function)
            if vmax is None:
                vmax = np.max(function)
            if function is not None:
                pts = mlab.points3d(x, y, z, resolution=32, opacity=self.opacity, transparent=self.opacity != 1.,
                                    scale_factor=1., vmax=vmax, vmin=vmin)
                self.ms = pts.mlab_source
                pts.glyph.scale_mode = 'scale_by_vector'
                self.ms.dataset.point_data.vectors = np.tile(2/np.sqrt(3)*r, (3, 1)).transpose()
                self.ms.dataset.point_data.scalars = function

            else:
                if self.ms is None:
                    self.ms = mlab.points3d(x, y, z, 2*r,
                                            color=self.color,
                                            resolution=32,
                                            scale_factor=1.,
                                            scale_mode='scalar',
                                            opacity=self.opacity,
                                            reset_zoom=False).mlab_source
                elif self.ms.points.shape[0] == data.shape[0]:
                    self.ms.set(x=x, y=y, z=z)
                else:
                    self.ms.reset(x=x, y=y, z=z, scalars=2*r)


def fulfill_bounding_box(bounding_box, x, y, z, time):
    values = bounding_box.values(time)
    x[x < values[0]] = values[0]
    x[x > values[1]] = values[1]
    y[y < values[2]] = values[2]
    y[y > values[3]] = values[3]
    z[z < values[4]] = values[4]
    z[z > values[5]] = values[5]


class PointSurfacePlotter:
    def __init__(self, bounding_box=None):
        self.ms = None
        self.bounding_box = bounding_box

    def plot(self, data, color, opacity, time=0.):
        n = int(data.shape[0]/3)
        x = data[0:3*n-2:3]
        y = data[1:3*n-1:3]
        z = data[2:3*n:3]
        if self.bounding_box:
            fulfill_bounding_box(self.bounding_box, x, y, z, time)

        # This orders the points so that a rectange is plotted
        pts = mlab.points3d(x, y, z, z)
        mesh = mlab.pipeline.delaunay2d(pts)

        pts.remove()
        if self.ms is None:
            self.ms = mlab.pipeline.surface(mesh, color=color, opacity=opacity, transparent=True,
                                            reset_zoom=False).mlab_source
        else:
            # Updating the pipeline with the new set of points
            # There is probably a cuter way to do this
            self.ms.points = mesh.mlab_source.points


class CylinderPlotter:
    def __init__(self, bounding_box=None):
        self.ms = None
        self.upper_plate = None
        self.lower_plate = None
        self.bounding_box = bounding_box
        self.length_extension = 0.
        self.closed = False

    def plot(self, data, color, opacity, time=0.):
        r = data[0]
        # axis = data[1:4]   # Todo use axis parameter
        point = data[4:7]
        length = data[7] + self.length_extension

        q, z = np.meshgrid(np.linspace(0, 2*pi, 100), np.linspace(point[2], point[2] + length, 100))

        x = r*np.cos(q) + point[0]
        y = r*np.sin(q) + point[1]
        if self.bounding_box:
            fulfill_bounding_box(self.bounding_box, x, y, z, time)
        if self.ms is None:
            self.ms = mlab.mesh(x, y, z, color=color, opacity=opacity, transparent=True, reset_zoom=False).mlab_source
        else:
            self.ms.set(x=x, y=y, z=z)

        if self.closed:
            rad, q = np.meshgrid(np.linspace(0, r, 100), np.linspace(0, 2*pi, 100))
            x = rad*np.cos(q) + point[0]
            y = rad*np.sin(q) + point[1]

            if self.upper_plate is None:
                self.upper_plate = mlab.mesh(x, y, 0*x + point[2], color=color, opacity=opacity, transparent=True,
                                             reset_zoom=False).mlab_source
            else:
                self.upper_plate.set(x=x, y=y, z=z)

            if self.lower_plate is None:
                self.lower_plate = mlab.mesh(x, y, 0*x, color=color, opacity=opacity, transparent=True,
                                             reset_zoom=False).mlab_source
            else:
                self.lower_plate.set(x=x, y=y, z=z)


PlotObject_ = namedtuple('PlotObject', ['start_idx', 'end_idx'])


class SurfacesPlotter:
    def __init__(self, surface_file_name, surfaces_colors=None, surfaces_opacities=None, plot_order=None,
                 bounding_boxes=None, visible_times=None):

        self.plotters = {}
        self.plotter_data = {}
        self.data = OrderedDict()
        self.counter = 0
        self.set_data_file(surface_file_name)

        self.surfaces_colors = {}
        self.surfaces_opacities = {}
        self.plot_order = plot_order
        self.bounding_boxes = {}
        self.visible_times = {}

        if surfaces_colors is None:
            surfaces_colors = {}
        if surfaces_opacities is None:
            surfaces_opacities = {}
        if bounding_boxes is None:
            bounding_boxes = {}
        if visible_times is None:
            visible_times = {}

        for surface_id in self.plotters:
            self.surfaces_colors[surface_id] = surfaces_colors.get(surface_id, (0., 0., 1.))
            self.surfaces_opacities[surface_id] = surfaces_opacities.get(surface_id, 0.5)
            self.bounding_boxes[surface_id] = bounding_boxes.get(surface_id, BoundingBox())
            self.visible_times[surface_id] = visible_times.get(surface_id, lambda t: True)

    def set_data_file(self, surface_file_name):
        with open(surface_file_name) as data_file:
            data_lines = data_file.readlines()

        # Inspect the first line
        line = data_lines[0]
        words = line.split(", ")

        id_idx = [i for i in range(len(words)) if words[i].upper().startswith('ID')]
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
            if len(data.shape) == 1:
                self.data[data[i]] = data[i]
            else:
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
                print("No surface data at time ", t)
                sys.exit(1)

        plot_order = self.plot_order
        if self.plot_order is None:
            plot_order = self.plotters.keys()

        for surface_id in plot_order:
            if self.visible_times[surface_id](t):
                plotter = self.plotters[surface_id]
                plotter.bounding_box = self.bounding_boxes[surface_id]
                plotter_data = self.plotter_data[surface_id]
                data = data_line[plotter_data.start_idx:plotter_data.end_idx]
                plotter.plot(data, self.surfaces_colors[surface_id], self.surfaces_opacities[surface_id], t)


if __name__ == '__main__':
    simulation_directory = '../results/cyclic_triaxial/test/'
    surfaces_plotter = SurfacesPlotter(simulation_directory + 'surface_positions.dat')
    bbox = BoundingBox()
    bbox.z_min = -0.01
    bbox.z_max = 0.05
    surfaces_plotter.plot()

    mlab.show()
