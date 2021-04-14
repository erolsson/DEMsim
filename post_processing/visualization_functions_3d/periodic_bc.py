from collections import OrderedDict
import sys

import numpy as np
from mayavi import mlab


class PeriodicBC:
    def __init__(self, data_file):
        self.color = (0, 0, 0)
        self.line = '--'
        self.side_width = None
        self.counter = 0
        self.data = OrderedDict()
        self.set_data_file(data_file)
        self.box = []
        self.xmin = None
        self.xmax = None
        self.ymin = None
        self.ymax = None
        self.zmin = None
        self.zmax = None

    def set_data_file(self, bc_file_name):
        data_lines = np.genfromtxt(bc_file_name, delimiter=',')
        times = data_lines[:, 0]
        self.data = OrderedDict(zip(times, data_lines[:, 1:]))

    def plot(self, t=0.):
        if t is None:
            data = self.data[list(self.data.keys())[self.counter]]
            self.counter += 1
        else:
            try:
                data = self.data[t]
                self.counter = 0
            except LookupError:
                print("No surface data at time ", t)
                sys.exit(1)
        self.draw(data)

    def draw(self, data):
        d = data[1:7:2] - data[0:6:2]
        points = np.zeros((8, 3))

        for i in range(3):
            fixed_positions = [self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax]
            if d[i] == 0:
                if fixed_positions[2*i] is None and fixed_positions[2*i + 1] is None:
                    dmin = np.min(data[[0, 2, 4] != 2*i])
                    dmax = np.max(data[[1, 3, 5] != 2*i+1])
                    data[2*i] = self.side_width[0] if self.side_width is not None else dmin
                    data[2*i+1] = self.side_width[1] if self.side_width is not None else dmax
                else:
                    data[2*i] = fixed_positions[2*i]
                    data[2*i + 1] = fixed_positions[2*i + 1]

        points[0, 0], points[3, 0], points[4, 0], points[7, 0] = 4*[data[0]]
        points[1, 0], points[2, 0], points[5, 0], points[6, 0] = 4*[data[1]]

        points[0, 1], points[1, 1], points[4, 1], points[5, 1] = 4*[data[2]]
        points[2, 1], points[3, 1], points[6, 1], points[7, 1] = 4*[data[3]]
        points[0:4, 2] = data[4]
        points[4:, 2] = data[5]

        lines = set()
        x_lines = [(0, 1), (2, 3), (4, 5), (6, 7)]
        y_lines = [(0, 3), (1, 2), (4, 7), (5, 6)]
        z_lines = [(0, 4), (1, 5), (2, 6), (3, 7)]

        if d[0] != 0:
            lines.update(y_lines)
            lines.update(z_lines)
        if d[1] != 0:
            lines.update(x_lines)
            lines.update(z_lines)
        if d[2] != 0:
            lines.update(x_lines)
            lines.update(y_lines)

        if len(self.box) == 0:
            for line in lines:
                self.box.append(mlab.plot3d(points[list(line), 0], points[list(line), 1], points[list(line), 2],
                                            color=self.color, tube_radius=0.01*np.max(d)).mlab_source)
        else:
            for i, line in enumerate(lines):
                self.box[i].set(x=points[list(line), 0], y=points[list(line), 1], z=points[list(line), 2])


if __name__ == '__main__':
    periodic_bc = PeriodicBC('../../results/periodic_bc_test/sim_1/periodic_bc.dou')
    periodic_bc.plot(0.01)
    mlab.show()
