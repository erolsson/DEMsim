from collections import namedtuple
import os
import numpy as np

from mayavi import mlab

from visualization_functions_3d import colors

Particle = namedtuple('Particle', ['position', 'radius'])


class BatteryContactPlotter:
    def __init__(self, directory, bounding_box=None):
        self.ms = None
        self.bounding_box = bounding_box
        self.directory = directory
        self.binder_radius = 0.01
        self.resolution = 20
        self.color = colors.blue

    def plot(self, time):
        contact_data = np.genfromtxt(self.directory + '/contacts/contacts_' + str(time) + '.dou', delimiter=',')
        particle_data = np.genfromtxt(self.directory + '/particles/particles_' + str(time) + '.dou', delimiter=',')
        periodic_bc_data = np.genfromtxt(self.directory + '/periodic_bc.dou', delimiter=',')
        surfaces = set()
        with open(self.directory + '/surface_positions.dou', 'r') as surface_data_file:
            for line in surface_data_file.readlines():
                data_line = line.split(',')
                line_time = float(data_line[-1])
                if time == line_time:
                    data_line = [word.strip() for word in data_line]
                    surface_ids = [int(word.split('=')[1]) for word in data_line if word.startswith('ID')]
                    surfaces.update(surface_ids)
        periodic_bc = periodic_bc_data[periodic_bc_data[:, 0] == time, 1:]
        print(periodic_bc)
        particles = {}
        for p in particle_data:
            particles[int(p[0])] = Particle(position=p[1:4], radius=p[7])
        for contact in contact_data:
            if int(contact[6]) != 0:   # contact[6] is force
                obj1 = int(contact[0])
                obj2 = int(contact[1])
                point_1 = particles[obj1].position
                n = -np.array([float(contact[2]), float(contact[3]), float(contact[4])])
                length = particles[obj1].radius - float(contact[5])    # contact[h] = overlap

                q, z0 = np.meshgrid(np.linspace(0, 2*np.pi, self.resolution),
                                   np.linspace(0, length, self.resolution))
                x0 = self.binder_radius*np.cos(q)
                y0 = self.binder_radius*np.sin(q)

                if n[0] != 0. and n[1] != 0:
                    ex = np.array([n[1], -n[0], 0])/np.sqrt(n[0]**2 + n[1]**2)
                else:
                    ex = np.array([1, 0, 0])
                ey = np.cross(n, ex)
                print(ex, ey, n)
                x = x0*ex[0] + y0*ey[0] + z0*n[0] + point_1[0]
                y = x0*ex[1] + y0*ey[1] + z0*n[1] + point_1[1]
                z = x0*ex[2] + y0*ey[2] + z0*n[2] + point_1[2]

                mlab.mesh(x, y, z, color=self.color, transparent=False, reset_zoom=False)


def main():
    contact_plotter = BatteryContactPlotter(os.path.expanduser('~/DEMsim/results/battery_rve/electrode_elaheh/'))
    contact_plotter.plot(16.4507)
    mlab.show()


if __name__ == '__main__':
    main()
