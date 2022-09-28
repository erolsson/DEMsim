from collections import namedtuple
import os
import numpy as np

from mayavi import mlab
from tvtk.tools import visual

Particle = namedtuple('Particle', ['position', 'radius'])


def arrow(p1, p2, radius, color, cone_radius=0.1, fraction_head=0.2):
    p1 = np.array(p1)
    p2 = np.array(p2)
    ar1 = visual.arrow(x=p1[0], y=p1[1], z=p1[2])
    arrow_length = np.linalg.norm(p1 - p2)
    ar1.length_cone = 0.2
    ar1.radius_shaft = radius
    ar1.radius_cone = cone_radius
    ar1.actor.scale = [arrow_length, arrow_length, arrow_length]
    ar1.pos = ar1.pos/arrow_length
    ar1.axis = p2 - p1
    ar1.color = color

    return ar1


class AsphaltContactPlotter:
    def __init__(self, directory):
        self.ms = None
        self.directory = directory
        self.binder_radius = 0.001
        self.resolution = 20
        self.mid_plane = 0.
        self.lt = 0.01

    def plot(self, time):
        contact_data = np.genfromtxt(self.directory / ('contacts/contacts_' + str(time) + '.dou'), delimiter=',')
        particle_data = np.genfromtxt(self.directory / ('particles/particles_' + str(time) + '.dou'), delimiter=',')

        surfaces = set()
        if os.path.isfile(self.directory / 'surface_positions.dou'):
            with open(self.directory / 'surface_positions.dou', 'r') as surface_data_file:
                for line in surface_data_file.readlines():
                    data_line = line.split(',')
                    line_time = float(data_line[-1])
                    if time == line_time:
                        data_line = [word.strip() for word in data_line]
                        surface_ids = [int(word.split('=')[1]) for word in data_line if word.startswith('ID')]
                        surfaces.update(surface_ids)

        particles = {}
        for p in particle_data:
            particles[int(p[0])] = Particle(position=p[1:4], radius=float(p[7]))
        points_1 = []
        points_2 = []
        forces = []
        tangential_forces = []
        tan_points_1 = []
        tan_points_2 = []
        for contact in contact_data:
            force = float(contact[6])/1000
            ft = float(contact[10])/1000
            if force != 0:   # contact[6] is force
                obj1 = int(contact[0])
                obj2 = int(contact[1])
                point_1 = particles[obj1].position
                n = -np.array([float(contact[2]), float(contact[3]), float(contact[4])])
                length = particles[obj1].radius - float(contact[5])    # contact[h] = overlap
                mid_point = point_1+n*length
                if obj2 in particles:
                    length += particles[obj2].radius
                    point_2 = point_1 + n*length
                    if point_1[2] < self.mid_plane < point_2[2] or point_1[2] > self.mid_plane > point_2[2]:
                        points_1.append(point_1)
                        points_2.append(point_2)
                        forces.append(force)
                        ft_vec = contact[7:10]/1000
                        nt = ft_vec/ft
                        tangential_forces.append(ft)
                        if point_1[2] < point_2[2]:
                            tan_points_1.append(mid_point - 0.5*nt*self.lt)
                            tan_points_2.append(mid_point + 0.5*nt*self.lt)
                        else:
                            tan_points_1.append(mid_point + 0.5*nt*self.lt)
                            tan_points_2.append(mid_point - 0.5*nt*self.lt)

        for point_1, point_2, tan_point_1, tan_point_2,  force, ft, in zip(points_1, points_2, tan_points_1,
                                                                           tan_points_2, forces, tangential_forces):
            mlab.plot3d([point_1[0], point_2[0]], [point_1[1], point_2[1]], [point_1[2], point_2[2]],
                        [force, force], tube_radius=1e-3, vmax=np.max(forces), vmin=np.min(forces),
                        line_width=2)
            pts = mlab.points3d([point_2[0]], [point_2[1]], [point_2[2]], [force], scale_factor=1,
                                vmax=np.max(forces), vmin=np.min(forces))
            pts.glyph.scale_mode = 'scale_by_vector'
            pts.mlab_source.dataset.point_data.vectors = np.tile(2/np.sqrt(3)*1e-3, (3, 1)).transpose()
            pts = mlab.points3d([point_1[0]], [point_1[1]], [point_1[2]], [force], scale_factor=1,
                                vmax=np.max(forces), vmin=np.min(forces))
            pts.glyph.scale_mode = 'scale_by_vector'
            pts.mlab_source.dataset.point_data.vectors = np.tile(2/np.sqrt(3)*1e-3, (3, 1)).transpose()
            interval_length = (np.max(forces) - np.min(forces))/2
            # Small forces,

            if ft < (interval_length + min(forces)):
                c = (ft - np.min(forces))/interval_length
                color = (0, 1-c, c)
            else:
                c = (ft - (interval_length + min(forces)))/interval_length
                color = (c, 1-c, 0)
            arrow(tan_point_1, tan_point_2, radius=5e-2, color=color)

        color_bar = mlab.colorbar(title="Contact force [kN]")
        color_bar.data_range = [np.min(forces), np.max(forces)]
