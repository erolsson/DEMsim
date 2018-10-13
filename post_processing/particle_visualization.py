import glob
import os
import re

import numpy as np
from mayavi import mlab


class SpheresPlotter:
    def __init__(self):
        self.ms = None

    def plot(self, data, color=(184./255, 115./255., 51./255.)):
        x = data[:, 1]
        y = data[:, 2]
        z = data[:, 3]
        R = data[:, 7]
        if self.ms is None:
            self.ms = mlab.points3d(x, y, z, 2 * R,
                                    color=color,
                                    resolution=32,
                                    scale_factor=1., scale_mode='scalar').mlab_source
        else:
            self.ms.set(x=x, y=y, z=z)


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
        f = mlab.gcf()
        n = len(frame_times)
        for i, t in enumerate(frame_times):
            particle_data = np.genfromtxt(simulation_directory + 'particles_' + str(t) + '.dat', delimiter=',')
            spheres_plotter.plot(particle_data)
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


if __name__ == '__main__':
    simulation_directory = '../results/gyratory_compaction/1/'
    mlab.figure(size=(1920, 1200))
    animate_simulation(simulation_directory, delay=10, save_frames=True, save_directory=simulation_directory+'figures')

    mlab.show()
