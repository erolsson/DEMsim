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
            

if __name__ == '__main__':
    simulation_directory = '../results/gyratory_compaction/1/'
    time = 0.4
    
    mlab.figure(size=(1920, 1200))
    plotter = SpheresPlotter()
    particle_file_name = simulation_directory + 'particles_' + str(time) + '.dat'
    
    particle_data = np.genfromtxt(particle_file_name, delimiter=',')
    plotter.plot(particle_data)
    mlab.show()
