from mayavi import mlab

from animation_functions import animate_simulation

if __name__ == '__main__':
    simulation_directory = '../results/gyratory_compaction/1/'
    mlab.figure(size=(1920, 1200))
    animate_simulation(simulation_directory, delay=10, save_frames=True, save_directory=simulation_directory+'figures')

    mlab.show()
