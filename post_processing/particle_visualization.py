from mayavi import mlab

from animation import Animation


if __name__ == '__main__':
    simulation_directory = '../results/gyratory_compaction/2/'
    mlab.figure(size=(1920, 1200))
    animation = Animation(simulation_directory)
    animation.save_directory = 'figures'
    animation.save_frames = True
    animation.delay = 0.01
    animation.run()

    mlab.show()
