import pathlib
from mayavi import mlab

from visualization_functions_3d.animation import Animation


def main():
    simulation_directory = str(pathlib.Path("~/DEMsim/results/periodic_bc_test/sim_1/").expanduser())
    mlab.figure(size=(1920, 1200), bgcolor=(1., 1., 1.))
    animation = Animation(simulation_directory)
    animation.save_directory = 'animation/imgs/'
    animation.save_frames = True
    animation.delay = 0.0001
    animation.plot_periodic_bc = True
    # animation.mirror_particles = True
    animation.start_time = 0.
    animation.view_surfaces = True
    animation.run()
    mlab.show()


if __name__ == '__main__':
    main()
