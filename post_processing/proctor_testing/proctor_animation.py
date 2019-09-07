from mayavi import mlab

from post_processing.animation import Animation

if __name__ == '__main__':
    simulation_directory = '../results/stone_compaction/8-16mm/animation/'

    mlab.figure(size=(1920, 1200), bgcolor=(1., 1., 1.))
    animation = Animation(simulation_directory)
    animation.save_directory = '../results/stone_compaction/8-16mm/fig'
    animation.save_frames = True
    animation.delay = 0.01
    animation.start_time = 3.4
    """
    animation.surfaces_colors[0] = (0., 1., 0.)
    animation.surfaces_colors[5001] = (1., 0., 0.)
    animation.surfaces_opacities[0] = 0.25
    animation.surfaces_opacities[5001] = 0.8
    animation.plot_order = [5001, 0, 5002]
    animation.visible_functions[5002] = plate_visible
    """
    animation.run()
    mlab.show()
