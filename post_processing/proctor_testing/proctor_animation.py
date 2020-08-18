from mayavi import mlab

from post_processing.animation import Animation

if __name__ == '__main__':
    simulation_directory = '../results/proctor/animation'

    mlab.figure(size=(1500, 800), bgcolor=(1., 1., 1.))
    animation = Animation(simulation_directory)
    animation.save_directory = '../post_processing/proctor/animation2/imgs/'
    animation.save_frames = True
    animation.delay = 0.01
    animation.start_time = 82.5
    animation.end_time = 83.7
    animation.surfaces_plotter.plotters[0].length_extension = 0.2
    animation.surfaces_colors[1] = (1., 0., 0.)
    animation.surfaces_opacities[1] = 1.

    def hammer_visible(t):
        return t > 83.27

    animation.visible_functions[1] = hammer_visible
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
