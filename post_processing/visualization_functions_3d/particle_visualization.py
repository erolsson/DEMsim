from mayavi import mlab

from animation import Animation
from colors import red

if __name__ == '__main__':
    simulation_directory = '../../results/asphalt_shear_box/big_small'

    # Doing some inspection to construct a good bounding_box figure
    # dimension_data = dimensions_cylinder(simulation_directory)
    # time = dimension_data[:, 0]
    # height = dimension_data[:, 3]

    # def plate_visible(t):
    #     return t >= time[height < height[0]][0]

    mlab.figure(size=(1920, 1200), bgcolor=(1., 1., 1.))
    animation = Animation(simulation_directory)
    animation.save_directory = '../../results/viscoelastic/cubic_box_animation'
    animation.save_frames = True
    animation.delay = 0.01
    animation.start_time = 0.43
    animation.surfaces_colors[2] = red
    animation.surfaces_colors[3] = red
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
