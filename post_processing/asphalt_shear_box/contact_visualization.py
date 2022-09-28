import pathlib

import numpy as np
import mayavi
from mayavi import mlab
from tvtk.tools import visual

from post_processing.visualization_functions_3d.plotting_functions import BoundingBox
from post_processing.visualization_functions_3d.snapshot import Snapshot
from asphalt_contact_plotter import AsphaltContactPlotter

fig = mlab.figure(size=(1024, 768), bgcolor=(1., 1., 1.), fgcolor=(0, 0., 0.))
visual.set_viewer(fig)

scene = fig.scene
scene.camera.position = [-0, -0.5, -0.5]
scene.camera.focal_point = [0., 0, 0]
# scene.camera.view_angle = 45.0
sim_name = "big_big_400kPa"
time = 23.477

directory = pathlib.Path("~/DEMsim/results/asphalt_shear_box/bonded").expanduser() / sim_name
snapshot = Snapshot(directory, contact_plotter_class=AsphaltContactPlotter)
surface_positions = np.genfromtxt(directory / "surface_positions.dou", delimiter=",")
mid_plane = surface_positions[-1, -3] + surface_positions[-1, -2]
snapshot.contact_plotter.mid_plane = mid_plane
snapshot.plot_periodic_bc = False
snapshot.particle_opacity = 1.
snapshot.particle_bounding_box.z_min = lambda t: mid_plane
bbox = BoundingBox()
bbox.z_max = lambda t: 0.03
bbox.z_min = lambda t: -0.03
snapshot.surface_bounding_boxes[2] = bbox
snapshot.surface_bounding_boxes[3] = bbox
snapshot.surfaces_opacities = {0: 0., 1: 0.5, 2: 0.5, 3: 0.0}
snapshot.plot(time)
scene.camera.view_up = [0, 0, -1]
mayavi.mlab.move(forward=0.05, right=None, up=-0.03)
mlab.savefig("contact_forces_" + sim_name + ".png", size=(1024, 1024))
mlab.show()
