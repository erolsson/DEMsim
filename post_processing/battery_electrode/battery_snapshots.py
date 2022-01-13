import pathlib

import numpy as np
from mayavi import mlab

import post_processing.visualization_functions_3d.colors as colors
from post_processing.visualization_functions_3d.snapshot import Snapshot
from post_processing.visualization_functions_3d.battery_contact_plotter import BatteryContactPlotter


time = 1.71
filename = "~/batteries/particles_2021/DEM_fill.png"
battery_contact_plotter = BatteryContactPlotter

mlab.figure(size=(2560, 1440), bgcolor=(1., 1., 1.), fgcolor=(0, 0., 0.))
directory = pathlib.Path('~/DEMsim/results/elaheh/mechanical-test-experimental-5-cykles-E34bt01Rbr05/').expanduser()
surface_positions = np.genfromtxt(directory / "surface_positions.dou", delimiter=",")
periodic_cell = np.genfromtxt(directory / "periodic_bc.dou", delimiter=",")
p_cell = periodic_cell[periodic_cell[:, 0] == time, 1:5][0]
particles = np.genfromtxt(directory / "particles" / "particles_10.515.dou", delimiter=",")
# z_max = np.max(particles[:, 3] + particles[:, 7])
z_max = surface_positions[surface_positions[:, -1] == time, -2]
if isinstance(z_max, np.ndarray):
    z_max = z_max[0]

snapshot = Snapshot(directory, contact_plotter_class=battery_contact_plotter)
snapshot.create_periodic_bc_plotter()
snapshot.plot_periodic_bc = True
snapshot.periodic_bc_plotter.zmin = 0
snapshot.periodic_bc_plotter.zmax = z_max
snapshot.mirror_particles = False
if battery_contact_plotter:
    snapshot.contact_plotter.color = colors.grey
    snapshot.contact_plotter.binder_radius = 0.5*0.03
# snapshot.surfaces_plotter.plotters.pop(1)

for z in [0, z_max]:
    for x in p_cell[0:2]:
        for y in p_cell[2:4]:
            mlab.plot3d([x, 1.3*x], [y, y], [z, z], color=colors.black, tube_radius=0.01*p_cell[1])
            mlab.plot3d([x, x], [y, 1.3*y], [z, z], color=colors.black, tube_radius=0.01*p_cell[1])

snapshot.plot(time)

mlab.savefig(filename=str(pathlib.Path(filename).expanduser()), size=(320, 320))

mlab.show()
