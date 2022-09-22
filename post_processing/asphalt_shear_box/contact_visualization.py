import pathlib

from mayavi import mlab

from post_processing.visualization_functions_3d.snapshot import Snapshot

mlab.figure(size=(1024, 768), bgcolor=(1., 1., 1.), fgcolor=(0, 0., 0.))
snapshot = Snapshot(pathlib.Path("~/DEMsim/results/asphalt_shear_box/bonded/big_big_100kPa").expanduser())
snapshot.plot_periodic_bc = False
snapshot.surfaces_plotter = None
snapshot.particle_opacity = 0.1
snapshot.plot(12.26)
mlab.show()

