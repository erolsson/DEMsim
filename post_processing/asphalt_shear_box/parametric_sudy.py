import pathlib
import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 24})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

main_directory = pathlib.Path("~/DEMsim/results/asphalt_shear_box/bonded_plane_wall_friction_2/").expanduser()

configurations = [(50, 100), (75, 100), (100, 100), (100, 75), (100, 50)]
pressures = [100, 400]
simulations = [1]

size_ratio = [c[0]/c[1] for c in configurations]

for p in pressures:
    max_f = np.zeros(len(configurations))
    for i, config in enumerate(configurations):
        for sim in simulations:
            directory = (main_directory / str(sim) / (str(config[0]) + "_" + str(config[1]) + "_" + str(p) + "kPa")
                         / "shear_test")
            surface_forces = np.genfromtxt(directory / "surface_forces.dou", delimiter=",")
            max_f[i] = np.max(surface_forces[:, -4])
    plt.plot(size_ratio, max_f, '-x', lw=2, ms=12)

plt.show()
