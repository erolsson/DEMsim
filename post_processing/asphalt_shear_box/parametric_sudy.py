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
simulations = [1, 2, 3]

size_ratio = [c[0]/c[1] for c in configurations]

for p in pressures:
    max_f = np.zeros((len(configurations), len(simulations)))
    for i, config in enumerate(configurations):
        for j, sim in enumerate(simulations):
            directory = (main_directory / str(sim) / (str(config[0]) + "_" + str(config[1]) + "_" + str(p) + "kPa")
                         / "shear_test")
            surface_forces = np.genfromtxt(directory / "surface_forces.dou", delimiter=",")
            max_f[i, j] = np.max(surface_forces[:, -4])
    plt.plot(size_ratio, np.mean(max_f, axis=1), '-x', lw=2, ms=12)

plt.show()
