import pathlib
import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

main_directory = pathlib.Path("~/DEMsim/results/asphalt_shear_box").expanduser()
for sim, c in zip(["big_big", "big_small"], ['r', 'b']):
    for p, line in zip(["100kPa", "400kPa"], ['--', '-']):
        directory = main_directory / (sim + "_" + p)
        surface_forces = np.genfromtxt(directory / "surface_forces.dou", delimiter=",")
        f = -surface_forces[:, -4]
        surface_positions = np.genfromtxt(directory / "surface_positions.dou", delimiter=",")
        d = surface_positions[:, -5]
        if d.shape[0] != f.shape[0]:
            n = min(d.shape[0], f.shape[0])
            f = f[0:n]
            d = d[0:n]
        if sim != "big_small" or p != "400kPa":
            plt.plot(d, f/1000, c + line, lw=2)

plt.ylim(0, 3)
plt.show()
