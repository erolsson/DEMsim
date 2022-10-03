import pathlib
import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 24})
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

main_directory = pathlib.Path("~/DEMsim/results/asphalt_shear_box/bonded_plane_wall_friction_2/").expanduser()

configurations = [(100, 100), (100, 75), (100, 50), (1000, 375), (100, 25)]
pressures = [100, 400]
simulations = [1, 2, 3]

size_ratio = [c[0]/c[1] for c in configurations]

for p, color in zip(pressures, 'br'):
    max_f = np.zeros((len(configurations), len(simulations)))
    for i, config in enumerate(configurations):
        for j, sim in enumerate(simulations):
            directory = (main_directory / str(sim) / (str(config[0]) + "_" + str(config[1]) + "_" + str(p) + "kPa")
                         / "shear_test")
            surface_forces = np.genfromtxt(directory / "surface_forces.dou", delimiter=",")
            max_f[i, j] = np.max(surface_forces[:, -4])
    plt.plot(size_ratio, np.mean(max_f, axis=1), '-' + color, lw=3, label=str(p) + " kPa")
    plt.errorbar(size_ratio, np.mean(max_f, axis=1), np.std(max_f, axis=1), fmt="none", elinewidth=2,
                 ecolor=color)

    for simulation, symbol in zip(["Small_Small", "Big_Small", "Big_Big"], ['o', 'x', 's']):
        max_f = np.zeros(3)
        for j, sim in enumerate(simulations):
            directory = main_directory / str(sim) / (simulation.lower() + "_" + str(p) + "kPa") / "shear_test"
            surface_forces = np.genfromtxt(directory / "surface_forces.dou", delimiter=",")
            max_f[j] = np.max(surface_forces[:, -4])
        ratio = 1 if simulation != "Big_Small" else 9.5/5.5
        plt.plot(ratio, np.mean(max_f), symbol + color, lw=3, mew=3, ms=12)
        plt.errorbar(ratio, np.mean(max_f), np.std(max_f), fmt="none", elinewidth=2,
                     ecolor=color)

plt.xlim(0.5, 4.5)
plt.xlabel("Size ratio $D_1/D_2$ [-]", fontsize=24)
plt.ylabel("Maximum force [kN]", fontsize=24)
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("parametric_study.png")

plt.show()
