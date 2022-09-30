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
configurations = ["Small_Small", "Big_Small", "Big_Big"]
pressures = ["0kPa", "100kPa", "400kPa", "800kPa"]
simulations = [1, 2, 3]
pressure_vec = [float(p.replace("kPa", "")) for p in pressures]

for simulation, c in zip(configurations, ['g', 'r', 'b', 'm']):
    max_f = np.zeros((len(pressures), len(simulations)))
    for i, p in enumerate(pressures):
        for j, sim in enumerate(simulations):
            directory = (main_directory / str(sim) / (simulation + "_" + str(p) + "kPa")
                         / "shear_test")
            surface_forces = np.genfromtxt(directory / "surface_forces.dou", delimiter=",")
            max_f[i, j] = np.max(surface_forces[:, -4])

    plt.plot(pressure_vec, np.mean(max_f, axis=1), '-' + c, lw=3, label=str(simulation).replace("_", "-") + " kPa")
    plt.errorbar(pressure_vec, np.mean(max_f, axis=1), np.std(max_f, axis=1), fmt="none", elinewidth=2,
                 ecolor=c)

plt.xlim(-50, 900)
plt.xlabel("Pressure$ [kPa]", fontsize=24)
plt.ylabel("Maximum force [kN]", fontsize=24)
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("pressure_dependency.png")

plt.show()

