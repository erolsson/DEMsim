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

fig = plt.figure(0)
fig.set_size_inches(13., 6., forward=True)
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([0.1, 0.15, 0.5, box.height])

configurations = [(100, 130), (100, 100), (100, 75), (100, 50), (1000, 375), (100, 25)]
pressures = [400]
simulations = [1, 2, 3]
plt.plot([-1, -2], [1, 1], 'w', label=r"\bf{Table 2}")
size_ratio = [c[1]/c[0] for c in configurations]
area = np.pi*50*50
for p, color in zip(pressures, 'br'):
    max_f = np.zeros((len(configurations), len(simulations)))
    for i, config in enumerate(configurations):
        for j, sim in enumerate(simulations):
            directory = (main_directory / str(sim) / (str(config[0]) + "_" + str(config[1]) + "_" + str(p) + "kPa")
                         / "shear_test")
            surface_forces = np.genfromtxt(directory / "surface_forces.dou", delimiter=",")
            max_f[i, j] = np.max(surface_forces[:, -4])
    plt.plot(size_ratio, np.mean(max_f, axis=1)/area, '-' + color, lw=3, label="$D_2/D_1=x/10$ mm")
    plt.errorbar(size_ratio, np.mean(max_f, axis=1)/area, np.std(max_f, axis=1)/area, fmt="none", elinewidth=2,
                 ecolor=color)

    for simulation, symbol in zip(["Small_Small", "Big_Small", "Big_Big"], ['o', 'x', 's']):
        max_f = np.zeros(3)
        for j, sim in enumerate(simulations):
            directory = main_directory / str(sim) / (simulation.lower() + "_" + str(p) + "kPa") / "shear_test"
            surface_forces = np.genfromtxt(directory / "surface_forces.dou", delimiter=",")
            max_f[j] = np.max(surface_forces[:, -4])
        ratio = 1 if simulation != "Big_Small" else 5.5/9.5
        plt.plot(ratio, np.mean(max_f)/area, symbol + color, lw=3, mew=3, ms=12)
        plt.errorbar(ratio, np.mean(max_f)/area, np.std(max_f)/area, fmt="none", elinewidth=2,
                     ecolor=color)

plt.plot([-1, -2], [-1, -2], '-b', lw=3, label="$D_2/D_1=x/10$ mm")
plt.xlim(0.2, 1.5)
plt.ylim(0)

plt.plot([-1, -2], [1, 1], 'w', label=r"white")
plt.plot([-1, -2], [1, 1], 'w', label=r"\bf{Table 3}")
for simulation, symbol in zip(["(5.5/5.5)", "(5.5/9.5)", "(9.5/9.5)"], ['o', 'x', 's']):
    plt.plot([-2], [1], 'b' + symbol, ms=12, mew=3, label=simulation.replace('_', '-'))
plt.xlabel("Size ratio $D_2/D_1$ [-]", fontsize=24)
plt.ylabel("Maximum shear stress [MPa]", fontsize=24)
legend = ax.legend(loc='upper left', bbox_to_anchor=(1., 1.035), numpoints=1)
legend.get_texts()[2].set_color("white")
plt.savefig("parametric_study.png")

plt.show()
