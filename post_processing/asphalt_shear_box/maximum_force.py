from copy import deepcopy
import pathlib

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import patches

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


fig = plt.figure(0)
fig.set_size_inches(11., 6., forward=True)
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([0.1, 0.15, 0.55, box.height])

# Bonded results
experimental_data_bonded = {"Big_Big": {100.: [], 400.: []},
                            "Small_Small": {100.: [], 400.: []},
                            "Big_Small": {100.: [], 400.: []}}
numerical_data_bonded = deepcopy(experimental_data_bonded)
simulations = [1, 2, 3]
simulation_directory = pathlib.Path("~/DEMsim/results/asphalt_shear_box/bonded_plane_wall_friction_2/").expanduser()

with open("experiments_bonded.dat", 'r') as experiment_file:
    for line in experiment_file.readlines():
        config, pressure, force, std = line.split(",")
        experimental_data_bonded[config][float(pressure)] = [float(force), float(std)]
        forces = []

        for sim in simulations:
            directory = simulation_directory / str(sim) / (config.lower() + "_" + str(int(pressure)) + "kPa") / "shear_test"
            surface_forces = np.genfromtxt(directory / "surface_forces.dou", delimiter=",")
            forces.append(np.max(surface_forces[:, -4])/1000)
        forces = np.array(forces)
        numerical_data_bonded[config][float(pressure)] = [np.mean(forces), np.std(forces, ddof=1)]

counter = 0
colors = {"Big_Big": 'b', "Small_Small": "g", "Big_Small": "r"}
f0_exp = experimental_data_bonded["Big_Big"][100][0]
f0_sim = numerical_data_bonded["Big_Big"][100][0]
l, = plt.plot([-1, -1], [-1, -1], 'w')
legend_entries = [None]*9
labels = [None]*9
legend_entries[0] = l
labels[0] = r"\textbf{Experiments}"
legend_entries[4] = l
labels[4] = "White"
legend_entries[5] = l
labels[5] = r"\textbf{Simulations}"

for p in [100, 400]:
    for config in ["Big_Big", "Small_Small", "Big_Small"]:
        data = experimental_data_bonded[config][p]
        label = config.replace("_", "-")
        bar = plt.bar(2*counter, data[0]/f0_exp, 1, yerr=data[1]/f0_exp, color=colors[config], linewidth=2, capsize=10,
                      align="edge", label=label)
        if p == 100:
            legend_entries[counter + 1] = bar
            labels[counter + 1] = label
        bar.errorbar.lines[2][0].set_lw(3)
        bar.errorbar.lines[1][0].set_mew(3)
        bar.errorbar.lines[1][1].set_mew(3)
        data = numerical_data_bonded[config][p]
        bar = plt.bar(2*counter + 1, data[0]/f0_sim, 1, yerr=data[1]/f0_sim, color=colors[config], linewidth=2,
                      capsize=10, alpha=0.5, align="edge")
        if p == 100:
            legend_entries[counter + 6] = bar
            labels[counter + 6] = label
        bar.errorbar.lines[2][0].set_lw(3)
        bar.errorbar.lines[1][0].set_mew(3)
        bar.errorbar.lines[1][1].set_mew(3)

        counter += 1

plt.plot([6, 6], [0, 2.5], ':k', lw=2)
plt.xlim(0, 12)
plt.ylim(0, 2.2)
plt.ylabel("Normalized maximum force", fontsize=24)
plt.text(1, 1.8, r'$\boldsymbol{p = 100}$ \textbf{kPa}', fontsize=24)
plt.text(7, 1.8, r'$\boldsymbol{p = 400}$ \textbf{kPa}', fontsize=24)
legend = ax.legend(legend_entries, labels, loc='upper left', bbox_to_anchor=(1., 1.035), numpoints=1)
legend.get_texts()[4].set_color("white")
plt.gca().add_artist(legend)
plt.xticks([], ())
plt.savefig("maximum_force")
plt.show()
