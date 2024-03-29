import pathlib
import numpy as np

from scipy.ndimage import uniform_filter1d

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
exp_directory = pathlib.Path("~/asphalt_bond_strength/experiments/bonded").expanduser()

fig = plt.figure(1)
fig.set_size_inches(11., 6., forward=True)
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([0.1, 0.15, 0.55, box.height])

fig = plt.figure(0)
fig.set_size_inches(8., 6., forward=True)
ax = plt.subplot(111)
ax.set_position([0.1*11/8, 0.15, 0.55*11/8, box.height])

area = np.pi*50*50
sizes = {"Small_Small": "(5.5/5.5)", "Big_Small": "(5.5/9.5)", "Big_Big": "(9.5/9.5)"}

for fig_number, (p, fig_label) in enumerate(zip(["100kPa", "400kPa"], ["a", "b"])):
    plt.figure(fig_number)
    ax = plt.subplot(111)
    for simulation, c in zip(sizes.keys(), ['g', 'r', 'b', 'm']):
        simulations = [1, 2, 3]
        data = np.genfromtxt(exp_directory / (simulation.lower() + "_" + p + ".dat"))
        plt.plot(data[:, 0], data[:, 1]*1000/area, c, lw=3, label=sizes[simulation])
        for sim in simulations:
            directory = main_directory / str(sim) / (simulation.lower() + "_" + p)/"shear_test"
            surface_forces = np.genfromtxt(directory / "surface_forces.dou", delimiter=",")
            surface_positions = np.genfromtxt(directory / "surface_positions.dou", delimiter=",")
            # plt.plot(surface_positions[:, -15]*1000, surface_forces[:, -4]/1000, c + '--')
            if sim == simulations[0]:
                d = surface_positions[:, -15]
                f = -surface_forces[:, -4]
                if d.shape != f.shape:
                    size = min(d.shape[0], f.shape[0])
                    d = d[0:size]
                    f = f[0:size]
            else:
                d += surface_positions[:, -15]
                f += -surface_forces[:, -4]

        plt.plot(d*1000/len(simulations), -uniform_filter1d(f, size=50)/len(simulations)/area, c + '--', lw=3)
        plt.text(0.5, 0.1, r"$\boldsymbol{\sigma_n = " + p.replace("kPa", r"\,\mathrm{kPa}") + "}$",
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform=ax.transAxes)

    plt.figure(fig_number)
    ax = plt.subplot(111)
    plt.ylim(0)
    plt.xlim(0, 7.5)

    if fig_number == 1:
        plt.plot([-2, -1], [0, 0], 'w', label='white')
        plt.plot([-2, -1], [0, 0], 'k', lw=3, label='Experiments')
        plt.plot([-2, -1], [0, 0], '--k', lw=3, label='Simulations')

        legend = ax.legend(loc='upper left', bbox_to_anchor=(1., 1.035), numpoints=1)
        legend.get_texts()[3].set_color("white")
    plt.xlabel("Displacement [mm]")
    plt.ylabel("Shear stress [MPa]")
    plt.text(0.08, 0.92, r'\bf{(' + fig_label + ')}', horizontalalignment='center', verticalalignment='center',
             transform=ax.transAxes)
    plt.savefig("bonded_" + p + ".png")
plt.show()
