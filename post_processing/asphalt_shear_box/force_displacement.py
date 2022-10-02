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


for fig_number in [0, 1]:
    fig = plt.figure(fig_number)
    fig.set_size_inches(11., 6., forward=True)
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([0.1, 0.15, 0.55, box.height])

main_directory = pathlib.Path("~/DEMsim/results/asphalt_shear_box/bonded_plane_wall_friction_2/").expanduser()
exp_directory = pathlib.Path("~/asphalt_bond_strength/experiments/bonded").expanduser()
for fig, p in enumerate(["100kPa", "400kPa"]):
    plt.figure(fig)
    for simulation, c in zip(["Small_Small", "Big_Small", "Big_Big"], ['g', 'r', 'b', 'm']):
        simulations = [1, 2, 3]
        data = np.genfromtxt(exp_directory / (simulation.lower() + "_" + p + ".dat"))
        plt.plot(data[:, 0], data[:, 1], c, lw=3, label=simulation.replace('_', '-'))
        for sim in simulations:
            directory = main_directory / str(sim) / (simulation.lower() + "_" + p)/"shear_test"
            surface_forces = np.genfromtxt(directory / "surface_forces.dou", delimiter=",")
            surface_positions = np.genfromtxt(directory / "surface_positions.dou", delimiter=",")
            kinetic_energy = np.genfromtxt(directory / "kinetic_energy.dou", delimiter=",")
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

        plt.plot(d[::10]*1000/len(simulations), -f[::10]/1000/len(simulations), c + '--', lw=3)

for fig in range(2):
    plt.figure(fig)
    plt.ylim(0)
    plt.xlim(0, 8)
    plt.plot([-2, -1], [0, 0], 'w', label='white')
    plt.plot([-2, -1], [0, 0], 'k', lw=3, label='Experiments')
    plt.plot([-2, -1], [0, 0], '--k', lw=3, label='Simulations')

    legend = ax.legend(loc='upper left', bbox_to_anchor=(1., 1.035), numpoints=1)
    legend.get_texts()[3].set_color("white")
    plt.xlabel("Displacement [mm]")
    plt.ylabel("Force [kN]")
plt.show()
