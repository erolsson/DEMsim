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

fig = plt.figure(0)
fig.set_size_inches(11., 6., forward=True)
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([0.1, 0.15, 0.55, box.height])

main_directory = pathlib.Path("~/DEMsim/results/asphalt_shear_box/bonded_plane/sim").expanduser()
for simulation, c in zip(["Small_Small", "Big_Small", "Big_Big"], ['g', 'r', 'b', 'm']):
    for p, line in zip(["100kPa", "400kPa"], ['--', '-']):
        n = 1
        for sim in range(1, n+1):
            directory = main_directory / str(sim) / (simulation.lower() + "_" + p)/"shear_test"
            surface_forces = np.genfromtxt(directory / "surface_forces.dou", delimiter=",")
            surface_positions = np.genfromtxt(directory / "surface_positions.dou", delimiter=",")
            kinetic_energy = np.genfromtxt(directory / "kinetic_energy.dou", delimiter=",")
            if sim == 1:
                d = surface_positions[:, -15]
                f = -surface_forces[:, -4]
            else:
                d += surface_positions[:, -15]
                f += -surface_forces[:, -4]
        if line == '-':
            label = simulation.replace('_', '-')
        else:
            label = None

        plt.figure(0)
        plt.plot(d*1000/n, -f/1000/n, c + line, lw=2, label=label)
plt.figure(0)
# plt.ylim(0, 3)
plt.xlim(0, 8)

plt.plot([-1, -2], [-1, -2], 'w', label='white')
plt.plot([-1, -2], [-1, -2], '--k', lw=2, label='100 kPa')
plt.plot([-1, -2], [-1, -2], '-k', lw=2, label='400 kPa')
plt.ylim(0)
legend = ax.legend(loc='upper left', bbox_to_anchor=(1., 1.035), numpoints=1)
legend.get_texts()[3].set_color("white")
plt.gca().add_artist(legend)
plt.xlabel("Displacement [mm]")
plt.ylabel("Force [kN]")
plt.savefig("no_bond.png")
plt.show()
