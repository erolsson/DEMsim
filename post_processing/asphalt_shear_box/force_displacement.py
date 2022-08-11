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

main_directory = pathlib.Path("~/DEMsim/results/asphalt_shear_box").expanduser()
for sim, c in zip(["Small_Small", "Big_Big", "Big_Small"], ['g', 'r', 'b']):
    for p, line in zip(["100kPa", "400kPa"], ['--', '-']):
        directory = main_directory / (sim.lower() + "_" + p)
        surface_forces = np.genfromtxt(directory / "surface_forces.dou", delimiter=",")
        f = -surface_forces[:, -4]
        surface_positions = np.genfromtxt(directory / "surface_positions.dou", delimiter=",")
        d = surface_positions[:, -5]
        if d.shape[0] != f.shape[0]:
            n = min(d.shape[0], f.shape[0])
            f = f[0:n]
            d = d[0:n]
        if line == '-':
            label = sim.replace('_', '-')
        else:
            label = None

        plt.plot(d*1000, f/1000, c + line, lw=2, label=label)

plt.ylim(0, 3)
plt.xlim(0, 8)

plt.plot([-1, -2], [-1, -2], 'w', label='white')
plt.plot([-1, -2], [-1, -2], '--k', lw=2, label='100 kPa')
plt.plot([-1, -2], [-1, -2], '-k', lw=2, label='400 kPa')

legend = ax.legend(loc='upper left', bbox_to_anchor=(1., 1.035), numpoints=1)
legend.get_texts()[3].set_color("white")
plt.gca().add_artist(legend)
plt.xlabel("Displacement [mm]")
plt.ylabel("Force [kN]")
plt.savefig("no_binder.png")
plt.show()
