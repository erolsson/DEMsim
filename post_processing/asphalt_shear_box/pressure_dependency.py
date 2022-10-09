import pathlib
import numpy as np

from scipy.optimize import fmin

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 24})
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath} \usepackage{gensymb}"
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

main_directory = pathlib.Path("~/DEMsim/results/asphalt_shear_box/bonded_plane_wall_friction_2/").expanduser()
exp_directory = pathlib.Path("~/asphalt_bond_strength/experiments").expanduser()

configurations = ["Small_Small", "Big_Small", "Big_Big"]
sizes = {"Small_Small": "(5.5/5.5)", "Big_Small": "(5.5/9.5)", "Big_Big": "(9.5/9.5)"}

pressures = [0, 100, 400, 800]
simulations = [1, 2, 3]

fig = plt.figure(0)
fig.set_size_inches(13., 6., forward=True)
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([0.1, 0.15, 0.5, box.height])
plt.plot([-1, -2], [-1, -1], 'w', label=r"\bf{DEM}")


def residual(par, data):
    r = 0
    for i, (p, f) in enumerate(data):
        r += np.sum(((par[0]*p + par[i+1] - f)**2))
    return r


def main():
    area = np.pi*50*50
    data = []

    for k, (simulation, c) in enumerate(zip(configurations, ['g', 'r', 'b', 'm'])):
        max_f = np.zeros((len(pressures), len(simulations)))

        for i, p in enumerate(pressures):
            for j, sim in enumerate(simulations):
                directory = (main_directory / str(sim) / (simulation.lower() + "_" + str(p) + "kPa")
                             / "shear_test")
                surface_forces = np.genfromtxt(directory / "surface_forces.dou", delimiter=",")
                max_f[i, j] = np.max(surface_forces[:, -4])

        f = np.mean(max_f, axis=1)
        plt.plot(pressures, f/area, '-x' + c, lw=3, label=sizes[simulation], ms=16, mew=2)
        data.append((np.array(pressures[1:]), np.array(f[1:]/area)))

    par = fmin(residual, [0.001, 0.7, 0.8, 0.5], args=(data, ), maxfun=1e6, maxiter=1e6)
    x = np.linspace(0, 800, 1000)
    for i, (simulation, c) in enumerate(zip(configurations, ['g', 'r', 'b', 'm'])):
        plt.plot(x, par[0]*x + par[i+1], '--' + c, lw=2)

    print(par)
    plt.plot([-1, -2], [-1, -1], 'w', label="white")
    plt.plot([-1, -2], [-1, -1], '--k', label=r"Eq. (11)")
    plt.plot([-1, -2], [-1, -1], 'ks', ms=8, label=r"Exp. data \n Raab et al (2011)")
    for k, (simulation, c) in enumerate(zip(configurations, ['g', 'r', 'b', 'm'])):
        exp_data = np.genfromtxt(exp_directory / ("shear_stress_0kPa_" + simulation.lower() + ".csv"), delimiter=",")
        exp_p = np.round(exp_data[:, 0], 1)*1e3
        plt.plot(exp_p, exp_data[:, 1], c + "s", ms=8)

    plt.xlim(-50, 900)
    plt.ylim(0, 2.2)
    plt.text(800, 0.75, r"$\boldsymbol{(5.5/5.5)\, c_0 = 1.10\, \mathrm{MPa}}$", horizontalalignment='right', color='g')
    plt.text(800, 0.55, r"$\boldsymbol{(5.5/9.5)\, c_0 = 1.27\, \mathrm{MPa}}$", horizontalalignment='right', color='r')
    plt.text(800, 0.35, r"$\boldsymbol{(9.5/9.5)\, c_0 = 0.85\, \mathrm{MPa}}$", horizontalalignment='right', color='b')
    plt.text(680, 0.15, r"$\boldsymbol{\varphi = 41.5 \degree}$", horizontalalignment='right')
    plt.xlabel(r"Confining stress $\sigma_n$ [kPa]", fontsize=24)
    plt.ylabel("Maximum shear stress [MPa]", fontsize=24)
    legend = ax.legend(loc='upper left', bbox_to_anchor=(1., 1.035), numpoints=1)
    legend.get_texts()[4].set_color("white")
    plt.gca().add_artist(legend)

    plt.savefig("pressure_dependency.png")

    plt.show()


if __name__ == '__main__':
    main()
