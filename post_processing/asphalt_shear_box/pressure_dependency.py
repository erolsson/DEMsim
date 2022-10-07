import pathlib
import numpy as np

from scipy.optimize import fmin

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
exp_directory = pathlib.Path("~/asphalt_bond_strength/experiments").expanduser()

configurations = ["Small_Small", "Big_Small", "Big_Big"]
pressures = [0, 100, 400, 800]
simulations = [1, 2, 3]


def residual(par, data):
    r = 0
    for i, (p, f) in enumerate(data):
        r += np.sum(((par[0]*p + par[i+1] - f)**2))
    return r


def main():
    area = np.pi*50*50
    data = []
    experiments = np.genfromtxt(exp_directory / "pressure_dependency.csv")
    for k, (simulation, c) in enumerate(zip(configurations, ['g', 'r', 'b', 'm'])):
        max_f = np.zeros((len(pressures), len(simulations)))

        for i, p in enumerate(pressures):
            for j, sim in enumerate(simulations):
                directory = (main_directory / str(sim) / (simulation.lower() + "_" + str(p) + "kPa")
                             / "shear_test")
                surface_forces = np.genfromtxt(directory / "surface_forces.dou", delimiter=",")
                max_f[i, j] = np.max(surface_forces[:, -4])

        f = np.mean(max_f, axis=1)
        plt.plot(pressures, f/area, '-o' + c, lw=3, label=str(simulation).replace("_", "-"), ms=10)
        data.append((np.array(pressures[1:]), np.array(f[1:])))

        mean = experiments[0, (1 + 2*k)::6]
        std = experiments[1, (1 + 2*k)::6]
        plt.plot([100, 400], mean, 'x' + c, lw=3, mew=3, ms=16 - 2*k)
        plt.errorbar([100, 400], mean, std, fmt="none", elinewidth=3,
                     ecolor=c)
    par = fmin(residual, [0, 0, 0, 0], args=(data, ), maxfun=1e6, maxiter=1e6)
    x = np.linspace(0, 800, 1000)
    for i, (simulation, c) in enumerate(zip(configurations, ['g', 'r', 'b', 'm'])):
        plt.plot(x, par[0]*x + par[i+1], '--' + c, lw=2)
    plt.xlim(-50, 900)
    plt.ylim(0, 2.2)
    plt.xlabel(r"Confining stress $\sigma_n$ [kPa]", fontsize=24)
    plt.ylabel("Maximum shear stress [MPa]", fontsize=24)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig("pressure_dependency.png")

    plt.show()


if __name__ == '__main__':
    main()
