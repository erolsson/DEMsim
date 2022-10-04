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
configurations = ["Small_Small", "Big_Small", "Big_Big"]
pressures = [0, 100, 400, 800]
simulations = [1, 2, 3]


def residual(par, data):
    r = 0
    for i, (p, f) in enumerate(data):
        r += np.sum(((par[0]*p + par[i+1] - f)**2))
    return r


def main():
    data = []
    for simulation, c in zip(configurations, ['g', 'r', 'b', 'm']):
        max_f = np.zeros((len(pressures), len(simulations)))
        for i, p in enumerate(pressures):
            for j, sim in enumerate(simulations):
                directory = (main_directory / str(sim) / (simulation.lower() + "_" + str(p) + "kPa")
                             / "shear_test")
                surface_forces = np.genfromtxt(directory / "surface_forces.dou", delimiter=",")
                max_f[i, j] = np.max(surface_forces[:, -4])

        f = np.mean(max_f, axis=1)/1e3
        plt.plot(pressures, f, '-' + c, lw=3, label=str(simulation).replace("_", "-"))
        plt.errorbar(pressures, f, np.std(max_f, axis=1)/1e3, fmt="none", elinewidth=2,
                     ecolor=c)
        data.append((pressures[1:], f[1:]))
        # a, b = np.polyfit(pressures[1:], f[1:], 1)
        # x = np.linspace(0, 800, 1000)
        # plt.plot(x, a*x + b, '--' + c, lw=2)

    print(fmin(residual, [0, 0, 0, 0], data=(data, )))

    plt.xlim(-50, 900)
    plt.xlabel("Pressure [kPa]", fontsize=24)
    plt.ylabel("Maximum force [kN]", fontsize=24)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig("pressure_dependency.png")

    plt.show()


if __name__ == '__main__':
    main()
