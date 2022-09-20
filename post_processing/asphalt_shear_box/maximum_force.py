from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


# Bonded results
experimental_data_bonded = {"Big_Big": {100.: [], 400.: []},
                            "Small_Small": {100.: [], 400.: []},
                            "Big_Small": {100.: [], 400.: []}}
numerical_data_bonded = deepcopy(experimental_data_bonded)

with open("experiments_bonded.dat", 'r') as experiment_file:
    for line in experiment_file.readlines():
        config, pressure, force, std = line.split(",")
        experimental_data_bonded[config][float(pressure)] = [float(force), float(std)]

counter = 0
colors = {"Big_Big": 'b', "Small_Small": "g", "Big_Small": "r"}
f0 = experimental_data_bonded["Big_Big"][100][0]
for p in [100, 400]:
    for config in ["Big_Big", "Small_Small", "Big_Small"]:
        data = experimental_data_bonded[config][p]
        bar = plt.bar(2*counter, data[0]/f0, 1, yerr=data[1]/f0, color=colors[config], linewidth=2, capsize=10)
        bar.errorbar.lines[2][0].set_lw(3)
        bar.errorbar.lines[1][0].set_mew(3)
        bar.errorbar.lines[1][1].set_mew(3)

        counter += 1

plt.show()
