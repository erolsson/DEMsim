import os

import matplotlib.pyplot as plt
import matplotlib

from pressure_density import relative_density_cylinder
from pressure_density import pressures_cylinder

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


simulation_directory = '../results/closed_die_compaction/4'
figure_save_directory = 'storakers_seminar/animation/graph_figures/'

d = relative_density_cylinder(simulation_directory)[:, 1]
p = pressures_cylinder(simulation_directory)[:, 1:]
total_time_points = d.shape[0]

p = p[d > d[0], :]/200e6
d = d[d > d[0]]*100
interesting_time_points = d.shape[0]
filling_tp = total_time_points - interesting_time_points

if not os.path.isdir(figure_save_directory):
    os.makedirs(figure_save_directory)

for tp in range(interesting_time_points):
    print "Generating figure", tp
    fig = plt.figure(0, (6, 6))
    plt.plot(d[:tp+1], p[:tp+1, 2], 'b', lw=3, label='Upper pressure')
    plt.plot(d[:tp+1], p[:tp+1, 1], 'r', lw=3, label='Lower pressure')
    plt.plot(d[:tp+1], p[:tp+1, 0], 'g', lw=3, label='Radial pressure')

    plt.plot(d[tp], p[tp, 2], 'bp', lw=3, ms=12)
    plt.plot(d[tp], p[tp, 1], 'rp', lw=3, ms=12)
    plt.plot(d[tp], p[tp, 0], 'gp', lw=3, ms=12)

    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([0.14, 0.12, box.width, box.height])

    plt.xlabel('Packing density [\%]')
    plt.ylabel(r'Pressure $p/\sigma_Y$')

    plt.xlim(d[0], 80)
    plt.ylim(0, 1.6)

    plt.grid(True)
    plt.legend(loc='upper left')
    # plt.tight_layout()
    idx_str = '0'*(len(str(total_time_points))-len(str(filling_tp + tp + 1))) + str(filling_tp + tp + 1)
    plt.savefig(figure_save_directory + 'graph' + idx_str + '.png')
    fig.clf()
