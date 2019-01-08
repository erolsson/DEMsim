import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

from post_processing.pressure_density import relative_density_cylinder

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


density = relative_density_cylinder('../../results/closed_die_compaction/4/')[:, 1]
total_timesteps = density.shape[0]
figure_save_directory = 'animation/animation_frames/'

if not os.path.isdir(figure_save_directory):
    os.makedirs(figure_save_directory)


def heading(increment):
    if density[increment] == density[0]:
        return "Filling"
    else:
        turning_idx = np.argmax(density)
        if increment <= turning_idx:
            return "Compaction"
        else:
            return "Unloading"


def create_frame(increment):
    fig = plt.figure(0, figsize=(11, 8))
    inc_str = '0'*(len(str(total_timesteps))-len(str(increment+1))) + str(increment+1)
    print "Generating figure", inc_str
    particle_image = plt.imread('animation/particle_figures/frame' + inc_str + '.png')
    image = np.ones((particle_image.shape[0], 1100, 4))
    image[:, :particle_image.shape[1], 0:3] = particle_image

    if density[increment] != density[0]:
        graph_image = plt.imread('animation/graph_figures/graph' + inc_str + '.png')
        image[particle_image.shape[0] - graph_image.shape[0]:, particle_image.shape[1]:, :] = graph_image
    plt.imshow(image)
    plt.text(800, 100, r'\bf{' + heading(increment) + '}', Fontsize=30,
             horizontalalignment='center', verticalalignment='center')
    plt.xticks([])
    plt.yticks([])
    plt.savefig(figure_save_directory + '/frame' + inc_str + '.png', dpi=300)
    fig.clf()


if __name__ == '__main__':
    for i in range(total_timesteps):
        create_frame(i)
