import glob
import pickle
import re
import sys

import numpy as np


def get_fractured_particles_at_stroke(directory, layer_idx, stroke_idx):
    filename = directory + '/' + 'layer_' + str(layer_idx) + '/stroke_' + str(stroke_idx) + '/particle_cracks.dou'
    data = np.genfromtxt(filename, delimiter=',')
    data = data[data[:, -1] == data[-1, -1]]
    return data.shape[0]


if __name__ == '__main__':
    output_directory = sys.argv[-1]
    layer_directories = glob.glob(output_directory + '/layer_*')
    layers = [int(re.findall(r'\d+', dirname)[-1]) for dirname in layer_directories]
    layer_count = 0
    fracture_data = [[0, 0, 0, 0]]
    for layer in sorted(layers):
        stroke_directories = glob.glob(output_directory + 'layer_' + str(layer) + '/stroke_*')
        strokes = [int(re.findall(r'\d+', dirname)[-1]) for dirname in stroke_directories]
        for stroke in sorted(strokes):
            dir_name = output_directory + 'layer_' + str(layer) + '/stroke_' + str(stroke)
            fracture_count = get_fractured_particles_at_stroke(output_directory, layer, stroke)
            fracture_data.append([layer_count, layer+1, stroke+1, fracture_count])

    fracture_data = np.array(fracture_data)
    with open(output_directory + '/fractured_particles.pkl', 'w') as fracture_pickle:
        pickle.dump(fracture_data, fracture_pickle)
