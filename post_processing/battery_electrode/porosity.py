import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def calculate_occupied_volume(time, simulation_directory):
    particles = np.genfromtxt(simulation_directory + 'particles/particles_' + str(time) + '.dou', delimiter=',')
    r = particles[:, 7]
    v_part = np.sum(4*np.pi*r**3/3)

    contacts = np.genfromtxt(simulation_directory + 'contacts/contacts_' + str(time) + '.dou', delimiter=',')
    binder_contacts = contacts[contacts[:, -2] == 1, :]
    max_overlap = binder_contacts[:, 10]
    bt = binder_contacts[:, -1]
    activated_binder_contacts = binder_contacts[max_overlap > bt, :]
    br = activated_binder_contacts[:, -1]
    # assuming a cylindric binder volume with radius
    rb = 0.3*0.03
    v_binder = np.sum(rb*rb*np.pi*bt)

    return v_binder + v_part


def main():
    time_for_test = 14.7
    directory = os.path.expanduser('/scratch/users/elaheh/DEMsim/results/viscoelastic/cubic_box-3200_h90-bt1-tryck')
    box_edges = np.genfromtxt(directory + 'periodic_bc.dou', delimiter=',')
    time = box_edges[:, 0]
    box_side = 2*box_edges[:, 1]
    surface_positions = np.genfromtxt(directory + 'surface_positions.dou', delimiter=',')
    box_height = surface_positions[surface_positions[:, -1] == time_for_test, -2]
    volume = box_side**2*box_height
    porosity = 1 - calculate_occupied_volume(time_for_test, directory)/volume
    plt.plot(time, porosity)

    plt.show()


if __name__ == '__main__':
    main()
