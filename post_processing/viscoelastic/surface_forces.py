import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from pandas import read_table

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def pressure_box(data_directory):
    force_data = np.genfromtxt(data_directory + '/surface_forces.dou', delimiter=',')
    with open(data_directory + '/surface_positions.dou', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
        id_index = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[index + 1][5:] for index in id_index]

    box_idx = [i for i, surface_type in enumerate(surface_types) if surface_type == 'PointSurface']

    surface_force = force_data[:, box_idx[1] * 5 + 1]
    surface_force_side = force_data[:, box_idx[4] * 5 + 2]

    p = surface_force / (0.15 * 0.15)
    p_side = surface_force_side
    return p


def indentation(data_directory):
    indentation_data = np.genfromtxt(data_directory + '/surface_positions.dou', delimiter=',')
    with open(data_directory + '/surface_positions.dou', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
        id_index = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[index + 1][5:] for index in id_index]

    box_idx = [i for i, surface_type in enumerate(surface_types) if surface_type == 'PointSurface']
    print(box_idx)
    b = box_idx[3] * 5 + 3
    print(b)
    surface_indentation = indentation_data[:, box_idx[1] * 15 + 5]
    surface_indentation_side = indentation_data[:, box_idx[4] * 15+3 ]

    indent = surface_indentation
    side = surface_indentation_side
    return indent


def binder(data_directory):

    rad = np.genfromtxt(data_directory + '/radie.dat', delimiter=',')
    binder_volume = 0
    n = 0
    bt = 0.5e-2
    for i in range(len(rad)):
        if binder_volume < 1.5 * 0.1877:
            binder_volume += (4 * 3.14 / 3) * ((rad[i] + bt) ** 3 - (rad[i]) ** 3)
            n += 1
    binder_surface = 0
    for i in range(n):
        binder_surface += 4*3.14*((rad[i]+bt)**2)
    particle_surface = 0
    for i in range(2553):
        particle_surface += 4*3.14*((rad[i])**2)
    contact_percent = binder_surface/particle_surface
    return contact_percent


if __name__ == '__main__':
    simulation_directory = 'C:/DEMsim/results/viscoelastic'
    pressure = pressure_box(simulation_directory)
    h = indentation(simulation_directory)
    h_side = indentation(simulation_directory)
    print (h_side)
    pressure_side = pressure_box(simulation_directory)
    # a use to draw relaxation
    a = np.linspace(0., 62.59, 6259)
    print(pressure_side)
    just_particle_volume = 1.33587
    box_width = 1.5
    #Prorosity= (1-((just_particle_volume)/ (box_width*box_width*h+0.075*(box_width*box_width*h))))*100;
    Prorosity = (1-(just_particle_volume/ (box_width*box_width*h)))*100;
    #print(h)
    #print(pressure)
    #ut = [i for i in np.arange(0, 7.64, 0.00001)]

    #contact_percent= binder(simulation_directory)
    #print(contact_percent)



    fig1, ax1 = plt.subplots()
    ax1.plot(pressure,h)
    ax1.set_title("Calendering")
    ax1.set_xlabel("pressure [Pa] ")
    ax1.set_ylabel("prosity  ")
    plt.show()

    mu, sigma = 0.05, 0.01 # mean and standard deviation
    s = np.random.normal(mu, sigma, 9000)
    np.savetxt('radius.out', s, delimiter=',')













#x = np.arange(5)
#y = np.exp(x)
#fig1, ax1 = plt.subplots()
#ax1.plot(x, y)
#ax1.set_title("Axis 1 title")
#ax1.set_xlabel("X-label for axis 1")

#z = np.sin(x)
#fig2, (ax2, ax3) = plt.subplots(nrows=2, ncols=1) # two axes on figure
#ax2.plot(x, z)
#ax3.plot(x, -z)

#w = np.cos(x)
#ax1.plot(x, w)

