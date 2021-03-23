import glob
import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib
from multiprocesser.multiprocesser import multi_processer

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})


def get_contact_output_times(directory):
    contact_files = glob.glob(directory + '/contacts/contacts*.dou')
    contact_times = []
    for contact_file in contact_files:
        contact_time = contact_file.replace(directory + '/contacts/contacts_', '').replace('.dou', '')
        contact_times.append(float(contact_time))
    return sorted(contact_times)


def calculate_particle_contacts(directory, time):
    contact_data = np.genfromtxt(directory + '/contacts/contacts_' + str(time) + '.dou', delimiter=',')
    particle_contacts = contact_data[np.logical_and(contact_data[:, -2] == 0, contact_data[:, 8] > 0), 1]
    surface_idx = {0, 1}
    values = particle_contacts*0 + 2  # Each contact belong to two particles
    for surface_id in surface_idx:
        values[particle_contacts == surface_id] -= 1    # ... except if it is a surface contact
    return np.sum(values)


def plot_mechanical_data_for_simulation(directory):
    periodic_bc = np.genfromtxt(directory + 'periodic_bc.dou', delimiter=',')
    surface_positions = np.genfromtxt(directory + 'surface_positions.dou', delimiter=',')
    force_fabric_tensor = np.genfromtxt(directory + 'force_fabric_tensor.dou', delimiter=',')
    time = surface_positions[:, -1]
    d = 2*periodic_bc[:, 2]
    d0 = d[0]
    w = 2*periodic_bc[:, 4]

    # Finding the point where the compaction force has decreased to zero after its maximum value
    # The thickness at this point will be the thickness of the electrode t0
    t0 = 0.880421

    t_start = time[d == d0][-1]
    volume = (d*w*t0)[time > t_start]
    sxx = force_fabric_tensor[time > t_start, 1]/volume
    syy = force_fabric_tensor[time > t_start, 5]/volume
    szz = force_fabric_tensor[time > t_start, 9]/volume

    sxx *= -1
    syy *= -1
    szz *= -1

    linear_strain = (d[time > t_start] - d0)/d0

    plt.figure(0)
    if sxx[-1] > 0:
        plt.plot(linear_strain, sxx/1e6, 'r', lw=2, label=r'$\sigma_{xx}$')
        plt.plot(linear_strain, syy/1e6, 'b', lw=2, label=r'$\sigma_{yy}$')
        plt.plot(linear_strain, szz/1e6, 'g', lw=2, label=r'$\sigma_{zz}$')
    else:
        plt.plot(linear_strain, sxx/1e6, 'r', lw=2)
        plt.plot(linear_strain, syy/1e6, 'b', lw=2)
        plt.plot(linear_strain, szz/1e6, 'g', lw=2)

    idx = [1]
    delta_strain = 3e-3
    while len(idx):
        # This finds the turning points
        bool_arr = np.logical_and(np.diff(np.abs(linear_strain)) < 0,
                                  np.abs(linear_strain[1:]) > np.abs(linear_strain[idx[0]]))

        idx = np.where(bool_arr)[0]
        if len(idx):
            e0 = linear_strain[idx[0]]
            e1 = e0 + np.abs(e0)/e0*delta_strain           # increasing the magnitude of the strain with delta_strain
            idx1 = np.argmin(np.abs(linear_strain - e1))   # Finding the index of the strain point closest to e1
            e1 = linear_strain[idx1]                       # Grabbing the correct e1
            dsxx = sxx[idx1] - sxx[idx[0]]
            dsyy = syy[idx1] - syy[idx[0]]
            dszz = szz[idx1] - szz[idx[0]]

            dexx = e1 - e0
            if dexx != 0:

                v = dsyy/(dsxx + dszz)
                plt.figure(1)
                plt.plot(e0, v, 'kx', ms=12, mew=2)

                E = (dsxx - v*(dsyy + dszz))/dexx
                plt.figure(2)
                plt.plot(e0, E/1e9, 'kx', ms=12, mew=2)

    contact_times = get_contact_output_times(directory)
    particles = np.genfromtxt(directory + '/particles/particles_' + str(contact_times[0]) + '.dou',
                              delimiter=',')
    no_particles = particles.shape[0]

    job_list = []
    contact_times = np.array(contact_times)
    contact_times = contact_times[contact_times > t_start]
    for c_time in contact_times:
        job_list.append((calculate_particle_contacts, [directory, c_time], {}))
    particle_contacts = np.array(multi_processer(job_list, delay=0., timeout=3600))
    particle_contact_per_particle = particle_contacts/no_particles
    e = np.interp(np.array(contact_times), time[time > t_start], linear_strain)

    plt.figure(3)
    plt.plot(e, particle_contact_per_particle, 'k', lw=2)


def main():
    directory = os.path.expanduser('C:/DEMsim/results/E34br05bt003/20210323compression/')
    plot_mechanical_data_for_simulation(directory)

    directory = os.path.expanduser('C:/DEMsim/results/tension-E34bt01Rbr05/unload_restart_file/')
    plot_mechanical_data_for_simulation(directory)

    plt.figure(0)
    plt.xlabel('Strain [-]')
    plt.ylabel('Stress [MPa]')
    plt.legend(loc='best')

    plt.figure(1)
    plt.xlabel('Strain [-]')
    plt.ylabel(r'$\nu$ [-]')

    plt.figure(2)
    plt.xlabel('Strain [-]')
    plt.ylabel('$E$ [GPa]')

    plt.figure(3)
    plt.xlabel('Strain [-]')
    plt.ylabel('Particle contacts / Particle [-]')

    plt.show()


if __name__ == '__main__':
    main()