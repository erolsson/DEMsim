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
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def plot_mechanical_data_for_contact(directory):
    contact_force = np.genfromtxt(directory + 'contact_force_control', delimiter=',')
    force_particle = np.genfromtxt(directory + 'contact_force_particle', delimiter=',')
    plt.figure(0)

    bt=0.003
    h = contact_force[:, 0]
    force = contact_force[:, 1]
    Fp= force_particle[:,1]
    E_0_part=1/(((1-(0.3**2))/142e9)*2)
    R= 3e-2
    F=4*E_0_part*R**(1/2)*h**(3/2)/3
    plt.plot(h/bt,Fp/E_0_part/R**2, 'b', ms=12 , linestyle='--', label='$F_{particle}$ ')
    plt.plot(h/bt,force/E_0_part/R**2, 'b', ms=12 ,  label='$F_{Tot}$ ')
    plt.plot(h/bt, F/E_0_part/R**2, 'k', ms=12, label='Hertz contact theory')
    plt.xlim(-1, max(h/bt))



def main():
    directory = os.path.expanduser(r'/../../DEMsim/results/viscoelastic/')
    plot_mechanical_data_for_contact(directory)


    plt.figure(0)

    plt.xlabel('$ h/ b{t}$ [-]')
    plt.ylabel('Contact Force [-] ')
    plt.legend(loc='best')


    plt.show()


if __name__ == '__main__':
    main()


