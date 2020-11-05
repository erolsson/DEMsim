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


def main():
    directory = os.path.expanduser("~/DEMsim/results/Battery")
    periodic_box = np.genfromtxt(directory + "/periodic_bc.dou", delimiter=',')
    idx = periodic_box[:, 1] != periodic_box[0, 1]
    force_fabric_tensor = np.genfromtxt(directory + "/force_fabric_tensor.dou", delimiter=',')
    particles = np.genfromtxt(directory + "/particles/particles_0.01.dou", delimiter=',')
    r = particles[:, 7]
    print(0.005/r)
    print(r+0.005/2)
    vp = sum(4*np.pi*r**3/3)
    time = force_fabric_tensor[idx, 0]
    lx = 2*periodic_box[idx, 2]
    ly = 2*periodic_box[idx, 4]
    lz = 0.540713
    print(vp/lx[0]/ly[0]/lz)
    plt.figure(10)
    plt.plot((lx - lx[0]))

    exx = (lx - lx[0])/lx[0]
    stresses = force_fabric_tensor[idx, 1:]/np.outer(lx, np.ones(9))/ly[0]/lz
    max_stress_idx = np.argmax(stresses[:, 0])
    print(time[max_stress_idx], stresses[max_stress_idx, 0]/1e6)
    plt.figure(0)
    plt.plot(stresses[:, 0]/1e6)
    plt.figure(1)
    plt.plot(stresses[:, 4]/1e6)
    plt.figure(2)
    plt.plot(stresses[:, 8]/1e6)

    plt.figure(3)
    plt.plot(exx)

    de = np.diff(exx)
    print(np.mean(np.diff(-stresses[:, 0])[de != 0]/de[de != 0])/1e6)

    plt.show()


if __name__ == '__main__':
    main()
