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


def prony(time, E0, a, tau):
    E = 0*time + 1.
    for ai, ti in zip(a, tau):
        E -= ai*(1-np.exp(-time/ti))
    return E0*E


def main():
    t = np.linspace(0, 100, 1000)
    plt.plot(t, 0.117 + 0.065*np.exp(-t/211) + 0.057*np.exp(-t/4807), 'k', lw=3)
    E0 = 0.117 + 0.065 + 0.057
    plt.plot(t, prony(t, E0, [0.065/E0, 0.057/E0], [211, 4807]))

    plt.plot(t, prony(t, 2*E0, [0.5, 0.065/E0/2, 0.057/E0/2], [1e-3, 211, 4807]))

    plt.show()



if __name__ == '__main__':
    main()
