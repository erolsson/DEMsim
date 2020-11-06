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
    R = 0.0425
    E = 1.7e9
    nu = 0.25
    h = np.linspace(0, 0.1*R, 1000)
    F = 4./3*E/(1-nu**2)*np.sqrt(R)*h**1.5
    plt.plot(h, F)

    F = (1-nu)/(1+nu)/(1-2*nu)/2/0.1/R*E*0.3*R*R*h
    plt.plot(h, F)
    plt.show()

o
if __name__ == '__main__':
    main()
