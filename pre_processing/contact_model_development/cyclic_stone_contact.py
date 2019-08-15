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


class StoneContact:
    def __init__(self, R0, E, v):
        self.R0 = R0
        self.E0 = E/(1 - v**2)

        self.n = 1.5
        self.m = 1.7

        self._h = 0
        self._F = 0

        # State variables
        self._hmax = 0
        self._hp = 0
        self._hl = 0

        # stiffnesses
        self.kp = 4./3*self.E0*self.R0**0.5
        self.ke = 4./3*self.E0*self.R0**0.5/0.95
        self.ku = 4./3*self.E0*self.R0**(2-self.m)/0.95
        self.kl = self.ke

        self._hs = 10e-3
        self._F0 = 100
        self.h1 = self._hs - (3*self._F0/4/self.E0/self.R0**0.5)**(2./3)

    def update(self, dh):
        self._h += dh
        if self._h >= self._hmax:
            if self._h < self._hs:
                self._F = self._F0/self._hs*self._h
            else:
                self._F = self.kp*(self._h - self.h1)**self.n
            self._hmax = self._h
            self._hp = self._hmax - (self._F/self.ke)**(1./self.n)
            self.ku = self._F/(self._h - self._hp)**self.m
        else:
            if dh >= 0:
                print "Reloading", dh, self._h - 10e-3, self.kl
                h1 = self._hl + (self._F0/self.kl)**(2./3)
                if self._h < h1 and self._h - h1 + self._hs > 0:
                    self._F = self._F0/self._hs*(self._h - h1 + self._hs)
                else:
                    self._F = self.kl*(self._h - self._hl)**self.n
                self.ku = self._F/(self._h - self._hp)**self.m
            elif self._h - self._hp > 0:
                self._F = self.ku*(self._h - self._hp)**self.m
                A = (self._F/(self.kp*(self._hmax - self.h1)**self.n))**(1./self.n)
                print "Unloading", A, self._hl, self.ku
                self._hl = (self._h - A*self._hmax)/(1-A)
                self.kl = self.kp*(self._hmax - self.h1)**self.n/(self._hmax - self._hl)**self.n
            else:
                self._F = 0

    def get_force(self):
        return self._F

    def get_plastic_indentation_depth(self):
        return self._hp

    def get_hl(self):
        return self._hl


def get_test_results(first_cycle=1, second_cycle=2):
    cycles_to_idx = {0: 0, 1: 1179, 2: 2343}
    data = np.genfromtxt(os.path.expanduser('~/work/stone_materials/exp_results/arlanda/2_0kN_test1.csv'),
                         delimiter=";", skip_header=2)
    return data[cycles_to_idx[first_cycle-1]:cycles_to_idx[second_cycle], [5, 4]]


if __name__ == '__main__':
    test_data = get_test_results(second_cycle=1)
    h = test_data[:, 0]
    plt.plot(h, test_data[:, 1]*1000, 'b', lw=2)
    h += 10e-3
    h[0] = 0
    print np.max(h)
    R = 6.25
    F = 0*h
    hp = 0*h
    hl = 0*h
    contact = StoneContact(R, 65E3, 0.15)
    for i in range(1, h.shape[0]):
        contact.update(h[i] - h[i-1])
        F[i] = contact.get_force()
        hp[i] = contact.get_plastic_indentation_depth()
        hl[i] = contact.get_hl()
    diff_h = 0*h
    diff_h[1:] = np.diff(h)
    plt.plot(h-10e-3, F, 'k', lw=2)

    plt.figure(2)
    plt.plot(h, hp)

    plt.figure(2)
    plt.plot(h, hl)
    plt.show()
