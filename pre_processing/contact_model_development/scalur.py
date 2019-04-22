from collections import namedtuple

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


class Contact2D:
    def __init__(self, particle1, particle2, overlap, mu):
        self.R0 = 1./(1/particle1.R + 1/particle2.R)
        h = overlap*self.R0
        self.mu = mu
        self.a = np.sqrt(self.R0*h)
        G1 = particle1.E/2/(1 + particle1.v)
        G2 = particle2.E/2/(1 + particle2.v)
        self.km = 8/((2 - particle1.v)/G1 + (2 - particle2.v)/G2)
        self.kh = 1./((1 - particle1.v**2)/particle1.E + (1 - particle2.v**2)/particle2.E)*4./3*np.sqrt(self.R0)

        self.Q1 = 0
        self.Q2 = 0
        self.Q = 0
        self.old_dt = 0
        self.delta = 0
        self.F = 4./3*self.kh*np.sqrt(R0)*overlap**1.5
        self.h = overlap
        self.load_var = 0

    def calc_tangential_force(self, dt, dh, mu):
        f1 = 4./3*self.kh*np.sqrt(R0)*(self.h+dh)**1.5
        df = f1 - self.F
        self.F = f1
        self.h += dh

        if dt < 0 < self.old_dt:   # sign shift
            self.Q1 = self.Q

        if self.old_dt < 0 < dt:   # sign shift
            self.Q2 = self.Q
        sgn = 1
        if dt < 0:
            sgn = -1

        self.Q1 += mu*df
        self.Q2 -= mu*df

        if dh <= 0 or self.load_var > 0:
            if dt > 0 and self.Q1 == 0. and self.Q2 == 0.:
                q = (1 - (self.Q + mu*df)/mu/self.F)
            elif dt < 0:
                q = (1 - (self.Q1 - self.Q + 2*mu*df)/(mu*self.F*2))
            else:
                q = (1 - (self.Q - self.Q2 + 2*mu*df)/mu/self.F/2)
            if dh <= 0:
                self.load_var = 0
            if q > 0:
                q = q**(1./3)
            else:
                q = 0
        else:
            self.load_var += self.km*np.sqrt(self.h*self.R0)*abs(dt) - mu*df
            q = 1

        kt = self.km*np.sqrt(self.h*self.R0)*q
        if dh != 0.:
            kt += sgn*mu*(1-q)*df/dt
        self.old_dt = dt
        self.Q += kt*dt
        self.delta += dt
        print "2d", kt, sgn*mu*(1-q)*df/dt, q


class Contact3D:
    def __init__(self, particle1, particle2, overlap, mu):
        self.R0 = 1./(1/particle1.R + 1/particle2.R)
        h = overlap*self.R0
        self.mu = mu
        self.a = np.sqrt(self.R0*h)
        G1 = particle1.E/2/(1 + particle1.v)
        G2 = particle2.E/2/(1 + particle2.v)
        self.km = 8/((2 - particle1.v)/G1 + (2 - particle2.v)/G2)
        self.kh = 1./((1 - particle1.v**2)/particle1.E + (1 - particle2.v**2)/particle2.E)*4./3*np.sqrt(self.R0)

        self.Q1 = np.zeros(3)
        self.Q = np.zeros(3)
        self.old_dt = np.zeros(3)
        self.delta = np.zeros(3)
        self.F = 4./3*self.kh*np.sqrt(R0)*overlap**1.5
        self.h = overlap
        self.load_var_tan = np.zeros(3)
        self.load_var_norm = 0

    def calc_tangential_force(self, dt, dh, mu):
        f1 = 4./3*self.kh*np.sqrt(R0)*(self.h+dh)**1.5
        df = f1 - self.F
        self.F = f1
        self.h += dh

        if np.dot(dt, self.old_dt) < 0:   # sign shift
            self.Q1 = np.copy(self.Q)

        sgn = 1

        if np.linalg.norm(self.Q1):
            self.Q1 += mu*df*self.Q1/np.linalg.norm(self.Q1)

        if dh <= 0 or np.linalg.norm(self.load_var_tan) - self.load_var_norm > 0:
            if np.dot(dt, self.old_dt) > 0 and np.linalg.norm(self.Q1) == 0.:
                q = (1 - (np.linalg.norm(self.Q) + mu*df)/mu/self.F)
            else:
                q = (1 - (np.linalg.norm(self.Q1 - self.Q) + 2*mu*df)/(mu*self.F*2))
            if dh <= 0:
                self.load_var_tan *= 0
                self.load_var_norm = 0
            if q > 0:
                q = q**(1./3)
            else:
                q = 0
        else:
            self.load_var_tan += self.km*np.sqrt(self.h*self.R0)*dt
            self.load_var_norm += mu*df
            q = 1

        kt = self.km*np.sqrt(self.h*self.R0)*q
        if dh != 0.:
            kt += sgn*mu*(1-q)*df/np.linalg.norm(dt)
        self.old_dt = dt
        self.Q += kt*dt
        self.delta += dt
        print "3d", kt, sgn*mu*(1-q)*df/np.linalg.norm(dt), q


Particle = namedtuple('Particle', ['R', 'E', 'v'])

N = 1000

p1 = Particle(R=6.25, E=65E3, v=0.15)
p2 = Particle(R=6.25, E=65E3, v=0.15)
R0 = 1./(1/p1.R + 1/p2.R)
c2d = Contact2D(p1, p2, 8e-2/R0, 0.3)
c3d = Contact3D(p1, p2, 8e-2/R0, 0.3)

dmax = 0.01*R0
delta2d = np.zeros(5*N)
ft2d = np.zeros(5*N)
delta3d = np.zeros((5*N, 3))
ft3d = np.zeros((5*N, 3))

inc = dmax/N
for i in range(0, N):
    c2d.calc_tangential_force(inc, dh=0, mu=0.3)
    c3d.calc_tangential_force(np.array([inc, 0, 0]), dh=0, mu=0.3)
    delta2d[i] = c2d.delta
    ft2d[i] = c2d.Q

    delta3d[i, :] = c3d.delta
    ft3d[i, :] = c3d.Q
print "Unloading"
for i in range(0, 2*N):
    dtx = -inc
    dh = 1E-5
    c2d.calc_tangential_force(dtx, dh=dh, mu=0.3)
    c3d.calc_tangential_force(np.array([dtx, 0, 0]), dh=dh, mu=0.3)
    delta2d[i + N] = c2d.delta
    ft2d[i + N] = c2d.Q

    delta3d[i + N, :] = c3d.delta
    ft3d[i + N, :] = c3d.Q

for i in range(0, 2*N):
    dtx = -inc
    dh = 1E-5
    c2d.calc_tangential_force(inc, dh=dh, mu=0.3)
    c3d.calc_tangential_force(np.array([inc, 0, 0]), dh=0, mu=0.3)
    delta2d[i + 3*N] = c2d.delta
    ft2d[i + 3*N] = c2d.Q

    delta3d[i + 3*N, :] = c3d.delta
    ft3d[i + 3*N, :] = c3d.Q


plt.plot(delta2d, ft2d)
plt.plot(delta3d[:, 0], ft3d[:, 0], 'r', lw=2)
print ft3d
plt.show()
