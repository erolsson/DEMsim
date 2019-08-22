import os

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import trapz

raw_data = np.genfromtxt('../size_distributions/proctor_distribution.dat', delimiter=',')
r = np.zeros(raw_data.shape[0] + 1)
Fv = 0*r
Fv[1:] = raw_data[:, 1] / 100
r[1:] = raw_data[:, 0]/2
r[0] = r[1]/2
plt.semilogx(r, Fv)

radii = np.exp(np.linspace(np.log(r[0]), np.log(r[-1]), 1000))
Fn = 0*radii
for i, r_val in enumerate(radii[1:], start=1):
    points = np.count_nonzero(r < r_val) + 1
    r_vec = radii[0:i]
    f_vec = np.interp(r_vec, r, Fv)

    Fn[i] = f_vec[-1]/r_val**3 + 3*trapz(f_vec/r_vec**4, r_vec)

Fn /= Fn[-1]
plt.semilogx(radii, Fn)
rmin = 4
rmax = 8
Fn_new = (Fn - np.interp(rmin, radii, Fn))/(np.interp(rmax, radii, Fn)-np.interp(rmin, radii, Fn))
plt.figure(2)
print Fn - np.interp(rmin, radii, Fn)
plt.plot(radii, Fn_new, '-*')
plt.xlim(rmin, rmax)
plt.ylim(0, 1)

p = np.random.rand(10000)
particle_radii = np.interp(p, Fn_new, radii)
print np.max(particle_radii), np.min(particle_radii), np.mean(particle_radii), np.median(particle_radii)

plt.show()
