import numpy as np

from number_distribution import NumberDistribution

rmax = 8.
rmin = 0.063/2

cdf_data = np.zeros((1000, 2))
cdf_data[:, 0] = np.linspace(rmin, rmax, 1000)
cdf_data[:, 1] = (cdf_data[:, 0]/rmax)**0.5

distr = NumberDistribution(cdf_data)
v = 116.43*(25.4*2)**2*np.pi*0.8
print distr.number_of_particles_for_volume(v, min_size=1)
sizes = distr.generate(30000, min_size=1., max_size=8)
np.savetxt('../simulations/proctor_test/fuller.dat', sizes)
