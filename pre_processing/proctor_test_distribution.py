import numpy as np
from number_distribution import NumberDistribution

raw_data = np.genfromtxt('../size_distributions/proctor_distribution.dat', delimiter=',')
raw_data[:, 0] /= 2
number_distr = NumberDistribution(raw_data, min_size=1.)
sizes = number_distr.generate(10000, min_size=4., max_size=8)
np.savetxt('../simulations/proctor_test/8-16mm.dat', sizes/1000)
