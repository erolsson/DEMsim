import numpy as np

r_min = 3
r_max = 10
number = 10000

radii = np.random.uniform(r_min, r_max, 10000)
np.savetxt('../simulations/closed_die_compaction/size_distribution.dat', radii/1000)
