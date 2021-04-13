import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import math

mu = 5e-2
variance = 0
sigma = math.sqrt(variance)
x = np.linspace(mu - 0*sigma, mu + 0*sigma, 5000)
plt.plot(x, stats.norm.pdf(x, mu, sigma))
print(x)

np.savetxt('radie.dat', (x))
plt.show()