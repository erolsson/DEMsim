import os

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import fmin
from scipy.stats import cumfreq
from scipy.stats import truncnorm
from scipy.stats import rv_continuous
from scipy.integrate import quad


class NumberDistribution(rv_continuous):

    def __init__(self, volume_data, base_distribution=truncnorm):
        self.base_distribution = base_distribution
        self.sizes = volume_data[:, 0]
        self.values = volume_data[:, -1]

        def residual(par):
            bounds = (self.sizes[0] - par[0])/par[1], (self.sizes[-1] - par[0])/par[1]
            return sum((self.values - self.base_distribution.cdf(self.sizes, bounds[0], bounds[1],
                                                                 loc=par[0], scale=par[1]))**2)

        self.parameters = fmin(residual, [3e-3, 10e-3])
        self.a_norm = (self.sizes[0] - self.parameters[0])/self.parameters[1]
        self.b_norm = (self.sizes[-1] - self.parameters[0])/self.parameters[1]
        super(NumberDistribution, self).__init__(a=self.a_norm, b=self.b_norm)

        self.B = quad(self._f, self.sizes[0], self.sizes[-1])[0]
        num_cdf_points = 100
        self._cdf_sizes = np.linspace(self.sizes[0], self.sizes[-1], num_cdf_points)
        self._cdf_values = np.zeros(num_cdf_points)
        for i, r in enumerate(self._cdf_sizes):
            self._cdf_values[i] = quad(self._f, self.sizes[0], r, args=(self.B, ))[0]

    def _f(self, x, scale=1.):
        return self.base_distribution.pdf(x, self.a_norm, self.b_norm,
                                          loc=self.parameters[0], scale=self.parameters[1])/x**3/scale

    def _cdf(self, x, *args, **kwargs):
        cdf_values = np.interp(x, self._cdf_sizes, self._cdf_values)
        cdf_values[x <= self.sizes[0]] = 0
        cdf_values[x >= self.sizes[-1]] = 1
        return cdf_values

    def volume_cdf(self, x):
        return self.base_distribution.cdf(x, self.a_norm, self.b_norm,
                                          loc=self.parameters[0], scale=self.parameters[1])


def generate_number_distribution_from_gradation(gradation_file_name, number_of_particles, size_factor=1.,
                                                plot_graphs=False):
    data = np.genfromtxt(gradation_file_name,  delimiter=',')
    data[:, 0] *= size_factor

    if data[-1, 1] > 10:
        data[:, 1] /= 100

    distribution = NumberDistribution(data)
    r = distribution.rvs(size=number_of_particles)
    if plot_graphs:
        plt.figure(0)
        plt.plot(data[:, 0], data[:, 1], '-ko', lw=2, label='Raw data')

        size = np.linspace(data[0, 0]-1e-3, data[-1, 0]+1e-3, 1000)
        plt.plot(size, distribution.volume_cdf(size), '--k', lw=2, label='Normal fit')
        plt.legend(loc='best')

        plt.figure(1)
        plt.plot(size, distribution.cdf(size), 'r', lw=2, label='Number distribution')
        plt.legend(loc='best')

        res = cumfreq(r, defaultreallimits=(np.min(r), np.max(r)))
        x = res.lowerlimit + np.linspace(0, res.binsize*res.cumcount.size, res.cumcount.size)
        plt.bar(x, res.cumcount/number_of_particles, width=res.binsize)

        volume = 4*np.pi*r**3/3
        volume = np.sort(volume)
        cdf_vol = np.cumsum(volume)/np.sum(volume)
        plt.figure(0)
        plt.plot(np.sort(r), cdf_vol)
    return r


if __name__ == '__main__':
    gradation_file = os.path.expanduser('~/work/railway_ballast/experimental_results/'
                                        'sun_et_al_16/size_distribution.dat')
    radii_file_name = '../simulations/cyclic_triaxial/size_distribution.dat'

    directory = os.path.dirname(os.path.abspath(radii_file_name))
    if not os.path.isdir(directory):
        os.makedirs(directory)
    radii = generate_number_distribution_from_gradation(gradation_file, 1000, 1.5e-4, plot_graphs=True)
    print radii
    np.savetxt(radii_file_name, radii)
    plt.show()
