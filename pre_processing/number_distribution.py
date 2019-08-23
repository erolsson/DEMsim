import numpy as np
from scipy.integrate import trapz


class NumberDistribution:
    def __init__(self, volume_data, num_points=1000, min_size=None, max_size=None):
        # if fraction is given in percent
        self.n = num_points
        if volume_data[-1, 1] > 10.:
            volume_data[:, 1] /= 100

        # If cdf does not start at zero
        if volume_data[0, 1] > 0:
            temp = np.zeros((volume_data.shape[0] + 1, 2))
            temp[1:] = volume_data
            temp[0, 0] = temp[1, 0]/2
            volume_data = temp

        rmin = volume_data[0, 0]
        rmax = volume_data[-1, 0]
        if min_size:
            rmin = min_size
        if max_size:
            rmax = max_size

        self.r = np.linspace(rmin, rmax, num_points)
        self.number_cdf = 0*self.r
        self.volume_cdf = np.interp(self.r, volume_data[:, 0], volume_data[:, 1])

        for i, r_val in enumerate(self.r[1:], start=1):
            r_vec = self.r[0:i]
            f_vec = np.interp(r_vec, volume_data[:, 0], volume_data[:, 1])

            self.number_cdf[i] = f_vec[-1]/r_val**3 + 3*trapz(f_vec/r_vec**4, r_vec)

        self.norm_c = self.number_cdf[-1]
        self.number_cdf /= self.norm_c

    def _create_transformed_cdf(self, cdf, min_size=None, max_size=None):
        # transforming the distribution
        fmin = 0
        fmax = 1
        rmin = self.r[0]
        rmax = self.r[-1]
        if min_size:
            fmin = np.interp(min_size, self.r, cdf)
            rmin = min_size
        if max_size:
            fmax = np.interp(max_size, self.r, cdf)
            rmax = max_size
        r = np.linspace(rmin, rmax, self.n)

        return r, np.interp(r, self.r, (cdf - fmin) / (fmax - fmin))

    def generate(self, number, min_size=None, max_size=None):
        r, transformed_cdf = self._create_transformed_cdf(self.number_cdf, min_size, max_size)
        p = np.random.rand(number)
        return np.interp(p, transformed_cdf, r)

    def number_of_particles_for_volume(self, volume, min_size=None, max_size=None):
        r, transformed_cdf = self._create_transformed_cdf(self.number_cdf, min_size, max_size)
        return volume/(4*np.pi*r[-1]**3/3 - 4*np.pi*trapz(transformed_cdf*r**2, r))


if __name__ == '__main__':
    raw_data = np.genfromtxt('../size_distributions/proctor_distribution.dat', delimiter=',')
    raw_data[:, 0] /= 2
    number_distr = NumberDistribution(raw_data, min_size=1.)
    sizes = number_distr.generate(10000, min_size=4., max_size=8)
    print np.min(sizes), np.max(sizes), np.mean(sizes), np.median(sizes)
    v = 116.43*(25.4*2)**2*np.pi*0.6
    n = int(number_distr.number_of_particles_for_volume(v, min_size=4, max_size=8))
    print n
    radii = number_distr.generate(n)
    print np.sum(4*np.pi*radii**3/3), v*0.6
