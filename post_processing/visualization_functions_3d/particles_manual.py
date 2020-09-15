import numpy as np
from mayavi import mlab

from visualization_functions_3d.periodic_bc import PeriodicBC


def main():
    x = [0.017302485751200292, 0.020766645421414792]
    y = [-0.017936378251383566, 0.024032796165419458]
    z = [0.07139545774200616, -0.068784743434672585]
    r = [0.0074175868457049586, 0.0035272133491102875]

    mlab.figure(size=(1920, 1200), bgcolor=(1., 1., 1.))
    mlab.points3d(x, y, z, 2*np.array(r),
                  color=(184./255, 115./255., 51./255.),
                  resolution=32,
                  scale_factor=1.,
                  scale_mode='scalar')

    x = [0.017302485751200292, 0.020766645421414792]
    y = [-0.017936378251383566, -0.029078223540450372]
    z = [-0.074872880649644466, -0.068784743434672585]
    r = [0.0074175868457049586, 0.0035272133491102875]

    mlab.points3d(x, y, z, 2*np.array(r),
                  color=(192./255, 192./255., 192./255.),
                  resolution=32,
                  scale_factor=1.,
                  scale_mode='scalar')

    periodic_bc = PeriodicBC("../../results/periodic_bc_test/sim_1/periodic_bc.dou")
    pb_data1 = np.array([-0.0265555, 0.0265555, -0.0265555, 0.0265555, -0.073134169195825313, 0.073134169195825313])
    periodic_bc.draw(pb_data1)

    mlab.show()


if __name__ == '__main__':
    main()
