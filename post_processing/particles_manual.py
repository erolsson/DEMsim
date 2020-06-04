import numpy as np
from mayavi import mlab


def main():
    x = [0., 0., 1.]
    y = [0., 0., 1.]
    z = [-.9, 0.9, 0.]
    r = [1., 1., 1.]
    mlab.figure(size=(1920, 1200), bgcolor=(1., 1., 1.))
    mlab.points3d(x, y, z, 2*np.array(r),
                  color=(184./255, 115./255., 51./255.),
                  resolution=32,
                  scale_factor=1.,
                  scale_mode='scalar')
    mlab.show()


if __name__ == '__main__':
    main()
