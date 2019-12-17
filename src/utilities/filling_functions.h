//
// Created by erolsson on 2018-09-10.
//

#ifndef DEMSIM_UTILITIES_H
#define DEMSIM_UTILITIES_H

#include <vector>

#include "vec3.h"

namespace DEM {
    std::vector<Vec3> random_fill_cylinder(double z0, double z1, double cylinder_radius,
            const std::vector<double>& radii);
}

namespace DEM {
    std::vector<Vec3>random_fill_box(double z0, double z1, double box_width,
                                     const std::vector<double>& radii, double bt);
}


#endif //DEMSIM_UTILITIES_H
