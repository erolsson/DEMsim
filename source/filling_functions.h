//
// Created by erolsson on 2018-09-10.
//

#ifndef DEMSIM_UTILITIES_H
#define DEMSIM_UTILITIES_H

#include <vector>

#include "vec3.h"

namespace DEM {
    std::vector<Vec3> random_fill_cylinder(double z0, double z1, double R, std::vector<double> radii);
    bool check_overlaps(const Vec3& point, double radius, const std::vector<Vec3>& pos,
                       const std::vector<double>& radii);
}

#endif //DEMSIM_UTILITIES_H
