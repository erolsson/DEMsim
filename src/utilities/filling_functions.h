//
// Created by erolsson on 2018-09-10.
//

#ifndef DEMSIM_UTILITIES_H
#define DEMSIM_UTILITIES_H

#include <vector>

#include "vec3.h"

namespace DEM {
    std::vector<Vec3> random_fill_cylinder(double z0, double z1, double cylinder_radius,
            const std::vector<double>& radii, double distance_between_objects=0);


    std::vector<Vec3>random_fill_box(double z0, double z1, double box_width,
                                     const std::vector<double>& radii, double bt);
    std::vector<Vec3>random_fill_box(double x1, double x2, double y1, double y2, double z1, double z2,
                                     const std::vector<double>& radii);


    bool check_overlaps(const Vec3 &point, double radius, const std::vector<DEM::Vec3> &particle_positions,
                    const std::vector<double> &radii);
}

#endif //DEMSIM_UTILITIES_H
