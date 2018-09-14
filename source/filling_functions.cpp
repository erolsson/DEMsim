//
// Created by erolsson on 2018-09-10.
//

#include "filling_functions.h"

#include <random>

std::vector<DEM::Vec3> DEM::random_fill_cylinder(double z0, double z1, double R, std::vector<double> radii)
{
    std::vector<Vec3> particle_positions;
    std::random_device random_device;
    std::default_random_engine rand_endine(random_device());
    std::uniform_real_distribution<double> dist_r(-R, R);
    std::uniform_real_distribution<double> dist_z(z0, z1);
    for (auto particle_radius : radii) {
        bool overlapping = true;
        Vec3 position {};
        while(overlapping) {
            position.x = dist_r(rand_endine);
            position.y = dist_r(rand_endine);
            position.z = dist_z(rand_endine);

            //Check if a particle at the chosen position overlaps with an other
            if (position.x*position.x+position.y*position.y < (R-particle_radius)*(R-particle_radius)) {
                overlapping = check_overlaps(position, particle_radius, particle_positions, radii);
            }
        }

        particle_positions.push_back(position);

    }
    return particle_positions;
}

bool DEM::check_overlaps(const DEM::Vec3& point, double radius, const std::vector<DEM::Vec3>& particle_positions,
                        const std::vector<double>& radii)
{
    for (unsigned i = 0; i != particle_positions.size(); ++i) {
        Vec3 position_i = particle_positions[i];
        if ((point-position_i).length() < radii[i] + radius) {
            return true;
        }
    }
    return false;
}
