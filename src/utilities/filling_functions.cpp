//
// Created by erolsson on 2018-09-10.
//

#include "filling_functions.h"

#include <array>
#include <vector>
#include <random>

#include "../engine/periodic_bc_handler.h"

std::vector<DEM::Vec3> DEM::random_fill_cylinder(double z0, double z1, double cylinder_radius,
        const std::vector<double>& radii, double distance_between_objects)
{
    std::vector<Vec3> particle_positions;
    std::random_device random_device;
    std::default_random_engine rand_engine(random_device());
    for (auto r : radii) {
        std::uniform_real_distribution<double> dist_r(-cylinder_radius+r + distance_between_objects,
                                                      cylinder_radius-r - distance_between_objects);
        std::uniform_real_distribution<double> dist_z(z0+r + distance_between_objects, z1-r - distance_between_objects);
        bool overlapping = true;
        Vec3 position {};
        while(overlapping) {
            position.x() = dist_r(rand_engine);
            position.y() = dist_r(rand_engine);
            position.z() = dist_z(rand_engine);

            //Check if a particle at the chosen position overlaps with an other
            if (position.x()*position.x()+position.y()*position.y() <
                   (cylinder_radius - r - distance_between_objects)*(cylinder_radius - r - distance_between_objects)) {
                overlapping = check_overlaps(position, r + distance_between_objects, particle_positions, radii);
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
        if ((radii[i] + radius - (point-position_i).length()) > 0.) {
            return true;
        }
    }
    return false;
}



std::vector<DEM::Vec3> DEM::random_fill_box(double x1, double x2, double y1, double y2, double z1, double z2,
                                       const std::vector<double>& radii, double delta) {
    std::vector<Vec3> particle_positions;
    std::random_device random_device;
    std::default_random_engine rand_engine(0); //(random_device());
    for (auto r : radii) {
        std::uniform_real_distribution<double> dist_x(x1 + r + delta, x2 - r - delta);
        std::uniform_real_distribution<double> dist_y(y1 + r + delta, y2 - r - delta);
        std::uniform_real_distribution<double> dist_z(z1 + r + delta, z2 - r - delta);

        bool overlapping = true;
        Vec3 position {};
        while(overlapping) {
            position.x() = dist_x(rand_engine);
            position.y() = dist_y(rand_engine);
            position.z() = dist_z(rand_engine);
            //Check if a particle at the chosen position overlaps with an other
            overlapping = check_overlaps(position, r + delta, particle_positions, radii);
        }

        particle_positions.push_back(position);
    }
    return particle_positions;
}

std::vector<DEM::Vec3> DEM::random_fill_box_periodic(double x1, double x2, double y1, double y2, double z1, double z2,
                                            const std::vector<double>& radii, double delta,
                                            const std::string& directions) {
    std::vector<Vec3> particle_positions;
    std::vector<Vec3> mirror_positions;
    std::random_device random_device;
    std::default_random_engine rand_engine(0); //(random_device());
    std::array<bool, 3> periodic_directions = {false, false, false};
    std::array<Interval, 3> periodic_boundaries;
    for(const auto& dir: directions) {
        if (dir == 'x') {
            periodic_directions[0] = true;
            periodic_boundaries[0].min = x1;
            periodic_boundaries[0].max = x2;
        }
        if (dir == 'y') {
            periodic_directions[1] = true;
            periodic_boundaries[1].min = y1;
            periodic_boundaries[1].max = y2;
        }
        if (dir == 'z') {
            periodic_directions[2] = true;
            periodic_boundaries[2].min = z1;
            periodic_boundaries[2].max = z2;
        }
    }
    for (auto r : radii) {
        std::uniform_real_distribution<double> dist_x(x1 + delta + !periodic_directions[0]*r,
                                                      x2 - delta - !periodic_directions[0]*r);
        std::uniform_real_distribution<double> dist_y(y1 + delta + !periodic_directions[1]*r,
                                                      y2 - delta - !periodic_directions[1]*r);
        std::uniform_real_distribution<double> dist_z(z1 + delta + !periodic_directions[2]*r,
                                                      z2 - delta - !periodic_directions[2]*r);

        bool overlapping = true;
        Vec3 position {};
        std::vector<Vec3> particle_mirror_positions {};
        while(overlapping) {
            particle_mirror_positions.clear();
            position.x() = dist_x(rand_engine);
            position.y() = dist_y(rand_engine);
            position.z() = dist_z(rand_engine);
            //Check if a particle at the chosen position overlaps with an other
            overlapping = check_overlaps(position, r + delta, particle_positions, radii);
            if (!overlapping) {
                overlapping = check_overlaps(position, r + delta, mirror_positions, radii);
            }

            if (!overlapping) {
                std::array<double, 3> d1 = {0, 0, 0};
                std::array<double, 3> d2 = {0, 0, 0};
                std::array<bool, 3> active_directions = {false, false, false};
                for (unsigned axis = 0; axis != 3; ++axis) {
                    if (periodic_directions[axis]) {
                        d1[axis] = position[axis] - (r + delta) - periodic_boundaries[axis].min;
                        d2[axis] = periodic_boundaries[axis].max - (r + delta) - position[axis];
                        if (d1[axis] < 0 || d2[axis] < 0) {
                            active_directions[axis] = true;
                            Vec3 new_position = position;
                            new_position[axis] -= (periodic_boundaries[axis].max
                                                   - periodic_boundaries[axis].min)*d1[axis]/(std::abs(d1[axis]));
                            particle_mirror_positions.push_back(new_position);

                        }
                    }
                }
                std::array<std::array<unsigned, 2>, 3> idxs{{{0, 1},
                                                                    {0, 2},
                                                                    {1, 2}}};
                for (const auto idx: idxs) {
                    if (active_directions[idx[0]] && active_directions[idx[1]]) {
                        Vec3 new_position = position;
                        new_position[idx[0]] -= (periodic_boundaries[idx[0]].max
                                                 - periodic_boundaries[idx[0]].min)*d1[idx[0]]/(std::abs(d1[idx[0]]));
                        new_position[idx[1]] -= (periodic_boundaries[idx[1]].max
                                                 - periodic_boundaries[idx[1]].min)*d1[idx[0]]/(std::abs(d1[idx[0]]));
                        particle_mirror_positions.push_back(new_position);
                    }
                }
                if (active_directions[0] && active_directions[1] && active_directions[2]) {
                    Vec3 new_position = position;
                    for (unsigned i = 0; i != 3; ++i) {
                        new_position[i] -= (periodic_boundaries[i].max
                                            - periodic_boundaries[i].min)*d1[i]/(std::abs(d1[i]));
                    }
                    particle_mirror_positions.push_back(new_position);
                }
            }
            for (const auto& mirror_p: particle_mirror_positions){
                overlapping = overlapping || check_overlaps(mirror_p, r + delta, mirror_positions, radii);
                if (!overlapping) {
                    overlapping = check_overlaps(mirror_p, r + delta, particle_positions, radii);
                }
            }
        }
        particle_positions.push_back(position);
        mirror_positions.insert(mirror_positions.end(), particle_mirror_positions.begin(),
                                particle_mirror_positions.end());
    }
    return particle_positions;
}



