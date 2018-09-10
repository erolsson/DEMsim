//
// Created by erolsson on 2018-09-10.
//

#include "utilities.h"

std::vector<DEM::Vec3> DEM::randomF_fill_cylinder(double z0, double z1, double R, std::vector<double> radii)
{
    std::vector<Vec3> particle_positions;
    for(auto particle_radius : radii){
        bool overlapping = true;
        double x;
        double y;
        double z;
        Vec3 position {};
        while(overlapping){
            double r = double(rand())/(double(RAND_MAX)+1.0);
            position.x = -R + particle_radius + 2*r*(R-particle_radius);
            r = double(rand())/(double(RAND_MAX)+1.0);
            position.y = -R + particle_radius + 2*r*(R-particle_radius);
            r = double(rand())/(double(RAND_MAX)+1.0);
            position.z = z0 + particle_radius + r*(z1-2*particle_radius);
            //Check if a particle at the chosen position overlaps with an other
            if(position.x*position.x+position.y*position.y < (R-particle_radius)*(R-particle_radius)) {
                overlapping = check_overlaps(position, particle_radius, particle_positions, radii);
            }
        }

        particle_positions.push_back(position);

    }
    return particle_positions;
}

bool DEM::check_overlaps(const DEM::Vec3& point, double radius, const std::vector<DEM::Vec3>& pos,
                        const std::vector<double>& radii)
{
    return false;
}
