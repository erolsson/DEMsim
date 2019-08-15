//
// Created by erolsson on 2018-09-07.
//

#include "spherical_particle.h"

#include <sstream>

#include "../utilities/vec3.h"

template<typename ForceModel>
DEM::SphericalParticle<ForceModel>::SphericalParticle(double radius, const DEM::Vec3& position,
                                                      const DEM::Vec3& velocity, MaterialBase* material, unsigned id):
        SphericalParticleBase<ForceModel>(radius, position, velocity, material, id)

{
    // Empty constructor
}





